#!/bin/bash
#SBATCH --job-name=rnaseq_autopilot
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
# launch from your project folder:
# sbatch --partition=serial --time=24:00:00 align_fastp_star.sh config.sh

source $1
FASTQ_DIR=${FASTQ_DIR%/}

module load $FASTP_MODULE $STAR_MODULE $SAMTOOLS_MODULE $MULTIQC_MODULE

mkdir -p $OUT/trimmed $OUT/fastp_reports $OUT/star
SHEET=$OUT/sample_sheet.tsv

# 1. make the sample sheet (skipped if one already exists, so you can edit names)
# each sample in its own folder -> folder name is the sample name
# all fastqs in one folder -> sample name comes from the file name
if [ ! -f $SHEET ]; then
  echo -e "sample\tR1\tR2" > $SHEET
  for r1 in $(find $FASTQ_DIR -name "*_R1*.fastq.gz" | sort); do
    f=$(basename $r1)
    r2=$(dirname $r1)/${f/_R1/_R2}
    rel=${r1#$FASTQ_DIR/}
    if [[ $rel == */* ]]; then
      sample=${rel%%/*}
    else
      sample=$(echo ${f%%_R1*} | sed 's/_S[0-9]*$//')
    fi
    echo -e "$sample\t$r1\t$r2" >> $SHEET
  done
fi

# 2. trim and align each sample
for sample in $(awk -F'\t' 'NR>1 {print $1}' $SHEET); do

  # skip finished samples so the job can be resubmitted
  if [ -f $OUT/star/$sample/${sample}_Aligned.sortedByCoord.out.bam.bai ]; then
    continue
  fi

  echo "$(date) $sample"

  R1=$(awk -F'\t' -v x=$sample '$1==x {print $2}' $SHEET)
  R2=$(awk -F'\t' -v x=$sample '$1==x {print $3}' $SHEET)

  fastp -i $R1 -I $R2 \
    -o $OUT/trimmed/${sample}_R1.fastq.gz \
    -O $OUT/trimmed/${sample}_R2.fastq.gz \
    --detect_adapter_for_pe \
    -w $THREADS \
    -j $OUT/fastp_reports/${sample}.json \
    -h $OUT/fastp_reports/${sample}.html

  mkdir -p $OUT/star/$sample
  STAR --runThreadN $THREADS \
    --genomeDir $INDEX \
    --readFilesIn $OUT/trimmed/${sample}_R1.fastq.gz $OUT/trimmed/${sample}_R2.fastq.gz \
    --readFilesCommand zcat \
    --twopassMode $TWOPASS_MODE \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 30000000000 \
    --quantMode GeneCounts \
    --outFileNamePrefix $OUT/star/$sample/${sample}_

  samtools index $OUT/star/$sample/${sample}_Aligned.sortedByCoord.out.bam

  rm $OUT/trimmed/${sample}_R1.fastq.gz $OUT/trimmed/${sample}_R2.fastq.gz
done

# 3. combine gene counts into one matrix
# strand column is picked from the first sample: col 3 or 4 if one is 3x bigger, else col 2 (unstranded)
samples=$(awk -F'\t' 'NR>1 {print $1}' $SHEET)
first=$(echo $samples | cut -d' ' -f1)
s3=$(awk 'NR>4 {a+=$3} END {print a}' $OUT/star/$first/${first}_ReadsPerGene.out.tab)
s4=$(awk 'NR>4 {a+=$4} END {print a}' $OUT/star/$first/${first}_ReadsPerGene.out.tab)
if [ $s3 -gt $((s4 * 3)) ]; then COL=3; elif [ $s4 -gt $((s3 * 3)) ]; then COL=4; else COL=2; fi
echo "using ReadsPerGene column $COL (col3 total $s3, col4 total $s4)"

files=""
for s in $samples; do files="$files $OUT/star/$s/${s}_ReadsPerGene.out.tab"; done
(
  echo -e "gene_id\t$(echo $samples | tr ' ' '\t')"
  awk -v c=$COL 'FNR>4 {g[FNR]=$1; v[FNR]=v[FNR] "\t" $c} END {for (i=5; i<=FNR; i++) print g[i] v[i]}' $files
) > $OUT/gene_counts.tsv

# 4. QC summary across all samples
multiqc $OUT/fastp_reports $OUT/star -o $OUT/multiqc -f
