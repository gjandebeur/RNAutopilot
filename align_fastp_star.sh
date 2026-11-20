#!/bin/bash
#SBATCH --job-name=rnaseq_autopilot
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
# set partition and time when you launch, e.g.
# sbatch --partition=serial --time=24:00:00 align_fastp_star.sh config.sh

source $1

module load $FASTP_MODULE $STAR_MODULE $SAMTOOLS_MODULE $MULTIQC_MODULE

mkdir -p $OUT/trimmed $OUT/fastp_reports $OUT/star

# find the sample, R1 and R2 columns in the sample sheet
col() { head -1 $SAMPLE_SHEET | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
S=$(col sample)
C1=$(col R1)
C2=$(col R2)

for sample in $(awk -F'\t' -v s=$S 'NR>1 {print $s}' $SAMPLE_SHEET); do

  # skip finished samples so the job can be resubmitted
  if [ -f $OUT/star/$sample/${sample}_Aligned.sortedByCoord.out.bam.bai ]; then
    continue
  fi

  echo "$(date) $sample"

  R1=$(awk -F'\t' -v s=$S -v c=$C1 -v x=$sample 'NR>1 && $s==x {print $c}' $SAMPLE_SHEET)
  R2=$(awk -F'\t' -v s=$S -v c=$C2 -v x=$sample 'NR>1 && $s==x {print $c}' $SAMPLE_SHEET)

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

multiqc $OUT/fastp_reports $OUT/star -o $OUT/multiqc
