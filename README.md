**RNAutopilot**
================

Config-driven bulk RNA-seq pipeline for SLURM.
fastq -> fastp (trimming) -> STAR (alignment + gene counts) -> MultiQC

Each sample is processed one at a time in a single job. Finished samples are
skipped, so if a job hits its time limit, just resubmit the same command.


FILES
-----
config.sh              paths, settings and tool versions (edit this)
align_fastp_star.sh    the pipeline


HOW TO RUN
----------
1. Download this repo into your project folder
   (Code -> Download ZIP, or git clone).

2. Make a sample sheet (tab-separated) with at least these columns:
     sample    sample name, used for all output names
     R1        full path to the R1 fastq.gz
     R2        full path to the R2 fastq.gz
   Any other columns (donor, condition, batch, ...) are ignored by the
   pipeline but are useful to keep for DESeq2.

   Example:
   sample    condition    R1                          R2
   S1_ctrl   ctrl         /path/S1_R1_001.fastq.gz    /path/S1_R2_001.fastq.gz

3. Edit config.sh:
     OUT            where outputs go
     SAMPLE_SHEET   path to your sample sheet
     INDEX          STAR index (built with sjdbOverhang ~ read length - 1)
     THREADS        should match --cpus-per-task in align_fastp_star.sh
     TWOPASS_MODE   Basic for novel junction/transcript work, None for
                    standard differential expression

4. Launch (set partition and time for your cluster):
   sbatch --partition=serial --time=24:00:00 align_fastp_star.sh config.sh


OUTPUTS (in OUT)
----------------
fastp_reports/<sample>.html/.json            trimming and read QC
star/<sample>/<sample>_Aligned.sortedByCoord.out.bam (+ .bai)
star/<sample>/<sample>_ReadsPerGene.out.tab  gene counts
star/<sample>/<sample>_Log.final.out         mapping stats
multiqc/multiqc_report.html                  QC summary across all samples

Trimmed fastqs are deleted after each sample aligns.


STRANDEDNESS
------------
ReadsPerGene.out.tab has 3 count columns after the gene ID:
  column 2 = unstranded
  column 3 = read 1 matches gene strand
  column 4 = read 2 matches gene strand (most stranded kits, e.g. TruSeq stranded)
Check which of columns 3/4 is much larger on a test sample and use that one
for DESeq2.
