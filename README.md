**RNAutopilot**
================


Bulk RNA-seq on SLURM, start to finish:
fastq -> fastp (trimming) -> STAR (alignment + gene counts) -> count matrix -> MultiQC

You edit two lines, launch one command, and wait.


FILES
-----
config.sh              settings (edit the top two lines)
align_fastp_star.sh    the pipeline


HOW TO RUN
----------
1. Download this repo into your project folder
   (Code -> Download ZIP, or git clone).

2. Open config.sh and set:
     FASTQ_DIR   the folder with your fastq.gz files
     OUT         where you want the results

3. Launch:
   sbatch --partition=serial --time=24:00:00 align_fastp_star.sh config.sh

That's it. If the job runs out of time, run the same command again.
Finished samples are skipped and it picks up where it stopped.


YOUR FASTQ FILES
----------------
Files need _R1 / _R2 in their names (e.g. S1_S5_R1_001.fastq.gz).
Sample names are made automatically:
  - each sample in its own folder  -> the folder name is used
  - all fastqs in one folder       -> the name comes from the file name
One R1/R2 pair per sample. If a sample was split across lanes, merge those
files first.


WHAT YOU GET (in OUT)
---------------------
gene_counts.tsv                  counts for every gene x sample, ready for DESeq2
sample_sheet.tsv                 the sample list the pipeline made
multiqc/multiqc_report.html      QC summary for all samples (open in a browser)
fastp_reports/                   trimming QC per sample
star/<sample>/                   BAM, BAM index, gene counts, mapping stats

Trimmed fastqs are deleted after each sample aligns to save space.

