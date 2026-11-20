# ---- edit these two ----
FASTQ_DIR=/path/to/your/fastqs     # folder with your fastq.gz files (subfolders are fine)
OUT=/path/to/output_folder         # where all results go

# ---- usually leave these alone ----
INDEX=/s/acbm/Gabe/reference/star_index_GRCh38_gencode48
THREADS=8            # keep matching --cpus-per-task in align_fastp_star.sh
TWOPASS_MODE=None    # None = standard DE, Basic = novel junction/transcript work (slower)

FASTP_MODULE=fastp/1.0.1
STAR_MODULE=star/2.7.11b
SAMTOOLS_MODULE=samtools/1.23.1
MULTIQC_MODULE=multiqc/1.35
