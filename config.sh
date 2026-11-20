# paths
OUT=/path/to/project_folder
SAMPLE_SHEET=/path/to/project_folder/sample_sheet.tsv   # needs columns: sample, R1, R2
INDEX=/s/acbm/Gabe/reference/star_index_GRCh38_gencode48

# run settings
THREADS=8
TWOPASS_MODE=Basic   # Basic = better novel junctions (slower), None = standard DE

# tool versions
FASTP_MODULE=fastp/1.0.1
STAR_MODULE=star/2.7.11b
SAMTOOLS_MODULE=samtools/1.23.1
MULTIQC_MODULE=multiqc/1.35
