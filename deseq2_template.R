# DESeq2 differential expression, run locally in RStudio
#
# before running:
#  1. copy gene_counts.tsv and sample_sheet.tsv from the pipeline OUT folder
#     into the same folder as this script
#  2. open sample_sheet.tsv in Excel, add a column called "condition" and fill
#     in each sample's group (e.g. ctrl / treated). Samples left blank are skipped.
#  3. in RStudio: Session -> Set Working Directory -> To Source File Location
#
# first time only: install.packages("BiocManager"); BiocManager::install("DESeq2")

library(DESeq2)

# ---- edit these ----
counts_file  <- "gene_counts.tsv"
sheet_file   <- "sample_sheet.tsv"
reference    <- "ctrl"        # baseline group everything is compared to
batch_column <- NULL          # e.g. "donor" for paired samples, otherwise NULL
out_dir      <- "deseq2_results"
# --------------------

dir.create(out_dir, showWarnings = FALSE)

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
sheet  <- read.delim(sheet_file)
sheet  <- sheet[!is.na(sheet$condition) & sheet$condition != "", ]
rownames(sheet) <- sheet$sample
counts <- counts[, sheet$sample]

sheet$condition <- relevel(factor(sheet$condition), ref = reference)

if (is.null(batch_column)) {
  design <- ~ condition
} else {
  sheet[[batch_column]] <- factor(sheet[[batch_column]])
  design <- as.formula(paste("~", batch_column, "+ condition"))
}

dds <- DESeqDataSetFromMatrix(counts, sheet, design)

# keep genes with at least 10 reads in as many samples as the smallest group
smallest_group <- min(table(sheet$condition))
dds <- dds[rowSums(counts(dds) >= 10) >= smallest_group, ]

dds <- DESeq(dds)

# PCA of all samples
vsd <- vst(dds)
pdf(file.path(out_dir, "PCA.pdf"))
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

# each group vs the reference
for (grp in setdiff(levels(sheet$condition), reference)) {
  name <- paste0(grp, "_vs_", reference)
  res <- results(dds, contrast = c("condition", grp, reference))
  res <- res[order(res$padj), ]
  write.csv(as.data.frame(res), file.path(out_dir, paste0(name, ".csv")))

  pdf(file.path(out_dir, paste0(name, "_MA.pdf")))
  plotMA(res, main = paste(grp, "vs", reference))
  dev.off()

  cat(name, ":", sum(res$padj < 0.05, na.rm = TRUE), "genes at padj < 0.05\n")
}
