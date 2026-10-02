#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 02_hgca_fetal_high.R <processed_dir> <output_dir> <min_log2fc>")
}

processed_dir <- args[[1]]
output_dir <- args[[2]]
min_log2fc <- as.numeric(args[[3]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_edger <- function(count_path, manifest_path) {
  counts <- read.delim(count_path, row.names = 1, check.names = FALSE)
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(counts), manifest$donor_id))
  stage <- factor(manifest$stage, levels = c("adult", "fetal"))
  design <- model.matrix(~ stage)
  y <- DGEList(counts = counts, group = stage)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "stagefetal")
  result <- topTags(test, n = Inf, sort.by = "none")$table
  result$FDR <- p.adjust(result$PValue, method = "BH")
  result$gene <- rownames(result)
  result
}

primary <- run_edger(
  file.path(processed_dir, "HGCA_primary_pseudobulk_counts.tsv.gz"),
  file.path(processed_dir, "HGCA_primary_sample_manifest.csv")
)
primary$HGCA_primary_pass <- primary$logFC >= min_log2fc & primary$FDR < 0.05
primary_out <- primary[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "HGCA_primary_pass")]
names(primary_out)[2:6] <- paste0("HGCA_", names(primary_out)[2:6])
names(primary_out)[names(primary_out) == "HGCA_logFC"] <- "HGCA_log2FC"

control <- run_edger(
  file.path(processed_dir, "HGCA_proliferative_pseudobulk_counts.tsv.gz"),
  file.path(processed_dir, "HGCA_proliferative_sample_manifest.csv")
)
control$HGCA_proliferative_direction_pass <- control$logFC > 0
control_out <- control[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "HGCA_proliferative_direction_pass")]
names(control_out)[2:6] <- paste0("HGCA_proliferative_", names(control_out)[2:6])
names(control_out)[names(control_out) == "HGCA_proliferative_logFC"] <- "HGCA_proliferative_log2FC"

merged <- merge(primary_out, control_out, by = "gene", all.x = TRUE, sort = FALSE)
gene_manifest <- read.csv(file.path(processed_dir, "HGCA_gene_manifest.csv"), stringsAsFactors = FALSE)
merged <- merge(gene_manifest[, c("gene", "gene_ids")], merged, by = "gene", all.y = TRUE, sort = FALSE)
names(merged)[names(merged) == "gene_ids"] <- "human_ensembl_gene_id"
merged$H1_pass <- merged$HGCA_primary_pass &
  !is.na(merged$HGCA_proliferative_direction_pass) &
  merged$HGCA_proliferative_direction_pass
merged <- merged[order(!merged$H1_pass, -merged$HGCA_log2FC), ]

write.csv(merged, file.path(output_dir, "HGCA_DE.csv"), row.names = FALSE)
write.csv(control_out, file.path(output_dir, "HGCA_proliferative_control.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "HGCA_sessionInfo.txt"))
message(sprintf(
  "Tested %d primary genes; %d pass primary and %d pass complete H1.",
  nrow(primary_out), sum(primary_out$HGCA_primary_pass), sum(merged$H1_pass)
))
