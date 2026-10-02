#!/usr/bin/env Rscript
# T1: TCGA COAD+READ primary tumour vs solid-tissue normal (edgeR QL,
# ~ tissue + project), plus two descriptive sensitivity analyses:
# paired tumour/normal (~ patient + tissue) and TCGA tumour vs GTEx
# transverse colon (~ tissue; source fully confounded, direction only).
#
# Usage: Rscript 03_tcga_tumor_vs_normal.R <recount3_processed_dir> <out_dir> <min_log2fc> <max_fdr>

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: 03_tcga_tumor_vs_normal.R <recount3_dir> <out_dir> <min_log2fc> <max_fdr>")
in_dir <- args[[1]]; out_dir <- args[[2]]
min_lfc <- as.numeric(args[[3]]); max_fdr <- as.numeric(args[[4]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rse <- readRDS(file.path(in_dir, "recount3_TCGA_COAD_READ_GTEx_COLON_gene_counts.rds"))
counts <- assay(rse, "counts")
genes <- data.frame(
  gene_id = sub("\\..*$", "", rownames(rse)),
  symbol = as.character(rowData(rse)$gene_name),
  gene_type = as.character(rowData(rse)$gene_type),
  stringsAsFactors = FALSE
)

read_cd <- function(f) read.csv(file.path(in_dir, f), row.names = 1, check.names = FALSE)
tcga <- rbind(read_cd("coldata_tcga_COAD.csv")[, c("study", "tcga.tcga_barcode",
                                                   "tcga.cgc_sample_sample_type",
                                                   "tcga.gdc_cases.samples.is_ffpe",
                                                   "tcga.gdc_cases.submitter_id")],
              read_cd("coldata_tcga_READ.csv")[, c("study", "tcga.tcga_barcode",
                                                   "tcga.cgc_sample_sample_type",
                                                   "tcga.gdc_cases.samples.is_ffpe",
                                                   "tcga.gdc_cases.submitter_id")])
colnames(tcga) <- c("project", "barcode", "sample_type", "is_ffpe", "patient")
tcga$sample_id <- rownames(tcga)
tcga$tissue <- c("Primary Tumor" = "Tumor", "Solid Tissue Normal" = "Normal")[tcga$sample_type]
tcga$lib_size <- colSums(counts[, tcga$sample_id])
tcga$exclusion <- ifelse(is.na(tcga$tissue), "sample_type_not_01_or_11",
                         ifelse(tolower(as.character(tcga$is_ffpe)) %in% c("true", "1"), "FFPE", ""))
# One aliquot per patient per tissue: keep the largest library.
elig <- tcga[tcga$exclusion == "", ]
elig <- elig[order(elig$patient, elig$tissue, -elig$lib_size), ]
dup <- duplicated(elig[, c("patient", "tissue")])
tcga$exclusion[match(elig$sample_id[dup], tcga$sample_id)] <- "duplicate_aliquot_smaller_library"
tcga$included <- tcga$exclusion == ""
write.csv(tcga, file.path(out_dir, "TCGA_sample_inclusion.csv"), row.names = FALSE)
inc <- tcga[tcga$included, ]
cat("Included TCGA samples:\n"); print(table(inc$tissue, inc$project))

fit_qlf <- function(y_counts, design, coef, group) {
  y <- DGEList(y_counts)
  keep <- filterByExpr(y, group = group)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
  res$row <- rownames(res)
  res
}

# ---- T1 primary --------------------------------------------------------------
inc$tissue <- factor(inc$tissue, levels = c("Normal", "Tumor"))
inc$project <- factor(inc$project)
d1 <- model.matrix(~ tissue + project, data = inc)
t1 <- fit_qlf(counts[, inc$sample_id], d1, "tissueTumor", inc$tissue)
out <- data.frame(row = rownames(counts), genes)
out <- merge(out, data.frame(row = t1$row, T1_log2FC = t1$logFC, T1_logCPM = t1$logCPM,
                             T1_PValue = t1$PValue, T1_FDR = t1$FDR), by = "row", all.x = TRUE)

# ---- Paired sensitivity ------------------------------------------------------
paired_pat <- names(which(table(inc$patient) == 2))
pi <- inc[inc$patient %in% paired_pat, ]
pi$patient <- factor(pi$patient)
d2 <- model.matrix(~ patient + tissue, data = pi)
t2 <- fit_qlf(counts[, pi$sample_id], d2, "tissueTumor", pi$tissue)
out <- merge(out, data.frame(row = t2$row, Paired_log2FC = t2$logFC, Paired_FDR = t2$FDR),
             by = "row", all.x = TRUE)

# ---- GTEx sensitivity --------------------------------------------------------
gtex <- read_cd("coldata_gtex_COLON.csv")
gtex_ids <- rownames(gtex)[gtex[["gtex.smtsd"]] == "Colon - Transverse"]
tum <- inc$sample_id[inc$tissue == "Tumor"]
gs <- data.frame(sample_id = c(gtex_ids, tum),
                 tissue = factor(rep(c("Normal", "Tumor"), c(length(gtex_ids), length(tum))),
                                 levels = c("Normal", "Tumor")))
d3 <- model.matrix(~ tissue, data = gs)
t3 <- fit_qlf(counts[, gs$sample_id], d3, "tissueTumor", gs$tissue)
out <- merge(out, data.frame(row = t3$row, GTEx_log2FC = t3$logFC, GTEx_FDR = t3$FDR),
             by = "row", all.x = TRUE)

out$T1_tested <- !is.na(out$T1_log2FC)
out$T1_pass <- out$T1_tested & out$T1_log2FC >= min_lfc & out$T1_FDR < max_fdr
out$Paired_concordant <- !is.na(out$Paired_log2FC) & sign(out$Paired_log2FC) == sign(out$T1_log2FC)
out$GTEx_concordant <- !is.na(out$GTEx_log2FC) & sign(out$GTEx_log2FC) == sign(out$T1_log2FC)
out <- out[order(out$T1_PValue), ]
out$row <- NULL
write.csv(out, file.path(out_dir, "TCGA_tumor_vs_normal_DEG.csv"), row.names = FALSE)

summ <- data.frame(
  contrast = c("T1_tumor_vs_normal", "Paired_tumor_vs_normal", "TCGA_tumor_vs_GTEx_transverse"),
  n_normal = c(sum(inc$tissue == "Normal"), length(paired_pat), length(gtex_ids)),
  n_tumor = c(sum(inc$tissue == "Tumor"), length(paired_pat), length(tum)),
  n_tested = c(nrow(t1), nrow(t2), nrow(t3)),
  n_up = c(sum(out$T1_pass, na.rm = TRUE),
           sum(out$Paired_log2FC >= min_lfc & out$Paired_FDR < max_fdr, na.rm = TRUE),
           sum(out$GTEx_log2FC >= min_lfc & out$GTEx_FDR < max_fdr, na.rm = TRUE))
)
write.csv(summ, file.path(out_dir, "TCGA_contrast_summary.csv"), row.names = FALSE)
print(summ)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_03_tcga.txt"))
