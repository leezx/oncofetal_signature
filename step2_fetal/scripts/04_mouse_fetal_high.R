#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 04_mouse_fetal_high.R <counts.tsv.gz> <output_dir>")
}

input_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("ENTREZID", "SYMBOL")
if (!all(required %in% names(dat))) stop("Missing ENTREZID or SYMBOL column")

nonempty_symbols <- dat$SYMBOL[!is.na(dat$SYMBOL) & dat$SYMBOL != ""]
n_duplicated_symbol_rows <- sum(duplicated(nonempty_symbols))
symbol_check <- data.frame(
  metric = c("total_rows", "nonempty_symbol_rows", "duplicated_nonempty_symbol_rows"),
  value = c(nrow(dat), length(nonempty_symbols), n_duplicated_symbol_rows)
)
write.csv(symbol_check, file.path(output_dir, "Mouse_symbol_duplicate_check.csv"), row.names = FALSE)
if (n_duplicated_symbol_rows > 0L) {
  stop("Duplicated non-empty mouse gene symbols detected; resolve by stable gene ID before analysis")
}

sample_columns <- setdiff(names(dat), required)
stage <- ifelse(grepl("^E16", sample_columns), "fetal",
                ifelse(grepl("^Adult", sample_columns), "adult", NA_character_))
if (anyNA(stage)) stop("Could not assign every sample to fetal or adult")
stage <- factor(stage, levels = c("adult", "fetal"))

counts <- as.matrix(dat[, sample_columns, drop = FALSE])
storage.mode(counts) <- "integer"
gene_label <- ifelse(is.na(dat$SYMBOL) | dat$SYMBOL == "", as.character(dat$ENTREZID), dat$SYMBOL)
rownames(counts) <- make.unique(gene_label)

y <- DGEList(counts = counts, genes = dat[, required, drop = FALSE], group = stage)
design <- model.matrix(~ stage)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "stagefetal")
result <- topTags(test, n = Inf, sort.by = "none")$table
result$FDR <- p.adjust(result$PValue, method = "BH")
result$Mouse_pass <- result$logFC >= 0.5 & result$FDR < 0.05

raw_cpm <- cpm(y$counts, normalized.lib.sizes = FALSE, log = FALSE)
extreme <- result[result$logFC > 5, , drop = FALSE]
extreme_cpm <- raw_cpm[rownames(extreme), , drop = FALSE]
colnames(extreme_cpm) <- c(
  "E16.5_1_raw_CPM", "E16.5_2_raw_CPM", "E16.5_3_raw_CPM",
  "Adult_1_raw_CPM", "Adult_2_raw_CPM", "Adult_3_raw_CPM"
)
extreme_out <- data.frame(
  mouse_entrez_id = extreme$ENTREZID,
  mouse_gene = extreme$SYMBOL,
  Mouse_log2FC = extreme$logFC,
  Mouse_logCPM = extreme$logCPM,
  Mouse_PValue = extreme$PValue,
  Mouse_FDR = extreme$FDR,
  extreme_cpm,
  fetal_mean_raw_CPM = rowMeans(extreme_cpm[, 1:3, drop = FALSE]),
  adult_mean_raw_CPM = rowMeans(extreme_cpm[, 4:6, drop = FALSE]),
  low_fetal_count_flag = rowMeans(extreme_cpm[, 1:3, drop = FALSE]) < 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
extreme_out <- extreme_out[order(-extreme_out$Mouse_log2FC), ]
write.csv(extreme_out, file.path(output_dir, "Mouse_log2FC_gt5_raw_CPM.csv"), row.names = FALSE)

result <- result[order(!result$Mouse_pass, -result$logFC), ]

out <- data.frame(
  mouse_entrez_id = result$ENTREZID,
  mouse_gene = result$SYMBOL,
  Mouse_log2FC = result$logFC,
  Mouse_logCPM = result$logCPM,
  Mouse_F = result$F,
  Mouse_PValue = result$PValue,
  Mouse_FDR = result$FDR,
  Mouse_pass = result$Mouse_pass,
  stringsAsFactors = FALSE
)
write.csv(out, file.path(output_dir, "Mouse_DE.csv"), row.names = FALSE)

sample_manifest <- data.frame(sample_id = sample_columns, stage = as.character(stage))
write.csv(sample_manifest, file.path(output_dir, "Mouse_sample_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Mouse_sessionInfo.txt"))
message(sprintf(
  "Wrote %d tested genes; %d pass M1; %d genes have log2FC > 5; duplicated symbols = %d.",
  nrow(out), sum(out$Mouse_pass), nrow(extreme_out), n_duplicated_symbol_rows
))
