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
message(sprintf("Wrote %d tested genes; %d pass M1.", nrow(out), sum(out$Mouse_pass)))
