#!/usr/bin/env Rscript
# Generic patient-level pseudobulk edgeR QL contrast used for S1, S2 and the
# proliferation control. Only eligible pseudobulks (>= min cells) are used.
#
# Usage:
#   Rscript 04_pseudobulk_edger.R <counts.csv.gz> <samples.csv> <case_group> \
#     <reference_group> <covariate|none> <out.csv>
# counts:  gene_id, symbol, then one column per pseudobulk
# samples: pseudobulk, patient, group, eligible[, covariate]
# Output log2FC is case / reference.

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: 04_pseudobulk_edger.R <counts> <samples> <case> <reference> <covariate|none> <out.csv>")
counts_f <- args[[1]]; samples_f <- args[[2]]; case <- args[[3]]; ref <- args[[4]]
covar <- args[[5]]; out_f <- args[[6]]

cnt <- read.csv(counts_f, check.names = FALSE)
smp <- read.csv(samples_f, stringsAsFactors = FALSE)
smp <- smp[smp$eligible %in% c(TRUE, "True", "TRUE") & smp$group %in% c(case, ref), ]
smp$group <- factor(smp$group, levels = c(ref, case))
n_per <- table(smp$group)
if (any(n_per < 3)) stop(sprintf("Hard stop: fewer than 3 eligible patients in a group (%s)",
                                 paste(names(n_per), n_per, collapse = ", ")))

mat <- as.matrix(cnt[, smp$pseudobulk])
rownames(mat) <- cnt$gene_id

form <- ~ group
if (covar != "none") {
  smp[[covar]] <- factor(smp[[covar]])
  if (nlevels(smp[[covar]]) > 1) {
    trial <- model.matrix(~ group + smp[[covar]], data = smp)
    if (qr(trial)$rank == ncol(trial)) {
      form <- as.formula(paste("~", covar, "+ group"))
    } else {
      message("Covariate ", covar, " not estimable; falling back to ~ group")
    }
  }
}
design <- model.matrix(form, data = smp)
y <- DGEList(mat)
keep <- filterByExpr(y, group = smp$group)
y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
coef <- paste0("group", case)
res <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "PValue")$table

out <- data.frame(gene_id = sub("\\..*$", "", rownames(res)),
                  symbol = cnt$symbol[match(rownames(res), cnt$gene_id)],
                  log2FC = res$logFC, logCPM = res$logCPM, PValue = res$PValue, FDR = res$FDR)
write.csv(out, out_f, row.names = FALSE)
meta <- data.frame(case = case, reference = ref, n_case = n_per[[case]], n_reference = n_per[[ref]],
                   design = deparse(form), n_tested = nrow(out))
write.csv(meta, sub("\\.csv$", "_model.csv", out_f), row.names = FALSE)
print(meta)
