#!/usr/bin/env Rscript
# Panel-aware edgeR QL contrast (CIOC data-QC amendment).
# Identical to the original gate contrast except that predefined panel genes
# with any counts are exempt from filterByExpr:
#   keep = filterByExpr(...) | (panel & rowSums(counts) > 0)
# Two modes reproduce the two original implementations:
#   generic       step3_cancer/04_pseudobulk_edger.R: filterByExpr(group), normLibSizes, glmQLFit
#   design_robust step2 HGCA / mouse scripts: filterByExpr(design), calcNormFactors, glmQLFit(robust=TRUE)
# Usage: Rscript panel_edger.R <counts.csv.gz> <samples.csv> <case> <ref> <covariate|none> <out.csv>
#                              <min_units> <panel_symbols.txt> <generic|design_robust>
# counts: gene_id, symbol, <pseudobulk columns>; samples: pseudobulk, group, eligible[, covariate]
suppressPackageStartupMessages(library(edgeR))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 9) stop("Usage: panel_edger.R counts samples case ref covariate out min_units panel mode")
cnt <- read.csv(a[1], check.names = FALSE)
smp <- read.csv(a[2], stringsAsFactors = FALSE)
case <- a[3]; ref <- a[4]; covar <- a[5]; out_f <- a[6]; min_units <- as.integer(a[7])
panel <- readLines(a[8]); mode <- a[9]
smp <- smp[smp$eligible %in% c(TRUE, "True", "TRUE") & smp$group %in% c(case, ref), ]
smp$group <- factor(smp$group, levels = c(ref, case))
n_per <- table(smp$group)
if (any(n_per < min_units)) stop("Hard stop: too few eligible units")
mat <- as.matrix(cnt[, smp$pseudobulk])
rownames(mat) <- make.unique(as.character(cnt$gene_id))
form <- ~ group
if (covar != "none") {
  smp[[covar]] <- factor(smp[[covar]])
  if (nlevels(smp[[covar]]) > 1) {
    trial <- model.matrix(~ group + smp[[covar]], data = smp)
    if (qr(trial)$rank == ncol(trial)) form <- as.formula(paste("~", covar, "+ group"))
  }
}
design <- model.matrix(form, data = smp)
y <- DGEList(mat)
is_panel <- cnt$symbol %in% panel
base_keep <- if (mode == "generic") filterByExpr(y, group = smp$group) else filterByExpr(y, design = design)
forced <- is_panel & rowSums(mat) > 0 & !base_keep
keep <- base_keep | forced
if (mode == "generic") {
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
} else {
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
}
coef <- paste0("group", case)
res <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
cpm_all <- cpm(y, normalized.lib.sizes = TRUE)
idx <- match(rownames(res), rownames(mat))
out <- data.frame(gene_id = sub("\\..*$", "", cnt$gene_id[idx]), symbol = cnt$symbol[idx],
                  log2FC = res$logFC, logCPM = res$logCPM, PValue = res$PValue,
                  FDR = p.adjust(res$PValue, "BH"),
                  panel_gene = is_panel[idx], forced_by_panel_exemption = forced[idx],
                  mean_CPM_case = rowMeans(cpm_all[rownames(res), smp$group == case, drop = FALSE]),
                  mean_CPM_ref = rowMeans(cpm_all[rownames(res), smp$group == ref, drop = FALSE]))
# panel genes with zero counts in all contrast samples (not testable; zero_expression)
zero <- cnt$symbol %in% panel & rowSums(mat) == 0
if (any(zero)) {
  out <- rbind(out, data.frame(gene_id = sub("\\..*$", "", cnt$gene_id[zero]), symbol = cnt$symbol[zero],
                               log2FC = NA, logCPM = NA, PValue = NA, FDR = NA, panel_gene = TRUE,
                               forced_by_panel_exemption = FALSE, mean_CPM_case = 0, mean_CPM_ref = 0))
}
write.csv(out, out_f, row.names = FALSE)
message(sprintf("%s: %s vs %s (%d vs %d), design %s, tested %d, panel forced %d, panel zero %d",
                basename(out_f), case, ref, n_per[[case]], n_per[[ref]], deparse(form), sum(keep), sum(forced), sum(zero)))
