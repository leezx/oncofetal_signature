#!/usr/bin/env Rscript
# Addendum v1.2: stage-resolved HGCA (H1) contrasts, breakpoint 9 PCW (frozen).
# Same donor pseudobulks, filters and edgeR QL model as step2_fetal/02_hgca_fetal_high.R;
# only the fetal arm is subset by donor age.
# Usage: Rscript 09_hgca_stage_resolved.R <hgca_processed_dir> <de_dir> <tables_dir>
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
proc <- args[[1]]; de <- args[[2]]; tab <- args[[3]]
counts <- read.delim(file.path(proc, "HGCA_primary_pseudobulk_counts.tsv.gz"), row.names = 1, check.names = FALSE)
man <- read.csv(file.path(proc, "HGCA_primary_sample_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), man$donor_id))
man$age_weeks <- ifelse(man$stage == "fetal", as.numeric(sub("Wk$", "", man$age)), NA)
man$stage_group <- ifelse(man$stage == "adult", "Adult", ifelse(man$age_weeks < 9, "Fetal_early", "Fetal_midlate"))
write.csv(man[, c("donor_id", "stage", "age", "age_weeks", "stage_group", "n_cells")],
          file.path(tab, "HGCA_stage_resolved_units.csv"), row.names = FALSE)
print(table(man$stage_group))

run <- function(case, ref, name) {
  sel <- man$stage_group %in% c(case, ref)
  grp <- factor(ifelse(man$stage_group[sel] %in% case, "case", "ref"), levels = c("ref", "case"))
  design <- model.matrix(~ grp)
  y <- DGEList(counts = counts[, sel], group = grp)
  y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- topTags(glmQLFTest(fit, coef = "grpcase"), n = Inf, sort.by = "none")$table
  out <- data.frame(symbol = rownames(res), log2FC = res$logFC, logCPM = res$logCPM, F = res$F,
                    PValue = res$PValue, FDR = p.adjust(res$PValue, "BH"),
                    n_case = sum(grp == "case"), n_ref = sum(grp == "ref"))
  write.csv(out, file.path(de, paste0(name, ".csv")), row.names = FALSE)
  tn <- out[out$symbol == "TNFRSF12A", ]
  message(sprintf("%s: %d vs %d donors; TNFRSF12A %+.2f P %.3g", name, sum(grp == "case"), sum(grp == "ref"),
                  tn$log2FC, tn$PValue))
}
run(c("Fetal_early", "Fetal_midlate"), "Adult", "HGCA_all_fetal_vs_adult")
run("Fetal_early", "Adult", "HGCA_early_fetal_vs_adult")
run("Fetal_midlate", "Adult", "HGCA_late_fetal_vs_adult")
run("Fetal_midlate", "Fetal_early", "HGCA_late_vs_early_fetal")
writeLines(capture.output(sessionInfo()), file.path(de, "sessionInfo_09_hgca.txt"))
