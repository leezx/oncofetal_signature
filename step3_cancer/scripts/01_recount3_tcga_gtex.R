#!/usr/bin/env Rscript
# Acquire uniformly processed (recount3 / Monorail) gene-level counts for
# TCGA COAD, TCGA READ, and GTEx COLON, and save one combined
# RangedSummarizedExperiment with raw read counts and sample metadata.
#
# Usage: Rscript 01_recount3_tcga_gtex.R <dataset_dir>
#   <dataset_dir>/raw/recount3_cache   recount3 BiocFileCache (immutable source files)
#   <dataset_dir>/processed/v0.1       combined counts object + sample table

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 01_recount3_tcga_gtex.R <dataset_dir>")
dataset_dir <- args[[1]]
cache_dir <- file.path(dataset_dir, "raw", "recount3_cache")
out_dir <- file.path(dataset_dir, "processed", "v0.1")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(recount3_url = "https://duffel.rail.bio/recount3")
bfc <- BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)

projects <- available_projects(bfc = bfc)
wanted <- subset(
  projects,
  (file_source == "tcga" & project %in% c("COAD", "READ")) |
    (file_source == "gtex" & project == "COLON")
)
stopifnot(nrow(wanted) == 3)

fetch <- function(i) {
  rse <- create_rse(wanted[i, ], type = "gene", annotation = "gencode_v26", bfc = bfc)
  assay(rse, "counts") <- transform_counts(rse)  # coverage -> read counts
  assays(rse) <- assays(rse)["counts"]
  rse
}
rse_list <- lapply(seq_len(nrow(wanted)), fetch)

# Keep only metadata columns shared by all three projects so cbind succeeds;
# source-specific clinical columns are written separately.
for (i in seq_along(rse_list)) {
  src <- wanted$file_source[i]
  cd <- as.data.frame(colData(rse_list[[i]]))
  write.csv(cd, file.path(out_dir, sprintf("coldata_%s_%s.csv", src, wanted$project[i])), row.names = TRUE)
}
common <- Reduce(intersect, lapply(rse_list, function(x) colnames(colData(x))))
rse_list <- lapply(rse_list, function(x) { colData(x) <- colData(x)[, common]; x })
rse <- do.call(cbind, rse_list)

saveRDS(rse, file.path(out_dir, "recount3_TCGA_COAD_READ_GTEx_COLON_gene_counts.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_01_recount3.txt"))
cat(sprintf("Saved %d genes x %d samples\n", nrow(rse), ncol(rse)))
print(table(rse$study))
