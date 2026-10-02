#!/usr/bin/env Rscript
# H-bulk1 / H-bulk2 (benchmark v2 addendum v1.1): primary human fetal small
# intestine (Roadmap, SRP001371) vs adult small intestine / duodenum (HPA,
# ERP003613), recount3 raw gene counts, edgeR QL ~ stage.
#
# Sample selection uses metadata only:
#   fetal = the six Roadmap GSMs used by Finkbeiner 2015 (fetal SI total RNA)
#   adult = HPA samples whose tissue is duodenum or small intestine
# Usage: Rscript 07_bulk_primary_tissue.R <dataset_dir> <work_dir> <tables_dir> <edger_script>

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: 07_bulk_primary_tissue.R <dataset_dir> <work_dir> <tables_dir> <edger_script>")
ds <- args[[1]]; work <- args[[2]]; tables <- args[[3]]; edger <- args[[4]]
for (d in c(file.path(ds, "raw", "recount3_cache"), file.path(ds, "processed", "v0.1"), work, tables))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(recount3_url = "https://duffel.rail.bio/recount3")
bfc <- BiocFileCache::BiocFileCache(file.path(ds, "raw", "recount3_cache"), ask = FALSE)
proj <- available_projects(bfc = bfc)

fetch <- function(id) {
  rse <- create_rse(subset(proj, project == id & file_source == "sra"), type = "gene",
                    annotation = "gencode_v26", bfc = bfc)
  assay(rse, "counts") <- transform_counts(rse)
  assays(rse) <- assays(rse)["counts"]
  rse
}
roadmap <- fetch("SRP001371")
hpa <- fetch("ERP003613")

# Fetal: the six Roadmap fetal SI total-RNA samples used by Finkbeiner 2015
# (six distinct donors), selected by SRA experiment accession.
fetal_map <- data.frame(
  gsm = c("GSM1059486", "GSM1059507", "GSM1059508", "GSM1059517", "GSM1059519", "GSM1059521"),
  srx = c("SRX213997", "SRX214018", "SRX214019", "SRX214028", "SRX214030", "SRX214032"),
  donor = c("H-23769", "H-23887", "H-23914", "H-23808", "H-23941", "H-23964"),
  age_days = c(108, 108, 91, 115, 120, 98))
fi <- match(fetal_map$srx, roadmap$sra.experiment_acc)
if (anyNA(fi)) stop("Fetal SRX missing from recount3 SRP001371: ", paste(fetal_map$srx[is.na(fi)], collapse = ","))
fetal_counts <- assay(roadmap, "counts")[, fi]
colnames(fetal_counts) <- paste0("Fetal_", fetal_map$donor)

# Adult: HPA E-MTAB-1733 duodenum and small-intestine samples; each sample was
# sequenced as two runs, which are summed (technical replicates).
alias <- sub("^E-MTAB-1733:", "", as.character(hpa$sra.sample_title))
keep <- grepl("^(duodenum|smallintestine)_", alias)
adult_counts <- sapply(split(which(keep), alias[keep]), function(i) rowSums(assay(hpa, "counts")[, i, drop = FALSE]))
runs_per_sample <- table(alias[keep])
colnames(adult_counts) <- paste0("Adult_", colnames(adult_counts))

common <- intersect(rownames(fetal_counts), rownames(adult_counts))
counts <- cbind(fetal_counts[common, ], adult_counts[common, ])
adult_tissue <- ifelse(grepl("duodenum", colnames(adult_counts)), "adult duodenum", "adult small intestine")
units <- data.frame(
  pseudobulk = colnames(counts),
  patient = colnames(counts),
  group = rep(c("Fetal", "Adult"), c(ncol(fetal_counts), ncol(adult_counts))),
  source = rep(c("Roadmap SRP001371", "HPA ERP003613 (E-MTAB-1733)"), c(ncol(fetal_counts), ncol(adult_counts))),
  sample_id = c(fetal_map$gsm, names(runs_per_sample)),
  runs_summed = c(rep(1L, nrow(fetal_map)), as.integer(runs_per_sample)),
  age = c(paste0("day ", fetal_map$age_days), rep("adult", ncol(adult_counts))),
  tissue = c(rep("fetal small intestine", ncol(fetal_counts)), adult_tissue),
  library_size = colSums(counts),
  n_cells = 1L, eligible = TRUE
)
cnt <- data.frame(gene_id = sub("\\..*$", "", common),
                  symbol = as.character(rowData(roadmap)[common, "gene_name"]), counts, check.names = FALSE)
cnt <- cnt[!duplicated(cnt$gene_id), ]
f_counts <- file.path(work, "Hbulk_primary_tissue_counts.csv.gz")
write.csv(cnt, gzfile(f_counts), row.names = FALSE)
write.csv(units, file.path(work, "Hbulk1_units.csv"), row.names = FALSE)
u2 <- units; u2$eligible <- u2$group == "Fetal" | u2$tissue == "adult duodenum"
write.csv(u2, file.path(work, "Hbulk2_units.csv"), row.names = FALSE)
write.csv(units, file.path(tables, "Hbulk_primary_tissue_units.csv"), row.names = FALSE)
print(units[, c("pseudobulk", "group", "tissue", "sample_id", "runs_summed", "age", "library_size")])

run <- function(u, out, min_n) {
  status <- system2("Rscript", c(edger, f_counts, u, "Fetal", "Adult", "none", out, min_n))
  if (status != 0) stop("edgeR failed for ", out)
}
run(file.path(work, "Hbulk1_units.csv"), file.path(work, "Hbulk1_Roadmap_vs_HPA_SI.csv"), 3)
run(file.path(work, "Hbulk2_units.csv"), file.path(work, "Hbulk2_Roadmap_vs_HPA_duodenum.csv"), 2)
writeLines(capture.output(sessionInfo()), file.path(work, "sessionInfo_07_bulk.txt"))
