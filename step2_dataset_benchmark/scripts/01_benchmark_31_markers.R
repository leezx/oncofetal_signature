#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1L) normalizePath(args[[1]]) else normalizePath(".")
data_root <- if (length(args) >= 2L) args[[2]] else
  "/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/dataset_qualification_v0.1"
out_dir <- file.path(repo, "step2_dataset_benchmark", "results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(data_root, "processed"), recursive = TRUE, showWarnings = FALSE)

candidates <- read.delim(file.path(repo, "step2_fetal", "config", "literature_candidates_31.tsv"),
                         stringsAsFactors = FALSE)
# Older human annotations (Senger RPKM, Fawkner 10x features) still use
# pre-2019 symbols (CYR61, CTGF). Look candidates up under their source symbol.
aliases <- read.delim(file.path(repo, "step3_cancer", "config", "symbol_aliases.tsv"), stringsAsFactors = FALSE)
source_symbol <- function(g) ifelse(g %in% aliases$current_hgnc_symbol,
                                    aliases$source_symbol[match(g, aliases$current_hgnc_symbol)], g)
orth <- read.delim(
  "/Volumes/Stelligen_SSD/Stelligen/DATA/1.Databases/Ensembl_orthologues/release_116/raw/human_mouse_orthologues.tsv",
  check.names = FALSE, stringsAsFactors = FALSE
)
orth <- orth[orth[["Mouse homology type"]] == "ortholog_one2one",
             c("Gene name", "Mouse gene name")]
names(orth) <- c("human_gene", "mouse_gene")
gene_info <- read.delim(gzfile(file.path(data_root, "raw", "Mus_musculus.gene_info.gz")),
                        check.names = FALSE, stringsAsFactors = FALSE)
gene_info <- unique(gene_info[c("Symbol", "GeneID")])
mouse_map <- merge(orth, gene_info, by.x = "mouse_gene", by.y = "Symbol", all.x = TRUE)
names(mouse_map)[names(mouse_map) == "GeneID"] <- "mouse_entrez_id"
mouse_map <- mouse_map[!duplicated(mouse_map$human_gene), ]

## Senger: six fetal and three adult enterosphere RPKM profiles.
senger_dir <- file.path(data_root, "raw", "GSE101531")
if (!length(list.files(senger_dir, pattern = "RPKM.txt.gz$"))) {
  untar(file.path(senger_dir, "GSE101531_RAW.tar"), exdir = senger_dir)
}
senger_files <- list.files(senger_dir, pattern = "RPKM.txt.gz$", full.names = TRUE)
read_senger <- function(path) {
  x <- read.delim(gzfile(path), check.names = FALSE)
  setNames(x[[3]], x[[2]])
}
senger_list <- lapply(senger_files, read_senger)
senger_genes <- Reduce(intersect, lapply(senger_list, names))
senger_genes <- unique(senger_genes[!is.na(senger_genes) & nzchar(senger_genes)])
senger_mat <- do.call(cbind, lapply(senger_list, function(x) x[senger_genes]))
colnames(senger_mat) <- sub(".RPKM.txt.gz", "", sub("^GSM[0-9]+_", "", basename(senger_files)))
senger_stage <- ifelse(grepl("^AD", colnames(senger_mat)), "adult", "fetal")
stopifnot(sum(senger_stage == "fetal") == 6L, sum(senger_stage == "adult") == 3L)
senger_fit <- eBayes(lmFit(log2(senger_mat + 1), model.matrix(~factor(senger_stage, levels = c("adult", "fetal")))))
senger_de <- topTable(senger_fit, coef = 2, number = Inf, sort.by = "none")
senger_de$gene <- rownames(senger_de)
to_current <- setNames(aliases$current_hgnc_symbol, aliases$source_symbol)
renamed <- senger_de$gene %in% names(to_current) & !(to_current[senger_de$gene] %in% senger_de$gene)
senger_de$gene[renamed] <- to_current[senger_de$gene[renamed]]
senger_de <- senger_de[c("gene", "logFC", "P.Value", "adj.P.Val")]
names(senger_de)[2:4] <- c("Human_Senger_log2FC", "Human_Senger_PValue", "Human_Senger_FDR")

## Fawkner-Corbett: fetal EpCAM+ expression validation only.
fawkner_files <- list.files(file.path(data_root, "raw", "GSE158702"), pattern = "EPI.*tar.gz$", full.names = TRUE)
fawkner_rows <- list()
for (archive in fawkner_files) {
  pool <- sub("^GSM[0-9]+_", "", sub(".tar.gz", "", basename(archive)))
  extract_dir <- file.path(data_root, "processed", "GSE158702", pool)
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE)
    untar(archive, exdir = extract_dir)
  }
  mtx_path <- list.files(extract_dir, "matrix.mtx.gz$", recursive = TRUE, full.names = TRUE)
  feat_path <- list.files(extract_dir, "features.tsv.gz$", recursive = TRUE, full.names = TRUE)
  x <- readMM(gzfile(mtx_path[[1]]))
  features <- read.delim(gzfile(feat_path[[1]]), header = FALSE, stringsAsFactors = FALSE)
  total <- Matrix::colSums(x)
  detected <- Matrix::colSums(x > 0)
  keep <- total >= 500 & detected >= 200
  x <- x[, keep, drop = FALSE]
  symbols <- features[[2]]
  for (g in candidates$gene) {
    idx <- which(symbols == g)
    if (!length(idx)) idx <- which(symbols == source_symbol(g))
    # A gene absent from the feature list is not measured (NA), not zero.
    counts <- if (length(idx)) Matrix::colSums(x[idx, , drop = FALSE]) else NULL
    fawkner_rows[[length(fawkner_rows) + 1L]] <- data.frame(
      gene = g, pool = pool, n_cells = ncol(x),
      fetal_mean_UMI_per_10k = if (is.null(counts)) NA_real_ else mean(counts / Matrix::colSums(x) * 1e4),
      fetal_detection_fraction = if (is.null(counts)) NA_real_ else mean(counts > 0)
    )
  }
}
fawkner_long <- do.call(rbind, fawkner_rows)
write.csv(fawkner_long, file.path(out_dir, "Fawkner_fetal_expression_by_pool.csv"), row.names = FALSE)
fawkner <- aggregate(cbind(fetal_mean_UMI_per_10k, fetal_detection_fraction) ~ gene,
                     fawkner_long, median, na.action = na.pass)
names(fawkner)[2:3] <- c("Human_Fawkner_fetal_median_UMI_per_10k",
                         "Human_Fawkner_fetal_median_detection_fraction")
fawkner$Human_Fawkner_log2FC <- NA_real_
fawkner$Human_Fawkner_status <- "fetal_expression_only_no_matched_adult"

## Pikkupeura: independent fetal/adult contrasts within LN and collagen arms.
pik <- read.delim(gzfile(file.path(data_root, "raw", "GSE160449", "GSE160449_org_ln_rnaseq_counts.txt.gz")),
                  check.names = FALSE)
rownames(pik) <- pik$ENTREZID
pik_counts <- as.matrix(pik[, -1])
run_edger <- function(cols, stage) {
  y <- DGEList(pik_counts[, cols, drop = FALSE])
  keep <- filterByExpr(y, group = stage)
  y <- estimateDisp(calcNormFactors(y[keep, , keep.lib.sizes = FALSE]), model.matrix(~stage))
  fit <- glmQLFit(y, model.matrix(~stage), robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  tt$mouse_entrez_id <- rownames(tt)
  tt
}
## Reorder columns so adult precedes fetal for the model while retaining positive=fetal.
run_arm <- function(arm) {
  cols <- c(paste0("adult_", arm, "_", 1:3), paste0("fetal_", arm, "_", 1:3))
  run_edger(cols, factor(rep(c("adult", "fetal"), each = 3), levels = c("adult", "fetal")))
}
pik_ln <- run_arm("ln")
pik_c <- run_arm("c")
pik_small <- merge(pik_ln[c("mouse_entrez_id", "logFC", "FDR")],
                   pik_c[c("mouse_entrez_id", "logFC", "FDR")], by = "mouse_entrez_id", all = TRUE,
                   suffixes = c("_LN", "_collagen"))
names(pik_small)[2:5] <- c("Mouse_Pikkupeura_LN_log2FC", "Mouse_Pikkupeura_LN_FDR",
                           "Mouse_Pikkupeura_collagen_log2FC", "Mouse_Pikkupeura_collagen_FDR")
pik_small$Mouse_Pikkupeura_log2FC <- rowMeans(pik_small[c("Mouse_Pikkupeura_LN_log2FC", "Mouse_Pikkupeura_collagen_log2FC")], na.rm = TRUE)
pik_small$Mouse_Pikkupeura_direction_consistent <- with(pik_small, Mouse_Pikkupeura_LN_log2FC * Mouse_Pikkupeura_collagen_log2FC > 0)

## GSE44433: normalized matrix; WT fetal five versus WT adult five only.
gse_file <- file.path(data_root, "raw", "GSE44433", "GSE44433_series_matrix.txt.gz")
lines <- readLines(gzfile(gse_file))
start <- grep("!series_matrix_table_begin", lines) + 1L
end <- grep("!series_matrix_table_end", lines) - 1L
gse <- read.delim(text = paste(lines[start:end], collapse = "\n"), check.names = FALSE)
rownames(gse) <- sub("_at$", "", gse$ID_REF)
fetal_ids <- paste0("GSM108507", 5:9)
adult_ids <- paste0("GSM108509", 5:9)
gse_mat <- log2(as.matrix(gse[c(fetal_ids, adult_ids)]) + 1)
stage <- factor(rep(c("fetal", "adult"), each = 5), levels = c("adult", "fetal"))
gse_fit <- eBayes(lmFit(gse_mat, model.matrix(~stage)))
gse_de <- topTable(gse_fit, coef = 2, number = Inf, sort.by = "none")
gse_de$mouse_entrez_id <- rownames(gse_de)
gse_de <- gse_de[c("mouse_entrez_id", "logFC", "P.Value", "adj.P.Val")]
names(gse_de)[2:4] <- c("Mouse_GSE44433_log2FC", "Mouse_GSE44433_PValue", "Mouse_GSE44433_FDR")

benchmark <- merge(candidates, senger_de, by = "gene", all.x = TRUE)
benchmark <- merge(benchmark, fawkner, by = "gene", all.x = TRUE)
benchmark <- merge(benchmark, mouse_map, by.x = "gene", by.y = "human_gene", all.x = TRUE)
benchmark$mouse_entrez_id <- as.character(as.integer(benchmark$mouse_entrez_id))
benchmark <- merge(benchmark, pik_small, by = "mouse_entrez_id", all.x = TRUE)
benchmark <- merge(benchmark, gse_de, by = "mouse_entrez_id", all.x = TRUE)
priority <- c("TACSTD2", "GJA1", "CLU", "ANXA1", "TNFRSF12A")
benchmark$priority_marker <- benchmark$gene %in% priority
benchmark <- benchmark[match(candidates$gene, benchmark$gene), ]

write.csv(benchmark, file.path(out_dir, "Literature_31_four_dataset_benchmark.csv"), row.names = FALSE, na = "NA")
write.csv(benchmark[benchmark$priority_marker, ], file.path(out_dir, "Priority_5_marker_benchmark.csv"), row.names = FALSE, na = "NA")

summary <- data.frame(
  dataset = c("Human_Senger", "Human_Fawkner", "Mouse_Pikkupeura", "Mouse_GSE44433"),
  comparison = c("fetal_vs_adult", "fetal_expression_only", "fetal_vs_adult_two_culture_arms", "WT_fetal_vs_WT_adult"),
  n_candidates_with_measure = c(sum(!is.na(benchmark$Human_Senger_log2FC)),
                                sum(!is.na(benchmark$Human_Fawkner_fetal_median_detection_fraction)),
                                sum(!is.na(benchmark$Mouse_Pikkupeura_log2FC)),
                                sum(!is.na(benchmark$Mouse_GSE44433_log2FC))),
  n_positive_direction = c(sum(benchmark$Human_Senger_log2FC > 0, na.rm = TRUE), NA,
                           sum(benchmark$Mouse_Pikkupeura_log2FC > 0, na.rm = TRUE),
                           sum(benchmark$Mouse_GSE44433_log2FC > 0, na.rm = TRUE)),
  qualification_role = c("direct_human_culture_benchmark", "fetal_expression_validation_only",
                         "mouse_culture_benchmark", "mouse_in_vivo_replication")
)
write.csv(summary, file.path(out_dir, "Dataset_qualification_summary.csv"), row.names = FALSE, na = "NA")

metadata_audit <- data.frame(
  dataset = c("GSE101531", "GSE158702", "GSE160449", "GSE44433"),
  species = c("human", "human", "mouse", "mouse"),
  assay = c("bulk_RNAseq_RPKM", "scRNAseq_raw_UMI", "bulk_RNAseq_counts", "microarray_normalized_matrix"),
  material = c("epithelial_enterospheres", "primary_fetal_EpCAM_positive_epithelium",
               "fetal_and_adult_epithelial_cultures", "laser_microdissected_distal_ileal_epithelium"),
  fetal_units = c(6, 4, 6, 5),
  adult_units = c(3, 0, 6, 5),
  adult_reference_status = c("matched_within_study", "absent", "matched_within_arm", "matched_WT_within_study"),
  benchmark_use = c("fetal_adult_FC", "fetal_expression_only", "two_arm_fetal_adult_FC", "WT_fetal_adult_FC"),
  primary_limitation = c("culture_state", "no_adult_and_pooled_libraries", "culture_state_and_matrix_arm",
                         "older_microarray_and_partial_marker_coverage")
)
write.csv(metadata_audit, file.path(out_dir, "Dataset_metadata_qualification_audit.csv"), row.names = FALSE)

priority_recovery <- data.frame(
  dataset = c("Human_Senger", "Human_Fawkner", "Mouse_Pikkupeura", "Mouse_GSE44433"),
  rule = c("log2FC_gt_0", "detection_fraction_gt_0", "mean_log2FC_gt_0", "log2FC_gt_0"),
  measured_all_31 = c(sum(!is.na(benchmark$Human_Senger_log2FC)),
                      sum(!is.na(benchmark$Human_Fawkner_fetal_median_detection_fraction)),
                      sum(!is.na(benchmark$Mouse_Pikkupeura_log2FC)),
                      sum(!is.na(benchmark$Mouse_GSE44433_log2FC))),
  positive_all_31 = c(sum(benchmark$Human_Senger_log2FC > 0, na.rm = TRUE),
                      sum(benchmark$Human_Fawkner_fetal_median_detection_fraction > 0, na.rm = TRUE),
                      sum(benchmark$Mouse_Pikkupeura_log2FC > 0, na.rm = TRUE),
                      sum(benchmark$Mouse_GSE44433_log2FC > 0, na.rm = TRUE)),
  measured_priority_5 = c(sum(!is.na(benchmark$Human_Senger_log2FC[benchmark$priority_marker])),
                          sum(!is.na(benchmark$Human_Fawkner_fetal_median_detection_fraction[benchmark$priority_marker])),
                          sum(!is.na(benchmark$Mouse_Pikkupeura_log2FC[benchmark$priority_marker])),
                          sum(!is.na(benchmark$Mouse_GSE44433_log2FC[benchmark$priority_marker]))),
  positive_priority_5 = c(sum(benchmark$Human_Senger_log2FC[benchmark$priority_marker] > 0, na.rm = TRUE),
                          sum(benchmark$Human_Fawkner_fetal_median_detection_fraction[benchmark$priority_marker] > 0, na.rm = TRUE),
                          sum(benchmark$Mouse_Pikkupeura_log2FC[benchmark$priority_marker] > 0, na.rm = TRUE),
                          sum(benchmark$Mouse_GSE44433_log2FC[benchmark$priority_marker] > 0, na.rm = TRUE))
)
write.csv(priority_recovery, file.path(out_dir, "Marker_recovery_summary.csv"), row.names = FALSE)
message("Wrote 31-marker dataset qualification benchmark to: ", out_dir)
