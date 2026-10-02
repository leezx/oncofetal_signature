#!/usr/bin/env Rscript
# Volcano plots for each gated contrast, the bulk-vs-epithelial effect-size
# scatter, and the dataset admission-QA panel. Writes PDF, 600-dpi PNG and
# exact plotting source data.
#
# Usage: Rscript 06_plot_cancer_axis.R <tables_dir> <figures_dir> <source_data_dir> <min_log2fc> <max_fdr> [public]
# `public` draws only TCGA/Pelka content (no Joanito values or Joanito-derived
# labels) into *_TCGA_Pelka files that may be committed to the public repo.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 5:6) stop("Usage: 06_plot_cancer_axis.R <tables> <figures> <source_data> <min_lfc> <max_fdr> [public]")
public <- length(args) == 6 && args[[6]] == "public"
tab <- args[[1]]; fig <- args[[2]]; src <- args[[3]]
min_lfc <- as.numeric(args[[4]]); max_fdr <- as.numeric(args[[5]])
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
dir.create(src, showWarnings = FALSE, recursive = TRUE)

ev <- read.csv(file.path(tab, "CRC_high_evidence.csv"), check.names = FALSE)
ev$symbol <- ev$hgnc_symbol  # current HGNC symbols in figures
qa <- read.csv(file.path(tab, "Dataset_admission_QA.csv"))
panel <- c("CDH3", "CLDN1", "FOXQ1", "KRT23", "LGR5", "ASCL2", "MYC", "TESC",
           "CA1", "CA2", "CA4", "GUCA2A", "GUCA2B", "CLCA4", "AQP8", "SLC26A3", "MS4A12", "CEACAM7")
focus <- c("TACSTD2", "ANXA1", "CLU", "TNFRSF12A", "EMP1", "CCN1", "CCN2", "GJA1", "LY6D")

save_plot <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 600, device = ragg::agg_png)
}

volcano <- function(lfc, fdr, title, name, up_label) {
  d <- ev[!is.na(ev[[lfc]]) & !is.na(ev[[fdr]]), c("gene_id", "symbol", lfc, fdr)]
  names(d) <- c("gene_id", "symbol", "log2FC", "FDR")
  d$neglog10FDR <- pmin(-log10(pmax(d$FDR, 1e-300)), 300)
  d$class <- ifelse(d$log2FC >= min_lfc & d$FDR < max_fdr, up_label,
                    ifelse(d$log2FC <= -min_lfc & d$FDR < max_fdr, "Normal-high", "NS"))
  d$label <- ifelse(d$symbol %in% focus, d$symbol, NA)
  write.csv(d, file.path(src, paste0(name, "_source_data.csv")), row.names = FALSE)
  n_up <- sum(d$class == up_label); n_dn <- sum(d$class == "Normal-high")
  cols <- setNames(c("#B2182B", "#2166AC", "grey75"), c(up_label, "Normal-high", "NS"))
  p <- ggplot(d, aes(log2FC, neglog10FDR, colour = class)) +
    geom_point(size = 0.4, alpha = 0.6, stroke = 0) +
    geom_vline(xintercept = c(-min_lfc, min_lfc), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(max_fdr), linetype = 2, linewidth = 0.3) +
    geom_text_repel(aes(label = label), colour = "black", size = 2.6, max.overlaps = Inf,
                    min.segment.length = 0, na.rm = TRUE) +
    scale_colour_manual(values = cols) +
    labs(title = title, subtitle = sprintf("%s: %d   Normal-high: %d", up_label, n_up, n_dn),
         x = "log2 fold change (tumour / normal)", y = expression(-log[10]~FDR), colour = NULL) +
    theme_classic(base_size = 9) + theme(legend.position = "bottom")
  save_plot(p, name, 4.2, 4.4)
}

volcano("T1_log2FC", "T1_FDR", "T1  TCGA COAD+READ: primary tumour vs normal",
        "TCGA_tumor_vs_normal_volcano", "Tumour-high")
volcano("Pelka_log2FC", "Pelka_FDR", "S2  Pelka: tumour vs normal epithelium",
        "Pelka_tumor_vs_normal_epithelium_volcano", "Tumour-high")
has_s1 <- !public && "Joanito_log2FC" %in% names(ev) && any(!is.na(ev$Joanito_log2FC))
sfx <- if (public) "_TCGA_Pelka" else ""
if (public) {
  # Colour by non-restricted evidence only.
  t1 <- as.character(ev$bulk_support) %in% c("TRUE", "True")
  s2 <- as.character(ev$S2_pass) %in% c("TRUE", "True")
  ev$final_label <- ifelse(t1 & s2, "T1 and S2 pass",
                    ifelse(s2, "S2 pass only", ifelse(t1, "T1 pass only", "Neither")))
}
if (has_s1) {
  volcano("Joanito_log2FC", "Joanito_FDR", "S1  Joanito: malignant vs normal epithelium",
          "Joanito_malignant_vs_normal_volcano", "Malignant-high")
}

# Bulk vs epithelial effect-size scatter over all jointly tested genes.
epi_col <- if (has_s1) "Joanito_log2FC" else "Pelka_log2FC"
epi_lab <- if (has_s1) "S1 Joanito malignant / normal epithelium" else "S2 Pelka tumour / normal epithelium"
d <- ev[!is.na(ev$T1_log2FC) & !is.na(ev[[epi_col]]), c("gene_id", "symbol", "final_label", "T1_log2FC", epi_col)]
names(d)[5] <- "epithelial_log2FC"
rho <- cor(d$T1_log2FC, d$epithelial_log2FC, method = "spearman")
d$label <- ifelse(d$symbol %in% focus, d$symbol, NA)
write.csv(d, file.path(src, paste0("TCGA_vs_epithelial_effect_scatter", sfx, "_source_data.csv")), row.names = FALSE)
p <- ggplot(d, aes(T1_log2FC, epithelial_log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(aes(colour = final_label), size = 0.4, alpha = 0.5, stroke = 0) +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = Inf, min.segment.length = 0, na.rm = TRUE) +
  labs(title = "Bulk tumour vs epithelial re-expression",
       subtitle = sprintf("%d jointly tested genes; Spearman rho = %.3f", nrow(d), rho),
       x = "T1 TCGA log2FC (tumour / normal)", y = epi_lab, colour = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 1))) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
save_plot(p, paste0("TCGA_vs_epithelial_effect_scatter", sfx), 5, 5.4)

# Admission-QA panel: log2FC of the frozen marker panel in each contrast.
cols <- c(T1_TCGA = "T1_log2FC", S2_Pelka = "Pelka_log2FC")
if (has_s1) cols <- c(cols, S1_Joanito = "Joanito_log2FC")
q <- do.call(rbind, lapply(names(cols), function(k) {
  x <- ev[ev$symbol %in% panel, c("symbol", cols[[k]])]
  x <- x[!duplicated(x$symbol), ]
  data.frame(contrast = k, symbol = x$symbol, log2FC = x[[2]])
}))
q$expected <- ifelse(q$symbol %in% panel[1:8], "Tumour-up panel", "Tumour-down panel")
q$symbol <- factor(q$symbol, levels = rev(panel))
write.csv(q, file.path(src, paste0("Dataset_admission_QA", sfx, "_source_data.csv")), row.names = FALSE)
p <- ggplot(q, aes(log2FC, symbol, fill = log2FC > 0)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(expected ~ contrast, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
  labs(title = "Dataset admission QA: canonical CRC-vs-normal markers", x = "log2FC (tumour / normal)", y = NULL) +
  theme_bw(base_size = 9)
save_plot(p, paste0("Dataset_admission_QA", sfx), 2.2 + 2 * length(cols), 4.6)
cat("Figures written to", fig, "\n")
