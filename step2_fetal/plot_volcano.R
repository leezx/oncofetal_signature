#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1]] else "step2_fetal"
output_dir <- if (length(args) >= 2L) args[[2]] else input_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  "Fetal-high" = "#C84B31",
  "Adult-high" = "#3976A8",
  "Not significant" = "#C7C7C7"
)

classify_de <- function(log2fc, fdr) {
  case_when(
    !is.na(fdr) & fdr < 0.05 & log2fc >= 0.5 ~ "Fetal-high",
    !is.na(fdr) & fdr < 0.05 & log2fc <= -0.5 ~ "Adult-high",
    TRUE ~ "Not significant"
  )
}

prepare_volcano <- function(df, gene_col, fc_col, fdr_col) {
  df %>%
    transmute(
      gene = .data[[gene_col]],
      log2FC = .data[[fc_col]],
      FDR = .data[[fdr_col]],
      neg_log10_FDR = -log10(pmax(.data[[fdr_col]], .Machine$double.xmin)),
      status = factor(
        classify_de(.data[[fc_col]], .data[[fdr_col]]),
        levels = c("Not significant", "Adult-high", "Fetal-high")
      )
    )
}

select_labels <- function(df, n_each = 10L) {
  bind_rows(
    df %>% filter(status == "Fetal-high") %>% arrange(FDR, desc(log2FC)) %>% slice_head(n = n_each),
    df %>% filter(status == "Adult-high") %>% arrange(FDR, log2FC) %>% slice_head(n = n_each)
  ) %>% distinct(gene, .keep_all = TRUE)
}

theme_volcano <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.5, colour = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6.2),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.5, colour = "#444444"),
      plot.margin = margin(5, 7, 5, 5)
    )
}

make_volcano <- function(df, title) {
  labels <- select_labels(df)
  counts <- table(df$status)
  subtitle <- sprintf(
    "Fetal-high: %s  |  Adult-high: %s  |  thresholds: |log2FC| ≥ 0.5, FDR < 0.05",
    format(counts[["Fetal-high"]], big.mark = ","),
    format(counts[["Adult-high"]], big.mark = ",")
  )
  ggplot(df, aes(log2FC, neg_log10_FDR)) +
    geom_point(aes(colour = status), size = 0.7, alpha = 0.65, stroke = 0) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.3, colour = "#606060") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "#606060") +
    geom_text_repel(
      data = labels,
      aes(label = gene),
      size = 2,
      min.segment.length = 0,
      segment.size = 0.22,
      segment.color = "#777777",
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = 20261002,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(log[2] * " fold change (fetal / adult)"),
      y = expression(-log[10] * " FDR")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    theme_volcano()
}

save_plot <- function(plot, filename, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(paste0(filename, ".pdf"), plot, width = width_in, height = height_in,
         device = cairo_pdf, units = "in", bg = "white")
  ggsave(paste0(filename, ".png"), plot, width = width_in, height = height_in,
         device = ragg::agg_png, units = "in", dpi = 600, bg = "white")
}

human_raw <- read.csv(file.path(input_dir, "human_HGCA_fetal_vs_adult_DEG.csv"), check.names = FALSE)
mouse_raw <- read.csv(file.path(input_dir, "mouse_in_vivo_fetal_vs_adult_DEG.csv"), check.names = FALSE)
mouse_raw$plot_gene <- ifelse(
  is.na(mouse_raw$mouse_gene) | mouse_raw$mouse_gene == "",
  paste0("Entrez:", mouse_raw$mouse_entrez_id),
  mouse_raw$mouse_gene
)

human <- prepare_volcano(human_raw, "gene", "HGCA_log2FC", "HGCA_FDR")
mouse <- prepare_volcano(mouse_raw, "plot_gene", "Mouse_log2FC", "Mouse_FDR")

write.csv(human, file.path(output_dir, "human_volcano_source_data.csv"), row.names = FALSE)
write.csv(mouse, file.path(output_dir, "mouse_volcano_source_data.csv"), row.names = FALSE)

p_human <- make_volcano(human, "Human HGCA: fetal vs adult intestinal epithelium")
p_mouse <- make_volcano(mouse, "Mouse in vivo: fetal epithelium vs adult crypt")
combined <- p_human + p_mouse +
  plot_layout(guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 8, face = "bold"))

save_plot(p_human, file.path(output_dir, "human_HGCA_fetal_vs_adult_volcano"), 120, 100)
save_plot(p_mouse, file.path(output_dir, "mouse_in_vivo_fetal_vs_adult_volcano"), 120, 100)
save_plot(combined, file.path(output_dir, "human_mouse_fetal_vs_adult_volcano"), 183, 100)

message("Volcano figures and source-data tables written to: ", normalizePath(output_dir))
