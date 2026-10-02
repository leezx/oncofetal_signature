#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 07_plot_human_mouse_effect_scatter.R <conserved_csv> <figures_dir> <source_data_dir>")
}

input_path <- args[[1]]
figures_dir <- args[[2]]
source_data_dir <- args[[3]]
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

requested_markers <- c("TACSTD2", "GJA1", "ANXA1", "CLU", "TNFRSF12A", "EMP1", "LAMC2")
dat <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
conserved_flag <- tolower(as.character(dat$Conserved_Fetal_High)) == "true"
stopifnot(nrow(dat) == 706L, all(conserved_flag))

plot_data <- dat %>%
  transmute(
    human_gene,
    mouse_gene,
    human_ensembl_gene_id,
    mouse_ensembl_gene_id,
    HGCA_log2FC,
    Mouse_log2FC,
    requested_marker = human_gene %in% requested_markers
  )

marker_status <- data.frame(human_gene = requested_markers) %>%
  left_join(
    plot_data %>%
      select(human_gene, mouse_gene, HGCA_log2FC, Mouse_log2FC) %>%
      mutate(in_final_706 = TRUE),
    by = "human_gene"
  ) %>%
  mutate(in_final_706 = ifelse(is.na(in_final_706), FALSE, in_final_706))

write.csv(
  plot_data,
  file.path(source_data_dir, "human_mouse_effect_scatter_source_data.csv"),
  row.names = FALSE
)
write.csv(
  marker_status,
  file.path(source_data_dir, "human_mouse_requested_marker_status.csv"),
  row.names = FALSE
)

rho <- cor(plot_data$HGCA_log2FC, plot_data$Mouse_log2FC, method = "spearman")
label_data <- plot_data %>% filter(requested_marker)

p <- ggplot(plot_data, aes(HGCA_log2FC, Mouse_log2FC)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = "dashed", colour = "#777777") +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed", colour = "#777777") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#B0B0B0") +
  geom_point(size = 1.4, alpha = 0.65, colour = "#3976A8") +
  geom_point(
    data = label_data,
    size = 2.1,
    colour = "#C84B31"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = human_gene),
    size = 2.5,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.25,
    seed = 20261002,
    show.legend = FALSE
  ) +
  labs(
    title = "Conserved fetal-high effect sizes across species",
    subtitle = sprintf(
      "706 deterministic-intersection genes; Spearman ρ = %.2f (descriptive only)",
      rho
    ),
    x = expression("Human HGCA " * log[2] * " fold change (fetal / adult)"),
    y = expression("Mouse in vivo " * log[2] * " fold change (fetal / adult)"),
    caption = "Dashed lines show the frozen log2FC = 0.5 gates; grey line is y = x."
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6.5, colour = "black"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, colour = "#444444"),
    plot.caption = element_text(size = 5.7, colour = "#555555", hjust = 0),
    plot.margin = margin(6, 8, 6, 6)
  )

width_in <- 120 / 25.4
height_in <- 105 / 25.4
base_path <- file.path(figures_dir, "human_mouse_706_effect_size_scatter")
ggsave(paste0(base_path, ".pdf"), p, width = width_in, height = height_in,
       device = cairo_pdf, units = "in", bg = "white")
ggsave(paste0(base_path, ".png"), p, width = width_in, height = height_in,
       device = ragg::agg_png, units = "in", dpi = 600, bg = "white")

message(sprintf(
  "Wrote 706-gene effect-size scatter; %d of %d requested markers are present; Spearman rho = %.3f.",
  nrow(label_data), length(requested_markers), rho
))
