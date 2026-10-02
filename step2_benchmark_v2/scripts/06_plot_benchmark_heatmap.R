#!/usr/bin/env Rscript
# 31-marker x contrast log2FC heatmap (fetal vs adult), priority controls on
# top, nominal P < 0.05 marked, qualification verdict in the column labels.
# Usage: Rscript 06_plot_benchmark_heatmap.R <tables_dir> <figures_dir> <source_data_dir>

suppressPackageStartupMessages({
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: 06_plot_benchmark_heatmap.R <tables> <figures> <source_data>")
tab <- args[[1]]; fig <- args[[2]]; src <- args[[3]]
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
dir.create(src, showWarnings = FALSE, recursive = TRUE)

m <- read.csv(file.path(tab, "Benchmark_v2_31_marker_matrix.csv"), check.names = FALSE)
q <- read.csv(file.path(tab, "Benchmark_v2_qualification.csv"))
keys <- q$contrast
long <- do.call(rbind, lapply(keys, function(k) {
  p <- m[[paste0(k, "__PValue")]]
  fdr <- m[[paste0(k, "__FDR")]]
  data.frame(marker = m$marker, contrast = k,
             log2FC = m[[paste0(k, "__log2FC")]],
             support = ifelse(is.na(p), fdr, p))
}))
verdict <- setNames(paste0(q$contrast, "\n", ifelse(q$qualifies, "QUALIFIES", "fails"),
                           " (", q$markers_fetal_positive, "/", q$markers_measured, " +)"), q$contrast)
long$contrast_label <- factor(verdict[long$contrast], levels = verdict[keys])
long$marker <- factor(long$marker, levels = rev(m$marker))
long$capped <- pmax(pmin(long$log2FC, 4), -4)
long$sig <- ifelse(!is.na(long$support) & long$support < 0.05, "*", "")
write.csv(long, file.path(src, "Benchmark_v2_heatmap_source_data.csv"), row.names = FALSE)

p <- ggplot(long, aes(contrast_label, marker, fill = capped)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3, vjust = 0.75) +
  geom_hline(yintercept = nrow(m) - 3.5, linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-4, 4), na.value = "grey85", name = "log2FC\n(fetal/adult)\ncapped ±4") +
  labs(title = "Step 2 human benchmark v2: 31 literature markers",
       subtitle = "* P < 0.05 (Pikkupeura: FDR); grey = not measured; top four = priority controls",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(fig, "Benchmark_v2_31_marker_heatmap.pdf"), p, width = 3 + 0.5 * length(keys), height = 9, device = cairo_pdf)
ggsave(file.path(fig, "Benchmark_v2_31_marker_heatmap.png"), p, width = 3 + 0.5 * length(keys), height = 9, dpi = 600,
       device = ragg::agg_png)
cat("Heatmap written\n")
