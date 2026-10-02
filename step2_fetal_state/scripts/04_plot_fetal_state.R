#!/usr/bin/env Rscript
# Figures for the fetal epithelial state analysis (plan v1, addenda v1.1-1.2).
# 1. Unit-level heatmap (rows = donor x region units; per-gene z within dataset).
# 2. TNFRSF12A by epithelial state x age; one point = one unit x state
#    pseudobulk (never a cell), with % positive cells alongside.
# 3. The same for TACSTD2, CLU, ANXA1 (Fawkner; Gao states are exploratory).
# Usage: Rscript 04_plot_fetal_state.R <tables_dir> <figures_dir> <source_data_dir>
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(ggplot2); library(grid)
})
args <- commandArgs(trailingOnly = TRUE)
tab <- args[[1]]; fig <- args[[2]]; src <- args[[3]]
genes <- c("TNFRSF12A", "TACSTD2", "CLU", "ANXA1", "EPCAM", "OLFM4", "LGR5", "ASCL2", "MKI67",
           "YAP1", "CCN2", "CCN1", "ANKRD1")
role <- c(rep("test", 4), rep("control", 5), rep("YAP context", 4))
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 600)
}

# ---- 1. heatmap ----
u <- read.csv(file.path(tab, "Unit_gene_values.csv"), check.names = FALSE)
u <- u[as.character(u$eligible) %in% c("TRUE", "True"), ]
u <- u[order(u$dataset, u$pcw, u$donor, u$region), ]
z <- do.call(rbind, lapply(split(u, u$dataset), function(d) {
  m <- as.matrix(d[, genes]); scale(m)
}))
z[is.nan(z)] <- 0
rownames(z) <- paste0(u$unit, "  (", u$pcw, " wk, ", u$region, ")")
write.csv(cbind(u[, c("dataset", "unit", "donor", "pcw", "region", "n_cells")], u[, genes],
                setNames(as.data.frame(z), paste0(genes, "_z"))),
          file.path(src, "Fig1_unit_heatmap_source.csv"), row.names = FALSE)
reg_col <- c("TI" = "#1b9e77", "SI" = "#1b9e77", "proximal colon" = "#d95f02",
             "distal colon" = "#7570b3", "LI" = "#7570b3", "hindgut" = "#e7298a")
ra <- rowAnnotation(
  dataset = u$dataset, age = u$pcw, region = u$region, cells = anno_barplot(log10(u$n_cells), width = unit(1, "cm")),
  col = list(dataset = c(Fawkner2021 = "#386cb0", Gao2018 = "#bf5b17"),
             age = colorRamp2(c(6, 13, 25), c("#fff7bc", "#fe9929", "#662506")), region = reg_col),
  annotation_label = c("dataset", "age (wk)", "region", "log10 cells"))
ca <- HeatmapAnnotation(role = role, col = list(role = c(test = "black", control = "grey60", `YAP context` = "#984ea3")))
hm <- Heatmap(z, name = "z (within\ndataset)", col = colorRamp2(c(-2.5, 0, 2.5), c("#2166ac", "white", "#b2182b")),
              cluster_rows = FALSE, cluster_columns = FALSE, row_split = u$dataset, row_gap = unit(3, "mm"),
              column_split = factor(role, levels = unique(role)), top_annotation = ca, left_annotation = ra,
              row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
              cell_fun = function(j, i, x, y, w, h, fill) NULL,
              column_title = "Unit-level expression (Fawkner: pseudobulk log2 CPM+1; Gao: log2 mean TPM+1)",
              column_title_gp = gpar(fontsize = 10))
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf(file.path(fig, "Fig1_unit_gene_heatmap.pdf"), width = 9.5, height = 10)
  else png(file.path(fig, "Fig1_unit_gene_heatmap.png"), width = 9.5, height = 10, units = "in", res = 600)
  draw(hm, merge_legend = TRUE); invisible(dev.off())
}

# ---- 2/3. state x age ----
us <- read.csv(file.path(tab, "Unit_state_gene_values.csv"), check.names = FALSE)
us <- us[as.character(us$eligible) %in% c("TRUE", "True"), ]
state_order <- c("Stem", "TA", "Fetal progenitor", "Early enterocyte", "Mature enterocyte",
                 "BEST4/OTOP2", "Goblet", "EEC/secretory progenitor", "Unresolved")
us$state <- factor(us$state, levels = state_order)
us$region_group <- ifelse(us$region %in% c("TI", "SI"), "small intestine (TI/SI)",
                          ifelse(us$region == "hindgut", "hindgut", "colon (prox/dist/LI)"))
us$dataset_label <- ifelse(us$dataset == "Gao2018", "Gao 2018 (states exploratory)", "Fawkner 2021")
long <- do.call(rbind, lapply(c("TNFRSF12A", "TACSTD2", "CLU", "ANXA1"), function(g) {
  rbind(data.frame(us[, c("dataset", "dataset_label", "unit", "donor", "pcw", "region", "region_group", "state",
                          "n_cells_state")], gene = g, measure = "level (unit x state pseudobulk)", value = us[[g]]),
        data.frame(us[, c("dataset", "dataset_label", "unit", "donor", "pcw", "region", "region_group", "state",
                          "n_cells_state")], gene = g, measure = "% positive cells", value = us[[paste0(g, "_pct_pos")]]))
}))
write.csv(long, file.path(src, "Fig2_Fig3_state_by_age_source.csv"), row.names = FALSE)
pal <- c("small intestine (TI/SI)" = "#1b9e77", "colon (prox/dist/LI)" = "#7570b3", hindgut = "#e7298a")
plot_gene <- function(d, title) {
  ggplot(d, aes(pcw, value)) +
    geom_point(aes(colour = region_group, size = n_cells_state), alpha = 0.85) +
    facet_grid(measure + dataset_label ~ state, scales = "free_y", switch = "y",
               labeller = labeller(state = label_wrap_gen(12), dataset_label = label_wrap_gen(14),
                                   measure = label_wrap_gen(14))) +
    scale_colour_manual(values = pal, name = "region") +
    scale_size_area(max_size = 3.2, breaks = c(10, 100, 1000), name = "cells") +
    labs(x = "age (weeks; Fawkner PCW, Gao as reported)", y = NULL, title = title,
         subtitle = "One point = one donor x region unit within a state (Fawkner >= 10 cells, Gao >= 5); no hashtag-failed or sex-discordant units") +
    theme_bw(base_size = 8) +
    theme(strip.placement = "outside", strip.background = element_rect(fill = "grey95"),
          legend.position = "bottom", panel.grid.minor = element_blank())
}
p2 <- plot_gene(long[long$gene == "TNFRSF12A", ],
                "TNFRSF12A by epithelial state and developmental age (Fawkner: log2 CPM+1; Gao: log2 mean TPM+1)")
save_both(p2, "Fig2_TNFRSF12A_state_by_age", 13, 8.5)
for (g in c("TACSTD2", "CLU", "ANXA1")) {
  d <- long[long$gene == g & long$dataset == "Fawkner2021", ]
  save_both(plot_gene(d, paste(g, "by epithelial state and age (Fawkner 2021)")),
            paste0("Fig3_", g, "_state_by_age_Fawkner"), 13, 5)
}
writeLines(capture.output(sessionInfo()), file.path(src, "sessionInfo_04_plot.txt"))
