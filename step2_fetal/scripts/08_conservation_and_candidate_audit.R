#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: 08_conservation_and_candidate_audit.R",
    "<tables_dir> <candidate_tsv> <figures_dir> <source_data_dir>"
  ))
}

tables_dir <- args[[1]]
candidate_path <- args[[2]]
figures_dir <- args[[3]]
source_data_dir <- args[[4]]
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

as_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  tolower(as.character(x)) == "true"
}

cross <- read.csv(
  file.path(tables_dir, "Cross_species_evidence.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
hgca <- read.csv(
  file.path(tables_dir, "human_HGCA_fetal_vs_adult_DEG.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gao <- read.csv(
  file.path(tables_dir, "human_Gao_validation.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
candidates <- read.delim(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(candidates) == 31L, !anyDuplicated(candidates$gene))

hgca_small <- hgca %>%
  transmute(
    gene,
    HGCA_log2FC,
    HGCA_FDR,
    HGCA_primary_pass = as_bool(HGCA_primary_pass),
    HGCA_proliferative_log2FC,
    HGCA_proliferative_direction_pass = as_bool(HGCA_proliferative_direction_pass),
    H1_pass = as_bool(H1_pass)
  )

gao_small <- gao %>%
  transmute(
    gene,
    Gao_log2FC,
    Gao_direction_consistent = as_bool(Gao_direction_consistent),
    Gao_pass = as_bool(Gao_pass)
  )

mouse_small <- cross %>%
  transmute(
    gene = human_gene,
    human_ensembl_gene_id,
    mouse_ensembl_gene_id,
    mouse_gene,
    Mouse_log2FC,
    Mouse_FDR,
    Mouse_pass = as_bool(Mouse_pass),
    Conserved_Fetal_High = as_bool(Conserved_Fetal_High)
  )

audit <- candidates %>%
  left_join(hgca_small, by = "gene") %>%
  left_join(gao_small, by = "gene") %>%
  left_join(mouse_small, by = "gene") %>%
  mutate(
    first_failed_gate = case_when(
      is.na(HGCA_log2FC) ~ "HGCA_not_tested",
      !HGCA_primary_pass ~ "HGCA_primary",
      is.na(HGCA_proliferative_log2FC) ~ "HGCA_proliferative_not_tested",
      !HGCA_proliferative_direction_pass ~ "HGCA_proliferative_direction",
      is.na(Gao_log2FC) ~ "Gao_not_tested",
      Gao_log2FC < 0.5 ~ "Gao_effect_size",
      !Gao_direction_consistent ~ "Gao_adult_donor_direction",
      is.na(mouse_ensembl_gene_id) ~ "no_one_to_one_mouse_orthologue",
      is.na(Mouse_log2FC) ~ "Mouse_not_tested",
      !Mouse_pass ~ "Mouse_M1",
      Conserved_Fetal_High ~ "PASS_FINAL",
      TRUE ~ "unclassified"
    )
  ) %>%
  select(
    Gene = gene,
    Literature_label = literature_label,
    HGCA_log2FC,
    HGCA_FDR,
    HGCA_primary_pass,
    HGCA_proliferative_log2FC,
    HGCA_proliferative_direction_pass,
    H1_pass,
    Gao_log2FC,
    Gao_direction_consistent,
    Gao_pass,
    mouse_gene,
    Mouse_log2FC,
    Mouse_FDR,
    Mouse_pass,
    Final = Conserved_Fetal_High,
    first_failed_gate
  )

write.csv(
  audit,
  file.path(tables_dir, "Literature_31_gate_failure_audit.csv"),
  row.names = FALSE,
  na = "NA"
)

tested <- cross %>%
  filter(!is.na(HGCA_log2FC), !is.na(Mouse_log2FC)) %>%
  mutate(
    Human_H1 = as_bool(H1_pass),
    Mouse_M1 = as_bool(Mouse_pass),
    overlap_class = case_when(
      Human_H1 & Mouse_M1 ~ "Human H1 + Mouse M1",
      Human_H1 & !Mouse_M1 ~ "Human H1 only",
      !Human_H1 & Mouse_M1 ~ "Mouse M1 only",
      TRUE ~ "Neither"
    )
  )
stopifnot(!anyDuplicated(tested$human_ensembl_gene_id))

overlap_levels <- c("Neither", "Human H1 only", "Mouse M1 only", "Human H1 + Mouse M1")
tested$overlap_class <- factor(tested$overlap_class, levels = overlap_levels)

contingency <- table(
  Human_fetal_high = factor(tested$Human_H1, levels = c(TRUE, FALSE)),
  Mouse_fetal_high = factor(tested$Mouse_M1, levels = c(TRUE, FALSE))
)
fisher_result <- fisher.test(contingency, alternative = "greater")
spearman_result <- cor.test(
  tested$HGCA_log2FC,
  tested$Mouse_log2FC,
  method = "spearman",
  exact = FALSE
)

both_count <- sum(tested$Human_H1 & tested$Mouse_M1)
human_count <- sum(tested$Human_H1)
mouse_count <- sum(tested$Mouse_M1)

stats <- data.frame(
  n_tested_one_to_one = nrow(tested),
  spearman_rho = unname(spearman_result$estimate),
  spearman_p_value = spearman_result$p.value,
  fisher_odds_ratio = unname(fisher_result$estimate),
  fisher_p_value = fisher_result$p.value,
  human_H1_count = human_count,
  mouse_M1_count = mouse_count,
  both_H1_M1_count = both_count,
  fraction_human_H1_also_mouse_M1 = both_count / human_count,
  fraction_mouse_M1_also_human_H1 = both_count / mouse_count
)

overlap_table <- as.data.frame(contingency) %>%
  rename(count = Freq) %>%
  mutate(
    Human_fetal_high = as.character(Human_fetal_high),
    Mouse_fetal_high = as.character(Mouse_fetal_high)
  )

write.csv(
  tested,
  file.path(source_data_dir, "human_mouse_all_tested_one2one_source_data.csv"),
  row.names = FALSE
)
write.csv(
  overlap_table,
  file.path(tables_dir, "Human_mouse_fetal_high_overlap_2x2.csv"),
  row.names = FALSE
)
write.csv(
  stats,
  file.path(tables_dir, "Human_mouse_conservation_statistics.csv"),
  row.names = FALSE
)

label_genes <- c("TACSTD2", "GJA1", "ANXA1", "CLU", "TNFRSF12A", "EMP1", "LAMC2")
label_data <- tested %>% filter(human_gene %in% label_genes)

palette <- c(
  "Neither" = "#C7C7C7",
  "Human H1 only" = "#E28E2C",
  "Mouse M1 only" = "#33A6A6",
  "Human H1 + Mouse M1" = "#C84B31"
)

p <- ggplot(tested, aes(HGCA_log2FC, Mouse_log2FC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888888") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
  geom_point(aes(colour = overlap_class), size = 0.75, alpha = 0.55, stroke = 0) +
  geom_point(
    data = label_data,
    shape = 21,
    size = 2.1,
    stroke = 0.5,
    colour = "black",
    fill = "white"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = human_gene),
    size = 2.2,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.size = 0.22,
    seed = 20261002,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(
    title = "Genome-wide human–mouse fetal effect-size concordance",
    subtitle = sprintf(
      "%s tested one-to-one genes; Spearman ρ = %.2f; Fisher OR = %.2f",
      format(nrow(tested), big.mark = ","), stats$spearman_rho, stats$fisher_odds_ratio
    ),
    x = expression("Human HGCA " * log[2] * " fold change (fetal / adult)"),
    y = expression("Mouse in vivo " * log[2] * " fold change (fetal / adult)"),
    colour = NULL,
    caption = "Colours use frozen H1 and M1 calls; correlation and enrichment are diagnostic, not gates."
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.1, alpha = 1))) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6.5, colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 6.5, colour = "#444444"),
    plot.caption = element_text(size = 5.7, colour = "#555555", hjust = 0),
    plot.margin = margin(6, 8, 6, 6)
  )

width_in <- 140 / 25.4
height_in <- 118 / 25.4
base_path <- file.path(figures_dir, "human_mouse_all_tested_conservation_scatter")
ggsave(paste0(base_path, ".pdf"), p, width = width_in, height = height_in,
       device = cairo_pdf, units = "in", bg = "white")
ggsave(paste0(base_path, ".png"), p, width = width_in, height = height_in,
       device = ragg::agg_png, units = "in", dpi = 600, bg = "white")

message(sprintf(
  paste(
    "Audited 31 candidates; tested universe = %d; Spearman rho = %.3f;",
    "Fisher OR = %.3f, P = %.3g; overlap = %d."
  ),
  nrow(tested), stats$spearman_rho, stats$fisher_odds_ratio,
  stats$fisher_p_value, both_count
))
