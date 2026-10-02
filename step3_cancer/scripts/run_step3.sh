#!/usr/bin/env bash
# Run the complete Step 3 cancer axis from frozen inputs. Joanito (S1) steps
# run only when its Synapse files are present; otherwise the integration
# writes S1-pending labels and makes no CRC_high call.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step_dir="$(cd "$script_dir/.." && pwd)"
source "$step_dir/config/step3.env"

de_dir="$DATA_RESULTS/de"
tables_dir="$step_dir/results/tables"
figures_dir="$step_dir/results/figures"
source_data_dir="$step_dir/results/source_data"
mkdir -p "$de_dir" "$tables_dir" "$figures_dir" "$source_data_dir"

# Acquisition (idempotent).
[[ -s "$RECOUNT3_DIR/processed/v0.1/recount3_TCGA_COAD_READ_GTEx_COLON_gene_counts.rds" ]] ||
  Rscript "$script_dir/01_recount3_tcga_gtex.R" "$RECOUNT3_DIR"
# shellcheck disable=SC2086
python3 "$script_dir/00_synapse_download.py" --out-dir "$JOANITO_RAW" $JOANITO_SYN_IDS ||
  echo "WARNING: Joanito download unavailable (Synapse access requirement?); S1 will be pending." >&2

# T1 — TCGA.
Rscript "$script_dir/03_tcga_tumor_vs_normal.R" \
  "$RECOUNT3_DIR/processed/v0.1" "$de_dir" "$MIN_LOG2FC" "$MAX_FDR"

# S2 and Pelka proliferation control.
python3 "$script_dir/02_pelka_pseudobulk.py" --raw-dir "$PELKA_RAW" --out-dir "$PELKA_PB" --min-cells "$MIN_CELLS"
for ref in Normal Normal_prolif; do
  Rscript "$script_dir/04_pseudobulk_edger.R" "$PELKA_PB/Pelka_epithelial_pseudobulk_counts.csv.gz" \
    "$PELKA_PB/Pelka_pseudobulk_samples.csv" Tumor "$ref" none "$de_dir/Pelka_Tumor_vs_${ref}.csv"
done

# S1 — Joanito (only when downloaded). Label map is frozen (plan v3).
if [[ -s "$JOANITO_RAW/Epithelial_Count_matrix.h5" ]]; then
  python3 "$script_dir/07_joanito_pseudobulk.py" --raw-dir "$JOANITO_RAW" \
    --label-map "$step_dir/config/joanito_label_map.tsv" --out-dir "$JOANITO_PB" --min-cells "$MIN_CELLS"
  python3 "$script_dir/08_joanito_cohort_check.py" prepare \
    --samples "$JOANITO_PB/Joanito_pseudobulk_samples.csv" \
    --contingency "$de_dir/Joanito_cohort_group_contingency.csv" \
    --sensitivity-samples "$JOANITO_PB/Joanito_pseudobulk_samples_paired_cohorts.csv"
  counts="$JOANITO_PB/Joanito_epithelial_pseudobulk_counts.csv.gz"
  Rscript "$script_dir/04_pseudobulk_edger.R" "$counts" "$JOANITO_PB/Joanito_pseudobulk_samples.csv" \
    Malignant Normal cohort "$de_dir/Joanito_Malignant_vs_Normal.csv"
  Rscript "$script_dir/04_pseudobulk_edger.R" "$counts" "$JOANITO_PB/Joanito_pseudobulk_samples_paired_cohorts.csv" \
    Malignant Normal cohort "$de_dir/Joanito_Malignant_vs_Normal_sensitivity_paired_cohorts.csv"
  python3 "$script_dir/08_joanito_cohort_check.py" compare \
    --full "$de_dir/Joanito_Malignant_vs_Normal.csv" \
    --sensitivity "$de_dir/Joanito_Malignant_vs_Normal_sensitivity_paired_cohorts.csv" \
    --out "$de_dir/Joanito_full_vs_sensitivity_concordance.csv" \
    --min-log2fc "$MIN_LOG2FC" --max-fdr "$MAX_FDR"
fi

# Integration, QA, audit, figures.
python3 "$script_dir/05_crc_high_integration.py" --de-dir "$de_dir" \
  --candidates "$step_dir/config/literature_candidates_31.tsv" \
  --aliases "$step_dir/config/symbol_aliases.tsv" \
  --orthology "$step_dir/config/nonhuman_candidate_orthology.tsv" --out-dir "$tables_dir" \
  --min-log2fc "$MIN_LOG2FC" --max-fdr "$MAX_FDR"
cp "$de_dir/TCGA_tumor_vs_normal_DEG.csv" "$de_dir/TCGA_sample_inclusion.csv" \
   "$de_dir/TCGA_contrast_summary.csv" "$tables_dir/"
cp "$de_dir/Pelka_Tumor_vs_Normal.csv" "$tables_dir/Pelka_tumor_vs_normal_epithelium_DEG.csv"
cp "$de_dir/Pelka_Tumor_vs_Normal_prolif.csv" "$tables_dir/Pelka_tumor_vs_normal_prolif_control_DEG.csv"
cp "$PELKA_PB/Pelka_pseudobulk_samples.csv" "$tables_dir/"
for f in "$de_dir"/Joanito_*.csv; do [[ -e "$f" ]] && cp "$f" "$tables_dir/"; done
[[ -e "$JOANITO_PB/Joanito_pseudobulk_samples.csv" ]] && cp "$JOANITO_PB/Joanito_pseudobulk_samples.csv" "$tables_dir/"

Rscript "$script_dir/06_plot_cancer_axis.R" "$tables_dir" "$figures_dir" "$source_data_dir" "$MIN_LOG2FC" "$MAX_FDR"
Rscript "$script_dir/06_plot_cancer_axis.R" "$tables_dir" "$figures_dir" "$source_data_dir" "$MIN_LOG2FC" "$MAX_FDR" public

# Joanito-derived outputs are git-ignored (Synapse data-use terms); mirror them
# into the DATA analysis root so they are not held only in the working tree.
restricted="$DATA_RESULTS/restricted_joanito"
mkdir -p "$restricted"/{tables,figures,source_data}
cp "$tables_dir"/{CRC_high_evidence.csv,Dataset_admission_QA.csv,Literature_31_CRC_axis_audit.csv} "$restricted/tables/"
for f in "$tables_dir"/Joanito_*; do [[ -e "$f" ]] && cp "$f" "$restricted/tables/"; done
for f in "$figures_dir"/{Joanito_*,TCGA_vs_epithelial_effect_scatter.*,Dataset_admission_QA.*}; do
  [[ -e "$f" ]] && cp "$f" "$restricted/figures/"; done
for f in "$source_data_dir"/{Joanito_*,TCGA_vs_epithelial_effect_scatter_source_data.csv,Dataset_admission_QA_source_data.csv}; do
  [[ -e "$f" ]] && cp "$f" "$restricted/source_data/"; done
printf 'Step 3 complete. Version-controlled outputs: %s/results\n' "$step_dir"
