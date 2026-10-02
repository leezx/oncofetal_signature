#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step_dir="$(cd "$script_dir/.." && pwd)"
source "$step_dir/config/step2.env"

tables_dir="$step_dir/results/tables"
figures_dir="$step_dir/results/figures"
source_data_dir="$step_dir/results/source_data"
mkdir -p "$HGCA_PROCESSED" "$DATA_RESULTS" "$tables_dir" "$figures_dir" "$source_data_dir"

python3 "$script_dir/01_hgca_pseudobulk.py" \
  --h5ad "$HGCA_H5AD" \
  --output-dir "$HGCA_PROCESSED" \
  --min-cells "$MIN_CELLS"

Rscript "$script_dir/02_hgca_fetal_high.R" \
  "$HGCA_PROCESSED" "$DATA_RESULTS" "$MIN_LOG2FC"

python3 "$script_dir/03_gao_validation.py" \
  --data-dir "$GAO_RAW" \
  --output-dir "$DATA_RESULTS"

Rscript "$script_dir/04_mouse_fetal_high.R" \
  "$MOUSE_COUNTS" "$DATA_RESULTS"

python3 "$script_dir/05_cross_species_intersection.py" \
  --tables-dir "$DATA_RESULTS" \
  --orthologues "$ORTHOLOGUES"

cp "$DATA_RESULTS/HGCA_DE.csv" "$tables_dir/human_HGCA_fetal_vs_adult_DEG.csv"
cp "$DATA_RESULTS/HGCA_proliferative_control.csv" "$tables_dir/human_HGCA_proliferative_control_DEG.csv"
cp "$DATA_RESULTS/Gao_validation.csv" "$tables_dir/human_Gao_validation.csv"
cp "$DATA_RESULTS/Mouse_DE.csv" "$tables_dir/mouse_in_vivo_fetal_vs_adult_DEG.csv"
cp "$DATA_RESULTS/Mouse_symbol_duplicate_check.csv" "$tables_dir/Mouse_symbol_duplicate_check.csv"
cp "$DATA_RESULTS/Mouse_log2FC_gt5_raw_CPM.csv" "$tables_dir/Mouse_log2FC_gt5_raw_CPM.csv"
cp "$DATA_RESULTS/Cross_species_evidence.csv" "$tables_dir/Cross_species_evidence.csv"
cp "$DATA_RESULTS/Conserved_Fetal_High.csv" "$tables_dir/Conserved_Fetal_High.csv"

Rscript "$script_dir/06_plot_volcano.R" \
  "$tables_dir" "$figures_dir" "$source_data_dir"

Rscript "$script_dir/07_plot_human_mouse_effect_scatter.R" \
  "$tables_dir/Conserved_Fetal_High.csv" "$figures_dir" "$source_data_dir"

printf 'Step 2 complete. Version-controlled outputs: %s/results\n' "$step_dir"
