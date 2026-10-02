#!/usr/bin/env bash
# Fetal epithelial state analysis (plan v1 + addenda v1.1-1.2).
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step_dir="$(cd "$script_dir/.." && pwd)"
repo="$(cd "$step_dir/.." && pwd)"
DATA_ROOT="${FETAL_STATE_DATA_ROOT:-/Volumes/Stelligen_SSD/Stelligen/DATA}"
work="$DATA_ROOT/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step2_fetal_state_v0.1/cells"
fawkner="$DATA_ROOT/scRNAseq/GSE158702_FawknerCorbett2021/raw"
tables="$step_dir/results/tables"
mkdir -p "$work" "$tables" "$step_dir/results/figures" "$step_dir/results/source_data"
python3 "$script_dir/01_build_cells.py" --repo "$repo" --fawkner-raw "$fawkner" \
  --gao-raw "$DATA_ROOT/scRNAseq/GSE103239_Gao2018/raw" --work-dir "$work" --tables-dir "$tables"
python3 "$script_dir/02_assign_states.py" --repo "$repo" \
  --markers-xlsx "$fawkner/Fawkner2021_Mendeley_supplementary.xlsx" --work-dir "$work" --tables-dir "$tables"
python3 "$script_dir/03_tabulate_genes.py" --repo "$repo" --work-dir "$work" --tables-dir "$tables"
Rscript "$script_dir/04_plot_fetal_state.R" "$tables" "$step_dir/results/figures" "$step_dir/results/source_data"
