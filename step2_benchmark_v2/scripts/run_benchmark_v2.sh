#!/usr/bin/env bash
# Step 2 human benchmark v2: acquisition, pseudobulks, contrasts, 31-marker
# matrix, frozen qualification verdicts and heatmap. No signature is built.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step_dir="$(cd "$script_dir/.." && pwd)"
repo="$(cd "$step_dir/.." && pwd)"
DATA_ROOT="${STEP2V2_DATA_ROOT:-/Volumes/Stelligen_SSD/Stelligen/DATA}"
work="$DATA_ROOT/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step2_human_benchmark_v0.2"
tables="$step_dir/results/tables"
edger="$repo/step3_cancer/scripts/04_pseudobulk_edger.R"
mkdir -p "$work"/{pseudobulk,de} "$tables" "$step_dir/results/figures" "$step_dir/results/source_data"

bash "$script_dir/00_download_inputs.sh"

# H-new1: Fawkner Visium (adult arm has 2 sections: minimum group size 2).
python3 "$script_dir/01_fawkner_spatial_pseudobulk.py" \
  --raw-dir "$DATA_ROOT/SpatialTranscriptomics/GSE158328_FawknerCorbett2021/raw" \
  --work-dir "$work/pseudobulk" --out-dir "$tables"
for v in units units_colon_only; do
  Rscript "$edger" "$work/pseudobulk/Fawkner_spatial_epithelial_pseudobulk_counts.csv.gz" \
    "$work/pseudobulk/Fawkner_spatial_$v.csv" Fetal Adult none "$work/de/Hnew1_Fawkner_spatial_$v.csv" 2
done

# H-new2: Fawkner fetal scRNA (hashtag samples) vs Burclaff adult donors.
python3 "$script_dir/02_fawkner_scrna_pseudobulk.py" \
  --raw-dir "$DATA_ROOT/scRNAseq/GSE158702_FawknerCorbett2021/raw" --work-dir "$work/pseudobulk"
h5ad="$DATA_ROOT/scRNAseq/GSE185224_Burclaff2022/processed/v0.1/burclaff_annotated.h5ad"
[[ -s "$h5ad" ]] || gunzip -c "$DATA_ROOT/scRNAseq/GSE185224_Burclaff2022/raw/GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad.gz" > "$h5ad"
python3 "$script_dir/03_burclaff_adult_pseudobulk.py" --h5ad "$h5ad" \
  --fetal-counts "$work/pseudobulk/Fawkner_scRNA_fetal_pseudobulk_counts.csv.gz" \
  --fetal-units "$work/pseudobulk/Fawkner_scRNA_fetal_units.csv" --work-dir "$work/pseudobulk"
Rscript "$edger" "$work/pseudobulk/Hnew2_Fawkner_fetal_Burclaff_adult_counts.csv.gz" \
  "$work/pseudobulk/Hnew2_units.csv" Fetal Adult none "$work/de/Hnew2_Fawkner_Burclaff.csv"

# H-new3, H-new3-debug, Gao-original (recomputed): cross-platform effect-only.
python3 "$script_dir/04_gao_cross_platform_effect.py" \
  --gao-dir "$DATA_ROOT/scRNAseq/GSE103239_Gao2018/raw" --wang-dir "$DATA_ROOT/scRNAseq/GSE125970_Wang2020/raw" \
  --candidates "$repo/step2_fetal/config/literature_candidates_31.tsv" \
  --aliases "$repo/step3_cancer/config/symbol_aliases.tsv" --de-dir "$work/de" --tables-dir "$tables"

# H-bulk1 / H-bulk2 (addendum v1.1): Roadmap fetal SI vs HPA adult SI / duodenum, recount3.
Rscript "$script_dir/07_bulk_primary_tissue.R" "$DATA_ROOT/bulkRNAseq/recount3_Roadmap_HPA_intestine" \
  "$work/pseudobulk" "$tables" "$edger"
for f in Hbulk1_Roadmap_vs_HPA_SI Hbulk2_Roadmap_vs_HPA_duodenum; do
  cp "$work/pseudobulk/$f.csv" "$work/de/$f.csv"; cp "$work/pseudobulk/${f}_model.csv" "$work/de/${f}_model.csv"
  cp "$work/de/$f.csv" "$tables/${f}_DE.csv"; cp "$work/de/${f}_model.csv" "$tables/${f}_DE_model.csv"
done

# Addendum v1.2: stage-resolved Gao and HGCA (breakpoint 9 weeks).
python3 "$script_dir/08_gao_stage_resolved.py" \
  --gao-dir "$DATA_ROOT/scRNAseq/GSE103239_Gao2018/raw" --wang-dir "$DATA_ROOT/scRNAseq/GSE125970_Wang2020/raw" \
  --de-dir "$work/de" --tables-dir "$tables"
Rscript "$script_dir/09_hgca_stage_resolved.R" "$DATA_ROOT/scRNAseq/HGCA_Elmentaite2021/processed/v0.1" "$work/de" "$tables"
for f in HGCA_early_fetal_vs_adult HGCA_late_fetal_vs_adult HGCA_late_vs_early_fetal HGCA_all_fetal_vs_adult \
         GaoOriginal_early_GaoLI_vs_GSE103154 GaoOriginal_late_GaoLI_vs_GSE103154 GaoLI_late_vs_early \
         Hnew3_early_GaoSILI_vs_Wang Hnew3_late_GaoSILI_vs_Wang GaoSILI_late_vs_early; do
  cp "$work/de/$f.csv" "$tables/${f}_DE.csv"
done

# Version-controlled copies.
cp "$work/de/GaoOriginal_GaoLI_vs_GSE103154.csv" "$tables/GaoOriginal_GaoLI_vs_GSE103154_DE.csv"
cp "$work/de/Hnew3_GaoSILI_vs_Wang.csv" "$tables/Hnew3_GaoSILI_vs_Wang_DE.csv"
cp "$work/de/Hnew3debug_GaoLI_vs_WangColonRectum.csv" "$tables/Hnew3debug_GaoLI_vs_WangColonRectum_DE.csv"
cp "$work/de/Hnew1_Fawkner_spatial_units.csv" "$tables/Hnew1_Fawkner_spatial_DE.csv"
cp "$work/de/Hnew1_Fawkner_spatial_units_model.csv" "$tables/Hnew1_Fawkner_spatial_DE_model.csv"
cp "$work/de/Hnew1_Fawkner_spatial_units_colon_only.csv" "$tables/Hnew1_Fawkner_spatial_colon_only_DE.csv"
cp "$work/de/Hnew1_Fawkner_spatial_units_colon_only_model.csv" "$tables/Hnew1_Fawkner_spatial_colon_only_DE_model.csv"
cp "$work/de/Hnew2_Fawkner_Burclaff.csv" "$tables/Hnew2_Fawkner_Burclaff_DE.csv"
cp "$work/de/Hnew2_Fawkner_Burclaff_model.csv" "$tables/Hnew2_Fawkner_Burclaff_DE_model.csv"
cp "$work/pseudobulk/Hnew2_units.csv" "$tables/Hnew2_units.csv"
cp "$work/pseudobulk/Fawkner_spatial_units.csv" "$tables/Hnew1_Fawkner_spatial_units.csv"
rm -f "$tables/Fawkner_spatial_units.csv"

python3 "$script_dir/05_benchmark_matrix.py" --repo "$repo" --de-dir "$work/de" --out-dir "$tables"
Rscript "$script_dir/06_plot_benchmark_heatmap.R" "$tables" "$step_dir/results/figures" "$step_dir/results/source_data"
echo "Benchmark v2 complete: $step_dir/results"
