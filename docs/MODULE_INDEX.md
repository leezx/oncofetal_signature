# Module index

## Step 2 — developmental axis

| Module | Purpose | Primary outputs |
|---|---|---|
| `step2_fetal/scripts/00_download_inputs.sh` | Acquire public inputs into DATA | raw dataset files and BioMart response |
| `step2_fetal/scripts/01_hgca_pseudobulk.py` | Inspect HGCA metadata and aggregate donor raw counts | HGCA pseudobulks and manifests |
| `step2_fetal/scripts/02_hgca_fetal_high.R` | Run HGCA primary and proliferative edgeR contrasts | human HGCA DEG tables |
| `step2_fetal/scripts/03_gao_validation.py` | Calculate Gao donor-level effect-only replication | Gao validation table |
| `step2_fetal/scripts/04_mouse_fetal_high.R` | Run in-vivo mouse edgeR contrast | mouse DEG table |
| `step2_fetal/scripts/05_cross_species_intersection.py` | Apply one-to-one mapping and H1 ∩ H2 ∩ M1 | evidence and 706-gene set |
| `step2_fetal/scripts/06_plot_volcano.R` | Produce volcano plots and source data | PDF/PNG figures and CSVs |
| `step2_fetal/scripts/07_plot_human_mouse_effect_scatter.R` | Plot cross-species effect sizes for the final 706 genes | effect-size scatter and marker-status table |
| `step2_fetal/scripts/08_conservation_and_candidate_audit.R` | Audit 31 literature candidates and quantify conservation across all genes tested in both species | gate audit, all-tested scatter, 2×2 overlap, Spearman and Fisher statistics |
| `step2_fetal/scripts/run_step2.sh` | Execute modules 01–08 in order | complete Step 2 package |

Configuration is centralized in `step2_fetal/config/step2.env`. See
`step2_fetal/docs/RUNBOOK.md` for commands and restart points.

**Validity status:** the modules remain reproducible, but the current Step 2
dataset/contrast failed literature-marker biological QA. Do not use the 706-gene
output downstream; see `docs/MAJOR_REVISION_LOG.md`.

## Replacement dataset qualification

| Module | Purpose | Primary outputs |
|---|---|---|
| `step2_dataset_benchmark/scripts/00_download_benchmark_inputs.sh` | Download four benchmark datasets and NCBI mapping into DATA | raw inputs and checksums |
| `step2_dataset_benchmark/scripts/01_benchmark_31_markers.R` | Evaluate 31 prespecified markers without constructing a signature | four-dataset benchmark, metadata audit, priority-five table |

See `step2_dataset_benchmark/docs/BENCHMARK_REPORT.md` for the review decision.
