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

## Step 3 — CRC-high cancer axis

| Module | Purpose | Primary outputs |
|---|---|---|
| `step3_cancer/scripts/00_synapse_download.py` | Joanito Synapse files via REST, MD5-checked | raw Joanito files |
| `step3_cancer/scripts/01_recount3_tcga_gtex.R` | recount3 TCGA COAD/READ + GTEx COLON counts | combined RSE + coldata |
| `step3_cancer/scripts/02_pelka_pseudobulk.py` | Pelka patient × group epithelial pseudobulks | counts + inclusion table |
| `step3_cancer/scripts/03_tcga_tumor_vs_normal.R` | T1 plus paired and GTEx sensitivity | TCGA DEG table |
| `step3_cancer/scripts/04_pseudobulk_edger.R` | generic pseudobulk edgeR QL | S1/S2/control DEG tables |
| `step3_cancer/scripts/05_crc_high_integration.py` | admission QA, evidence table, labels, 31-gene audit | `CRC_high_evidence.csv` |
| `step3_cancer/scripts/06_plot_cancer_axis.R` | volcanoes, bulk-vs-epithelial scatter, QA panel | figures + source data |
| `step3_cancer/scripts/run_step3.sh` | run the full step | complete package |

`07_joanito_pseudobulk.py` and `config/joanito_label_map.tsv` are written once
the Joanito metadata is accessible. See `step3_cancer/docs/RUNBOOK.md`.

## Cross-step marker summary

| Module | Purpose | Primary outputs |
|---|---|---|
| `marker_summary/scripts/build_marker_summary.py` | log2FC / P / FDR of the 31 Step 1 markers in every Step 2 and Step 3 dataset (descriptive) | `marker_summary/results/Literature_31_marker_summary_statistics.{csv,xlsx}`, measurement-status table |
