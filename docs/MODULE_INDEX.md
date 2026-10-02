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

## Step 2 human benchmark v2

| Module | Purpose | Primary outputs |
|---|---|---|
| `step2_benchmark_v2/scripts/00_download_inputs.sh` | GSE158328, GSE158702 (+HTO), GSE185224, GSE125970 into DATA | raw files + checksums |
| `step2_benchmark_v2/scripts/01_fawkner_spatial_pseudobulk.py` | H-new1 Visium epithelial-spot pseudobulk per section | counts + section table |
| `step2_benchmark_v2/scripts/02_fawkner_scrna_pseudobulk.py` | H-new2 fetal: hashtag demultiplexing, per-sample pseudobulk | counts + unit table |
| `step2_benchmark_v2/scripts/03_burclaff_adult_pseudobulk.py` | H-new2 adult: Burclaff donor pseudobulks, joined to fetal | merged counts |
| `step2_benchmark_v2/scripts/04_gao_cross_platform_effect.py` | H-new3, H-new3-debug, Gao-original recomputed (effect-only) | DE + per-unit marker tables |
| `step2_benchmark_v2/scripts/07_bulk_primary_tissue.R` | H-bulk1/H-bulk2: recount3 Roadmap fetal SI vs HPA adult SI/duodenum, edgeR | DE tables + unit table |
| `step2_benchmark_v2/scripts/05_benchmark_matrix.py` | 31-marker × contrast matrix and frozen verdicts | matrix + qualification |
| `step2_benchmark_v2/scripts/06_plot_benchmark_heatmap.R` | heatmap | PDF/PNG + source data |
| `step2_benchmark_v2/scripts/run_benchmark_v2.sh` | run all | complete package |

edgeR contrasts reuse `step3_cancer/scripts/04_pseudobulk_edger.R` (now with an
optional minimum-units argument; default 3, H-new1 uses 2).

## Step 2 fetal epithelial state (manual)

| Module | Purpose | Primary outputs |
|---|---|---|
| `step2_fetal_state/scripts/01_build_cells.py` | Fawkner hashtag → sample (Mendeley key), sex/region/discordance QC; Gao SI/LI epithelial cells | cell h5ad (DATA) + hashtag check table |
| `step2_fetal_state/scripts/02_assign_states.py` | Leiden + author-marker state labels (analysed genes excluded), confidence rule | frozen labels + label checks |
| `step2_fetal_state/scripts/03_tabulate_genes.py` | unit and unit × state values, age correlations, composition decomposition, region differences | tables |
| `step2_fetal_state/scripts/04_plot_fetal_state.R` | unit heatmap; gene by state × age | PDF/PNG + source data |
| `step2_fetal_state/scripts/run_fetal_state.sh` | run all | complete package |

