# Step 3 runbook

## Purpose

Reproduce the cancer axis from frozen public inputs through the gene-level
`CRC_high_evidence.csv`, admission QA, literature audit, and figures. Rules are
defined in [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md).

## Requirements

- R ≥ 4.5 with `recount3`, `SummarizedExperiment`, `edgeR` (≥ 4), `ggplot2`,
  `ggrepel`, `dplyr`, `ragg`
- Python 3 with `h5py`, `numpy`, `pandas`, `scipy`
- Synapse personal access token in `~/.synapseConfig` (`authtoken = ...`) or
  `SYNAPSE_AUTH_TOKEN`, **and** the syn26844071 self-sign terms of use accepted
  once in the Synapse web UI by the account owner

## Run

From the repository root:

```bash
bash step3_cancer/scripts/run_step3.sh
```

| Module | Role |
|---|---|
| `00_synapse_download.py` | Joanito files via Synapse REST, MD5-checked, idempotent |
| `01_recount3_tcga_gtex.R` | recount3 TCGA COAD/READ + GTEx COLON gene counts |
| `02_pelka_pseudobulk.py` | Pelka patient × group epithelial pseudobulks |
| `03_tcga_tumor_vs_normal.R` | T1 plus paired and GTEx sensitivity |
| `04_pseudobulk_edger.R` | generic pseudobulk edgeR QL (S1, S2, proliferation controls) |
| `05_crc_high_integration.py` | admission QA, evidence table, final labels, 31-gene audit |
| `06_plot_cancer_axis.R` | volcanoes, bulk-vs-epithelial scatter, QA panel |
| `07_joanito_pseudobulk.py` | Joanito pseudobulks (written once metadata is accessible) |

Large intermediates go to `DATA` (see `config/step3.env`); only tables,
source data, and figures are copied into `step3_cancer/results/`.

## Frozen gates

| Gate | Rule |
|---|---|
| Admission QA | ≥ 80% of measured tumour-up and tumour-down panel genes in the expected direction |
| T1 TCGA | log2FC ≥ 0.5, FDR < 0.05, `~ tissue + project` |
| S1 Joanito discovery | log2FC ≥ 0.5, FDR < 0.05, `~ group + cohort` |
| S1 proliferation control | malignant / normal stem-TA-cycling log2FC > 0 |
| S2 Pelka | log2FC > 0, FDR < 0.05, `~ group` |
| Final | `CRC_high` = S1 ∩ S2; `bulk_support` = T1 (no veto) |
