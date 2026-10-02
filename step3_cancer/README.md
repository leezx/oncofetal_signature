# Step 3 — CRC-high cancer axis

Genome-wide tumour re-expression at two resolutions: TCGA COAD+READ bulk (T1,
population level) and CRC malignant epithelial scRNA (S1 Joanito, S2 Pelka,
epithelial-intrinsic). `CRC_high = S1 ∩ S2` passing the independent progenitor check P; TCGA is reported as `bulk_support`. Rules:
[`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md); commands:
[`docs/RUNBOOK.md`](docs/RUNBOOK.md).

**Status (2026-10-02):** all contrasts run; T1, S1 and S2 pass the
admission QA. Independent of the invalidated Step 2; the fetal ∩ CRC
intersection is deferred.

**Restricted outputs.** Joanito (Synapse syn26844071) may not be disclosed in
any derived form, so everything computed from it — Joanito DE tables,
`CRC_high_evidence.csv`, the final labels, the S1 admission-QA rows, the
literature audit, and the Joanito/bulk-vs-epithelial figures — is git-ignored.
These files are written locally by `run_step3.sh` and mirrored to
`$DATA_RESULTS/restricted_joanito/`.

## Key outputs

Committed (no Joanito content):

- `results/tables/TCGA_tumor_vs_normal_DEG.csv`, `Pelka_*_DEG.csv` — full DE tables
- `results/tables/Dataset_admission_QA_TCGA_Pelka.csv` — QA for T1, S2 and P
- `results/figures/TCGA_*`, `Pelka_*` — volcanoes, TCGA-vs-Pelka scatter, QA panel

Restricted (git-ignored; local and DATA only): `CRC_high_evidence.csv`,
`Dataset_admission_QA.csv`, `Literature_31_CRC_axis_audit.csv`, `Joanito_*`,
and the Joanito volcano, bulk-vs-Joanito scatter and full QA figures.

Positive log2FC always means higher in tumour.
