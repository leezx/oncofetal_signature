# Step 3 — CRC-high cancer axis

Genome-wide tumour re-expression at two resolutions: TCGA COAD+READ bulk (T1,
population level) and CRC malignant epithelial scRNA (S1 Joanito, S2 Pelka,
epithelial-intrinsic). `CRC_high = T1 ∩ S1 ∩ S2`. Rules:
[`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md); commands:
[`docs/RUNBOOK.md`](docs/RUNBOOK.md).

**Status (2026-10-02):** T1 and S2 complete and admitted by QA; S1 pending
Synapse access, so no `CRC_high` call exists yet. Independent of the
invalidated Step 2; the fetal ∩ CRC intersection is deferred.

## Key outputs

- `results/tables/CRC_high_evidence.csv` — one row per gene, all effects, gates, labels
- `results/tables/Dataset_admission_QA.csv` — canonical-marker direction check
- `results/tables/Literature_31_CRC_axis_audit.csv` — descriptive candidate audit
- `results/tables/TCGA_tumor_vs_normal_DEG.csv`, `Pelka_*_DEG.csv` — full DE tables
- `results/figures/` — volcanoes, bulk-vs-epithelial scatter, QA panel (PDF + 600-dpi PNG)

Positive log2FC always means higher in tumour.
