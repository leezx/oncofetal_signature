# Step 2 human benchmark v2

Dataset qualification for a replacement human developmental contrast. Three new
human fetal-vs-adult contrasts (H-new1 Visium, H-new2 Fawkner + Burclaff,
H-new3 Gao + Wang, plus a Gao debug pair) are benchmarked on the 31 Step 1
markers alongside HGCA, Gao-original, Senger and Pikkupeura. **No signature is
built.**

- Frozen rules: [`docs/QUALIFICATION_PLAN.md`](docs/QUALIFICATION_PLAN.md)
- Results and interpretation: [`docs/BENCHMARK_V2_REPORT.md`](docs/BENCHMARK_V2_REPORT.md)
- Run everything: `bash step2_benchmark_v2/scripts/run_benchmark_v2.sh`

**Outcome:** only H-new2 (Fawkner fetal EPCAM+ scRNA vs Burclaff adult
epithelium) qualifies among human contrasts (TNFRSF12A +1.83, P < 0.001;
18/24 markers fetal-positive). Visium is composition-confounded; the Gao debug
pair places the Gao failure in the fetal arm, not the adult reference.

| Path | Content |
|---|---|
| `scripts/00–06`, `run_benchmark_v2.sh` | acquisition, pseudobulks, contrasts, matrix, heatmap |
| `results/tables/Benchmark_v2_31_marker_matrix.csv` | 31 markers × 10 contrasts: log2FC / P / FDR |
| `results/tables/Benchmark_v2_qualification.csv` | rule 1, rule 2, verdict per contrast |
| `results/tables/*_DE.csv` | full gene tables of the new contrasts |
| `results/tables/*_units.csv`, `*_marker_unit_values.csv` | per-unit inclusion and marker values |
| `results/figures/Benchmark_v2_31_marker_heatmap.*` | heatmap (PDF + 600-dpi PNG) |

Intermediates: `DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step2_human_benchmark_v0.2/`.
