# Step 2 fetal differential expression

This is the complete, reproducible Step 2 analysis package. Raw and processed
expression matrices remain outside GitHub under the managed `DATA` hierarchy.

## Directory layout

```text
config/              frozen paths, thresholds, and Ensembl BioMart query
docs/                analysis plan, runbook, execution report, QA, HTML report
scripts/             00–06 ordered workflow
results/tables/      DEG and final cross-species tables
results/source_data/ exact volcano plotting data
results/figures/     PDF and 600-dpi PNG plots
```

## DEG tables

- `results/tables/human_HGCA_fetal_vs_adult_DEG.csv`: HGCA donor-pseudobulk edgeR results.
  The primary fetal-high rule is `HGCA_log2FC >= 0.5` and `HGCA_FDR < 0.05`.
  The table also retains the adult proliferative-control result and final `H1_pass`.
- `results/tables/mouse_in_vivo_fetal_vs_adult_DEG.csv`: mouse in-vivo fetal epithelium versus
  adult crypt epithelium edgeR results. `Mouse_pass` uses `Mouse_log2FC >= 0.5`
  and `Mouse_FDR < 0.05`.

Positive log2 fold change means higher expression in fetal epithelium. The full
tested-gene tables are included, not only the threshold-passing subset.

## Volcano plots

- `results/figures/human_HGCA_fetal_vs_adult_volcano.{pdf,png}`
- `results/figures/mouse_in_vivo_fetal_vs_adult_volcano.{pdf,png}`
- `results/figures/human_mouse_fetal_vs_adult_volcano.{pdf,png}`

Red marks fetal-high genes, blue marks adult-high genes, and grey marks genes
that do not meet both `|log2FC| >= 0.5` and `FDR < 0.05`. Labels show the ten
most statistically significant genes in each direction. PDF is the vector
version; PNG is a 600-dpi review preview.

`results/source_data/human_volcano_source_data.csv` and
`results/source_data/mouse_volcano_source_data.csv` provide the
exact plotting columns. Reproduce all figures with:

```bash
bash step2_fetal/scripts/run_step2.sh
```

For commands, inputs, outputs, and restart points, see
[`docs/RUNBOOK.md`](docs/RUNBOOK.md). The rendered checkpoint report is
[`docs/STEP2_REPORT.html`](docs/STEP2_REPORT.html).
