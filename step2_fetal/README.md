# Step 2 fetal differential expression

This directory contains the requested human and mouse fetal-versus-adult DEG
tables and volcano plots. Raw and processed expression matrices remain outside
GitHub under the managed `DATA` hierarchy.

## DEG tables

- `human_HGCA_fetal_vs_adult_DEG.csv`: HGCA donor-pseudobulk edgeR results.
  The primary fetal-high rule is `HGCA_log2FC >= 0.5` and `HGCA_FDR < 0.05`.
  The table also retains the adult proliferative-control result and final `H1_pass`.
- `mouse_in_vivo_fetal_vs_adult_DEG.csv`: mouse in-vivo fetal epithelium versus
  adult crypt epithelium edgeR results. `Mouse_pass` uses `Mouse_log2FC >= 0.5`
  and `Mouse_FDR < 0.05`.

Positive log2 fold change means higher expression in fetal epithelium. The full
tested-gene tables are included, not only the threshold-passing subset.

## Volcano plots

- `human_HGCA_fetal_vs_adult_volcano.{pdf,png}`
- `mouse_in_vivo_fetal_vs_adult_volcano.{pdf,png}`
- `human_mouse_fetal_vs_adult_volcano.{pdf,png}`: two-panel comparison.

Red marks fetal-high genes, blue marks adult-high genes, and grey marks genes
that do not meet both `|log2FC| >= 0.5` and `FDR < 0.05`. Labels show the ten
most statistically significant genes in each direction. PDF is the vector
version; PNG is a 600-dpi review preview.

`human_volcano_source_data.csv` and `mouse_volcano_source_data.csv` provide the
exact plotting columns. Reproduce all figures with:

```bash
Rscript step2_fetal/plot_volcano.R step2_fetal step2_fetal
```
