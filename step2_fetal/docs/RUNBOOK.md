# Step 2 runbook

## Purpose

Reproduce the developmental-axis analysis from frozen source files through the
human/mouse DEG tables, volcano plots, and 706-gene conserved fetal-high set.

## Requirements

- Python 3 with `anndata`, `numpy`, `pandas`, and `scipy`
- R with `edgeR`, `ggplot2`, `ggrepel`, `dplyr`, `patchwork`, and `ragg`
- `curl` for optional acquisition
- writable `/Volumes/Stelligen_SSD/Stelligen/DATA`, or set `STEP2_DATA_ROOT`

## Data acquisition

Raw data never enter Git. To download any missing public inputs and copy the
pre-existing mouse count matrix into its managed dataset directory:

```bash
bash step2_fetal/scripts/00_download_inputs.sh
```

The script is idempotent: an existing non-empty source file is not overwritten.
Source URLs, the frozen Ensembl release 116 XML query, and DATA destinations are
version controlled.

## Complete analysis

From the repository root:

```bash
bash step2_fetal/scripts/run_step2.sh
```

The wrapper reads `config/step2.env` and runs:

1. HGCA metadata filtering and donor-level raw-count pseudobulk.
2. HGCA edgeR fetal/adult and fetal/adult-proliferative contrasts.
3. Gao large-intestinal effect-only replication.
4. Mouse in-vivo edgeR fetal/adult contrast.
5. Ensembl one-to-one mapping and strict H1 ∩ H2 ∩ M1 intersection.
6. Human and mouse volcano plots.

Large intermediates are written to DATA. Stable DEG tables and figures are
copied into `step2_fetal/results/`.

## Restart points

Each numbered script is independently runnable. Its CLI is printed when invoked
without required arguments. If only figure styling changes, rerun module 06:

```bash
Rscript step2_fetal/scripts/06_plot_volcano.R \
  step2_fetal/results/tables \
  step2_fetal/results/figures \
  step2_fetal/results/source_data
```

## Frozen gates

| Gate | Rule |
|---|---|
| HGCA primary | fetal/adult log2FC ≥ 0.5 and FDR < 0.05 |
| HGCA proliferative control | fetal/adult-proliferative log2FC > 0 |
| Gao H2 | log2FC ≥ 0.5 and fetal median above each adult donor |
| Mouse M1 | fetal/adult log2FC ≥ 0.5 and FDR < 0.05 |
| Orthology | Ensembl release 116 `ortholog_one2one` |
| Final | H1 ∩ H2 ∩ M1 |

Do not optimize these thresholds using CRC or TWEAKR results.
