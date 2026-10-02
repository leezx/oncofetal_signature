# Step 2 dataset qualification benchmark

This package evaluates 31 literature-curated intestinal oncofetal markers before
any replacement developmental dataset is allowed to define a hard gate.

Raw matrices are stored only under:

`/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/dataset_qualification_v0.1`

The benchmark deliberately does not construct a genome-wide signature.

## Qualification contrasts

- Human Senger, GSE101531: six fetal versus three adult epithelial
  enterosphere samples; limma on `log2(RPKM + 1)`.
- Human Fawkner-Corbett, GSE158702: fetal EpCAM-positive scRNA-seq only. No
  matched adult arm exists, so this dataset reports fetal expression and
  detection, not fetal/adult fold change.
- Mouse Pikkupeura, GSE160449: fetal versus adult epithelial cultures are
  evaluated separately in the LN and collagen arms with edgeR; their median
  log2FC and directional agreement are reported.
- Mouse GSE44433: WT fetal E17.5 versus WT adult eight-week laser-microdissected
  distal ileal epithelium, five versus five; Tnf mutant and maternal-genotype
  arms are excluded.

Run:

```bash
bash step2_dataset_benchmark/scripts/00_download_benchmark_inputs.sh
Rscript step2_dataset_benchmark/scripts/01_benchmark_31_markers.R
```
