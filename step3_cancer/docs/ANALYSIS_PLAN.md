# Step 3 — CRC-High Cancer Axis

## Material Passport

- Origin: user task brief "Step 3 — Cancer axis: TCGA + malignant epithelial scRNA" (2026-10-02)
- Version label: step3_plan_v1
- Status: frozen before any tumour/normal expression result was examined.
  Only sample-level metadata (specimen types, cluster labels, patient counts)
  was inspected while writing this plan.

## Objective

Define, genome-wide and independently of the developmental axis, which genes
are re-expressed in colorectal cancer at two levels of resolution:

```text
Level 1  TCGA COAD+READ primary tumour vs normal colon      → population-level tumour re-expression (T1)
Level 2  CRC malignant epithelium vs normal epithelium (sc) → epithelial-intrinsic re-expression     (S1, S2)

CRC-High = T1 ∩ S1 ∩ S2
```

A gene is called **CRC-high** only when both levels pass. TCGA alone cannot
distinguish epithelial re-expression from fibroblast, myeloid, or endothelial
contributions to bulk tumour RNA; scRNA alone has fewer patients.

## Relationship to Step 2 (decision recorded 2026-10-02)

Step 2 is invalidated pending dataset review (`docs/MAJOR_REVISION_LOG.md`).
At the user's explicit instruction Step 3 starts now, under these constraints:

- Step 3 is computed **genome-wide** and does **not** use the provisional
  706-gene Step 2 output as a filter, prior, or universe.
- The fetal ∩ CRC intersection (Fig 1D) is **deferred** until a replacement
  Step 2 run passes biological QA. Step 3 outputs must not be used to choose
  or tune the replacement Step 2 datasets, contrasts, or thresholds.
- The 31 literature candidates are reported descriptively (audit table); they
  are never used to set thresholds or to admit datasets.

## Dataset roles

| Role | Dataset | Comparison | Gate |
|---|---|---|---|
| Bulk, population level | TCGA COAD + READ via recount3 (Monorail, GENCODE v26) | primary tumour (01) vs solid-tissue normal (11) | **T1** |
| Bulk sensitivity | GTEx v8 COLON via recount3 (same pipeline) | TCGA primary tumour vs GTEx transverse colon | descriptive only |
| Bulk sensitivity | TCGA paired subset | tumour vs matched normal, same patient | descriptive only |
| scRNA primary | Joanito et al. 2022 Nat Genet, Synapse syn26844071 (5 cohorts) | author-called malignant epithelial vs normal epithelial | **S1** |
| scRNA replication | Pelka et al. 2021 Cell, GSE178341 | tumour-specimen epithelium vs normal-specimen epithelium | **S2** |
| Proliferation control | normal stem/TA/cycling epithelium (see below) | malignant vs normal proliferative epithelium | directional, part of S1 |

Joanito and Pelka do not share patients or cohorts. Deferred (not in this
step): the in-house CRC Atlas (large-scale replication / gene–gene network,
Fig 1F) and the FAP → adenoma → CRC progression data (disease-evolution
validation). Neither is used for discovery.

## Dataset admission QA (prespecified, oncofetal-independent)

Lesson from Step 2: a reproducible contrast is not necessarily a valid one.
Before any gate is interpreted, every contrast must recover the expected
direction for canonical CRC-vs-normal colon markers **that are not oncofetal
candidates**:

- Tumour-up panel: `CDH3`, `CLDN1`, `FOXQ1`, `KRT23`, `LGR5`, `ASCL2`, `MYC`, `TESC`
- Tumour-down panel: `CA1`, `CA2`, `CA4`, `GUCA2A`, `GUCA2B`, `CLCA4`, `AQP8`,
  `SLC26A3`, `MS4A12`, `CEACAM7`

Admission rule per contrast (T1, S1, S2): among measured panel genes, at least
80% of the tumour-up panel has log2FC > 0 **and** at least 80% of the
tumour-down panel has log2FC < 0. A contrast that fails is stopped and
reviewed; thresholds are not changed to rescue it.

## Level 1 — T1, TCGA COAD + READ

### Samples

- Source: recount3 `tcga` projects `COAD` and `READ`, gene-level read counts
  (`transform_counts`), GENCODE v26.
- Tumour: GDC sample type `Primary Tumor` (barcode `-01`).
- Normal: `Solid Tissue Normal` (barcode `-11`).
- Exclude: FFPE samples; metastatic/recurrent samples; any other sample type.
- One sample per patient per group: if a patient has several eligible
  aliquots in a group, keep the one with the largest total counts.

### Model

edgeR quasi-likelihood on raw counts, `filterByExpr(group = tissue)`, TMM:

```text
~ tissue + project        (tissue: Normal / Tumor; project: COAD / READ)
```

T1 passes when tumour/normal **log2FC ≥ 0.5** and **BH FDR < 0.05**.

### Sensitivity analyses (descriptive, never gates)

1. **Paired**: patients with both a tumour and a normal sample;
   `~ patient + tissue`.
2. **GTEx**: TCGA primary tumour vs GTEx `Colon - Transverse` (mucosa-bearing),
   identical recount3 processing; `~ tissue`. Because source and tissue are
   perfectly confounded, this contrast is reported as directional concordance
   only (`GTEx_concordant = sign agrees with T1`).

## Level 2 — S1, Joanito 2022 malignant epithelium

### Cells

- Use the authors' epithelial count matrix and metadata (49,155 epithelial cells).
- Malignant: cells the authors classify as tumour/malignant epithelium
  (CNV-informed; iCMS2/iCMS3 or equivalent label) from tumour samples.
- Normal: cells the authors classify as normal epithelium from
  adjacent-normal samples.
- Normal-like epithelial cells found inside tumour samples are excluded from
  both groups.
- Exact label mapping is recorded in `config/joanito_label_map.tsv` after
  metadata inspection and before differential expression.

### Pseudobulk and model

- Aggregate raw counts per **patient × group**. A patient contributes at most
  one malignant and one normal pseudobulk.
- Eligibility: ≥ 50 cells per pseudobulk (same pragmatic rule as Step 2).
- edgeR QL, `filterByExpr(group)`, TMM, design `~ group + cohort` (cohort =
  the authors' dataset of origin). If cohort is not estimable, fall back to
  `~ group` and record it.

S1 discovery passes when malignant/normal **log2FC ≥ 0.5** and **FDR < 0.05**.

### Proliferation control (part of S1)

Tumour epithelium is depleted of differentiated colonocytes and enriched in
stem/TA-like cells, so malignant vs whole-normal-epithelium partly measures
differentiation and proliferation. Compare malignant epithelium with normal
**stem/TA/cycling** epithelium only (author labels; patient-level pseudobulk,
same eligibility and model). Required: log2FC **> 0** (directional, no FDR).

If Joanito lacks normal epithelial subtype labels, the control is computed in
Pelka using normal-specimen clusters `cE01`, `cE02`, `cE03` (Stem/TA-like,
Stem/TA-like/Immature Goblet, Stem/TA-like prolif), with the same rule.

**S1 = discovery ∩ proliferation control.**

## Level 2 — S2, Pelka 2021 replication

- Epithelial cells: `clTopLevel == "Epi"`.
- Tumour group: epithelial cells from `SPECIMEN_TYPE == "T"` (62 patients).
  The GEO release has no per-cell malignancy call, so this group is
  "tumour-specimen epithelium, predominantly malignant". This is why Pelka is
  replication and not the primary dataset.
- Normal group: epithelial cells from `SPECIMEN_TYPE == "N"` (36 patients).
- Pseudobulk per patient × specimen type, ≥ 50 cells, edgeR QL, `~ group`.

S2 passes when log2FC **> 0** and **FDR < 0.05** (direction plus significance;
the discovery effect-size threshold is applied once, in S1).

## Final labels

| Label | Rule |
|---|---|
| `CRC_high` | T1 ∩ S1 ∩ S2 |
| `CRC_high_epithelial_unreplicated` | T1 ∩ S1, fails S2 |
| `Tumour_level_only` | T1 but fails S1 discovery — possible microenvironment-driven bulk signal |
| `Epithelial_only` | S1 ∩ S2 but fails T1 — epithelial re-expression diluted in bulk |
| `Proliferation_associated_reject` | passes S1 discovery and T1 but malignant/normal-proliferative log2FC ≤ 0 |
| `Not_CRC_high` | everything else tested |

Gene matching across datasets uses Ensembl gene ID (version stripped) where
both sources provide it, otherwise the HGNC symbol; the matching key is
recorded per row.

## Hard stops

- Fewer than 3 eligible patients per group in any gated contrast.
- Dataset admission QA fails for a gated contrast.
- Joanito malignant/normal labels cannot be identified from the released
  metadata without using candidate genes.
- No cell-level Wilcoxon substitute for a stopped pseudobulk analysis.

## Required outputs

`step3_cancer/results/tables/`:

- `TCGA_tumor_vs_normal_DEG.csv` (T1 + paired + GTEx sensitivity columns)
- `Joanito_malignant_vs_normal_DEG.csv` (S1 discovery + proliferation control)
- `Pelka_tumor_vs_normal_epithelium_DEG.csv` (S2)
- `CRC_high_evidence.csv` — one row per gene with every effect, FDR, gate
  flag and final label
- `Dataset_admission_QA.csv`
- `Literature_31_CRC_axis_audit.csv` (descriptive)
- sample/pseudobulk inclusion tables for each dataset

Figures: one volcano per gated contrast, a T1 vs S1 effect-size scatter over
all jointly tested genes, and the admission-QA panel; PDF + 600-dpi PNG with
plotting source data.

## Storage

- Raw/processed data: `DATA/bulkRNAseq/recount3_TCGA_COAD_READ_GTEx_COLON/`,
  `DATA/scRNAseq/Joanito2022_syn26844071/`, `DATA/scRNAseq/GSE178341/`.
- Working results: `DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/`.
- Git receives code, configs, docs, final tables, source data, and figures only.

## Sources

- Joanito I. et al. Single-cell and bulk transcriptome sequencing identifies
  two epithelial tumor cell states and refines the consensus molecular
  classification of colorectal cancer. Nat Genet 54, 963–975 (2022).
  <https://doi.org/10.1038/s41588-022-01100-4>; data Synapse syn26844071.
- Pelka K. et al. Spatially organized multicellular immune hubs in human
  colorectal cancer. Cell 184, 4734–4752 (2021); GEO GSE178341.
- Wilks C. et al. recount3: summaries and queries for large-scale RNA-seq
  expression and splicing. Genome Biol 22, 323 (2021).
