# Pan-cancer epithelial reactivation of CIOC genes: plan (FROZEN 2026-10-05, before any statistics)

## Question (narrow, by review decision)

> In which carcinomas is each gene reactivated in malignant epithelial cells
> relative to the corresponding normal epithelial compartment?

This is a supplementary **annotation**. It adds no gate and makes no
membership change, and it does not address fetal tissue specificity.

**It replaces the bulk TCGA/GTEx breadth annotation as the primary evidence.**
That annotation (breadth plan v2) is downgraded to a secondary annotation
because bulk data are composition-confounded, which matters especially for
SPP1-, CLU- and CCN2-like genes.

**Interpretation.** CRC-biased and broadly reactivated components are both
CIOC components. The CIOC defines an intestinal oncofetal state, not CRC
exclusivity.

## Data

Kang et al., *Nat Commun* 2024, tumour–normal pan-cancer single-cell atlas
(Zenodo 10.5281/zenodo.10651059, `atlas_dataset.tar`).
- **Scale:** 104 datasets, 30 cancer types, 1,070 tumour and 493 normal
  samples.
- **Matrices:** log1p(CP10k) per dataset.
- **Cell metadata:** `Dataset`, `Organ_origin`, `Sample`, `Patient`, `Tissue`
  (tumour, adjacent normal or healthy normal), `Cancer type`, `cnv_status`
  (authors' inferCNVpy call) and `Celltype`.
- **Extraction (no statistics):** `scripts/pancancer_01_extract_kang.py`
  builds mean-CP10k pseudobulks per Sample × Tissue × Cancer type × Organ ×
  Celltype × cnv_status, for the 338 + Literature-31 + calibrator genes.

## Independence from signature construction

These Kang datasets overlap data used to build the signatures and are
**excluded from every comparison**:

| Dataset | Overlap |
|---|---|
| `crc_HOLee_GSE132465` | SMC; a Joanito cohort |
| `crc_HOLee_GSE144735` | KUL3; a Joanito cohort |
| `crc_HOLee_GSE132257` | SMC group; possible overlap |
| `crc_meta_GSE178318` | Che 2021; Extended CIOC / ECOS-38 refinement |
| `intestine_gca` | HGCA; Gate H |
| `crc_pan_blueprint` | Qian 2020, Leuven; possible KUL overlap with Joanito; excluded conservatively |

The CRC column therefore uses only independent Kang CRC datasets.

## Metadata scan (counts only; no expression values inspected)

The extraction processed 103 dataset files and resolved 369 of the 371
universe genes. LY6A and REG3B have no human orthologue.

**Vocabulary:**
- **Tissue:** Tumor, Normal or Metaplasia. Adjacent and healthy normal are
  not separate labels; healthy normal = Normal samples in normal-only
  datasets.
- **Celltype:** a single "Epithelial" label covers all epithelial and
  parenchymal cells, including hepatocytes and renal epithelium.
- **cnv_status:** tumor or normal (authors' inferCNVpy call).
- **Cancer types:** "LC" pools lung cancers, so LUAD and LUSC cannot be
  separated.

## Frozen analysis

**Cells:**
- **Malignant arm:** Celltype = Epithelial, cnv_status = tumor, Tissue = Tumor.
  These are CNV-inferred malignant epithelial cells, by the authors'
  inference.
- **Normal arm:** Celltype = Epithelial, cnv_status = normal, Tissue = Normal,
  with Organ_origin matching the cancer type.
- Metaplasia is excluded, and the six construction-overlapping datasets are
  excluded (see above).

**Unit.** Patient (within dataset) × arm. Cells are pooled across that
patient's samples by a cell-weighted mean of the sample pseudobulks. A unit
needs ≥ 20 cells. The value is y = log2(mean CP10k + 1). No cell-level
tests are used.

**Carcinomas and matched normal organ:**

| Cancer type | Normal organ |
|---|---|
| CRC | Colon |
| STAD | Stomach |
| PAAD | Pancreas |
| HCC | Liver |
| CHOL | Bile Duct |
| LC | Lung |
| BRCA | Breast |
| RCC | Kidney |
| OV | Ovary + Fallopian Tube |
| THCA | Thyroid |
| PRAD | Prostate |
| BLCA | Bladder |
| HNSC | Head and Neck |
| UCEC | Uterus |
| SSCC | Skin |

Non-carcinomas are excluded: GBM, LGG, sarcomas, UVM, MEL, NHL, MM, NB,
ALL, CLL, NET and WILM.

**Tiers, per cancer type:**
- **Tier 1 (within-study):** datasets with ≥ 2 units in each arm. The model
  is OLS y ~ arm + dataset on those datasets only, with ≥ 3 units per arm
  in total. Pairing within patients is not modelled, which is conservative.
- **Tier 2 (cross-study; weaker):** used if Tier 1 is unavailable. All
  tumour units vs all normal units of the matched organ, OLS y ~ arm, with
  ≥ 3 units per arm.
- Otherwise the cancer type is **not evaluable**.

**Units available:**

| Cancer | Tier | Malignant / normal units |
|---|---|---|
| CRC | 1 (`crc_GSE166555`, independent) | 12 / 11 |
| PAAD | 1 | 12 / 3 |
| HCC | 1 | 8 / 13 |
| CHOL | 1 | 3 / 3 |
| LC | 1 | 47 / 33 |
| BRCA | 1 | 30 / 22 |
| RCC | 1 | 13 / 24 |
| OV | 1 | 6 / 6 |
| SSCC | 1 | 7 / 10 |
| STAD | 2 | 4 / 11 |
| PRAD | 2 | 5 / 6 |
| BLCA | 2 | 4 / 3 |
| HNSC | 2 (normals: gingiva, salivary, nasopharynx) | 35 / 11 |
| THCA | not evaluable | 1 normal unit |
| UCEC | not evaluable | no normal |

**Reactivation.** The arm coefficient (log2FC) must be ≥ 0.5 with BH
FDR < 0.05. BH is applied within each cancer type across the resolved gene
universe (338 + Literature-31 + calibrators). The thresholds match the CIOC
CRC gate.

**Breadth:**
- **Primary:** k / n over the 12 evaluable non-CRC carcinomas (8 Tier 1 and
  4 Tier 2). The CRC column is reported separately.
- **Sensitivity:** k / n over the 8 Tier 1 non-CRC carcinomas.
- Counts are reported, with no pan-cancer dichotomy and no
  "pan-oncofetal" wording.
- **Descriptive labels:** "CRC-biased" when CRC is reactivated and
  k / n ≤ 2 / 12; otherwise "reactivated in k/n other carcinomas". These
  are descriptive only.

**Calibrators** (reported, not used for tuning): EPCAM, KRT8, CDH1, MKI67,
TOP2A, CEACAM5, CDX2, PTPRC, COL1A2, LYZ.

**Caveats (disclosed):**
- Kang matrices are log-normalised, so a pseudobulk mean of CP10k is used
  rather than raw counts.
- CNV calls come from the authors.
- Lung cancers are pooled.
- Tier 2 comparisons are cross-study.
- Several types have small n (CHOL 3/3, BLCA 4/3, STAD 4).

**Outputs:**
- **Restricted:** `results/Pancancer_epithelial_annotation.xlsx` (sheets
  CIOC_8, ECOS_38, Extended_222, Candidates_338, Calibrators, Cancer_types,
  Legend) and `results/Pancancer_epithelial_CIOC8.pdf`.
- **Script:** `scripts/pancancer_02_reactivation.py`.
- **Relationship to the bulk TCGA/GTEx annotation (breadth plan v2):** the
  bulk annotation is retained as a secondary annotation only.

## Amendment A1: expression scale (frozen 2026-10-05, after the v1 run, before recomputation)

**Problem.** The frozen transform y = log2(mean CP10k + 1) uses a pseudocount
of 1 CP10k, which equals **100 CPM**. That is about 100× the CPM + 1 scale
of the edgeR CRC gate, yet the plan states that its thresholds "match the
CIOC CRC gate". Genes below about 100 CPM are compressed towards zero, and
this covers most CIOC genes in epithelium.

**How it was found.** The v1 run returned no CIOC gene reactivated in the
independent CRC dataset. Inspection of the CRC arm means showed, for
example, TACSTD2 at 0.05 vs 0.80 CP10k (5 vs 80 CPM, about 16-fold) in
normal vs malignant epithelium, which v1 scored as about +0.8 on its
compressed scale and as non-significant.

This is a units error, not a threshold choice. It is disclosed as an
amendment made after seeing the v1 output.

**Change.** y = log2(CPM + 1), with CPM = mean CP10k × 100. Nothing else
changes: cells, units, tiers, models, thresholds (log2FC ≥ 0.5, BH
FDR < 0.05), breadth definitions and exclusions are as frozen.

**Record.**
- The v1 outputs are archived unmodified in DATA restricted
  `pancancer_v1_scale_error/`, with md5 values in WORKLOG.
- The v1 summary is in WORKLOG: no CIOC gene reactivated in independent
  CRC, and at most 2/12 other carcinomas for any CIOC gene.
- The amended analysis is applied once.
