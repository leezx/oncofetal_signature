# Pan-cancer epithelial reactivation of CIOC genes: plan (DRAFT; freeze pending the metadata scan)

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

## Analysis (to be frozen)

- **Tumour arm:** CNV-inferred malignant epithelial (or parenchymal) cells
  (`cnv_status` = tumour) from primary-tumour samples. These are "malignant"
  by the authors' CNV inference; "tumour-derived" is used where no CNV call
  exists.
- **Normal arm:** epithelial (or parenchymal) cells with `cnv_status` = normal,
  from adjacent-normal or healthy-normal samples of the matching organ.
- **Unit:** patient × arm pseudobulk, pooling that patient's samples. A unit
  needs ≥ 20 cells. Values are log2(mean CP10k + 1). **Never cell-level
  tests.**
- **Evidence tier per cancer type:**
  - **Tier 1 (within-study):** datasets with both arms (≥ 2 units each).
    The model is y ~ arm + dataset (OLS), with ≥ 3 units per arm in total.
  - **Tier 2 (cross-study; weaker):** if Tier 1 is unavailable, tumour
    units vs normal units of the same organ from any dataset, using y ~ arm
    and ≥ 3 units per arm.
  - Otherwise the cancer type is not evaluable.
- **Reactivation:** arm coefficient ≥ 0.5 and BH FDR < 0.05 within the cancer
  type, across the gene universe. This matches the CIOC CRC gate threshold.
- **Breadth:** k / n = the number of evaluable non-CRC carcinomas with
  reactivation over the number evaluable, reported as a count. No pan-cancer
  dichotomy is imposed. The CRC column is reported separately, and a
  Tier 1-only breadth is reported as a sensitivity analysis.
- **Carcinomas only:** non-epithelial malignancies are excluded (glioma,
  melanoma, uveal melanoma, lymphoma, myeloma, leukaemia, sarcoma,
  neuroblastoma, neuroendocrine tumours).

## Pending before freeze (metadata only)

1. Tissue values that define primary tumour, adjacent normal and healthy
   normal.
2. Celltype labels counted as epithelial/parenchymal lineage (for example
   hepatocytes for HCC and renal epithelium for RCC).
3. The Cancer type → normal organ mapping and the carcinoma list.
4. Units available per cancer type and tier (counts only).
