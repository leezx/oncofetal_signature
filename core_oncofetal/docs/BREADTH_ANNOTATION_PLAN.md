# Tissue and cancer breadth annotation: plan (frozen 2026-10-05, before computation)

## Purpose and scope

This is a supplementary **annotation** of each marker, not signature
construction. It answers one question:

> Is a gene a broad pan-tissue developmental / pan-cancer marker, or is it
> intestine- and CRC-biased?

It adds no gate. It does not change the CIOC (8), the 338 candidates, the
Extended CIOC (222) or ECOS-38. No membership changes are made after viewing
the results.

- **Genes annotated:** the 338 cross-species fetal–CRC candidates, which
  contain the CIOC, the Extended CIOC and ECOS-38, plus all Literature-31
  genes as context.
- **Membership labels** are restricted (Joanito-derived), so per-gene
  outputs are restricted.

One dataset per axis, by review decision (no further datasets):

| Axis | Dataset | Question |
|---|---|---|
| Fetal tissue breadth | Cao et al., Science 2020, fetal atlas (GSE156793) | Is fetal epithelial expression intestine-biased or shared across organs? |
| Cancer reactivation breadth | UCSC Toil TCGA + GTEx (uniform RSEM TPM) | Is tumour-vs-normal gain CRC-biased or shared across carcinomas? |

## A. Fetal tissue breadth (Cao 2020)

**Data.** The GEO aggregated tables are used:
- S6: pseudobulk CPM per organ × cell type (each column sums to 10⁶);
- S7: the fraction of cells detecting the gene;
- S1: cell counts, fetus counts, ages and assay type.

The ages are 72–129 days post-conception (about 10–18 weeks), which matches
the ≥ 9-week developmental window. Ensembl IDs are mapped to symbols through
S2, and current HGNC symbols are resolved through previous symbols.

**Organ epithelia (prespecified).** These are the author cell types that
are epithelial. An organ is used if it has ≥ 200 epithelial cells; eye
corneal/conjunctival epithelium (104 cells) is excluded.

| Organ | Epithelial cell types | Assay |
|---|---|---|
| Intestine | Intestinal epithelial cells | whole cell |
| Stomach | Ciliated epithelial; Goblet; MUC13_DMBT1-positive; Neuroendocrine; Parietal and chief; Squamous epithelial | whole cell |
| Pancreas | Acinar; Ductal; Islet endocrine | whole cell |
| Kidney | Metanephric; Ureteric bud | whole cell |
| Lung | Bronchiolar and alveolar epithelial; Ciliated epithelial; Neuroendocrine; Squamous epithelial | nuclei |
| Liver | Hepatoblasts | nuclei |
| Placenta | Extravillous trophoblasts; Syncytiotrophoblasts and villous cytotrophoblasts; Trophoblast giant cells | nuclei |
| Thymus | Thymic epithelial cells | nuclei |

**Organ-epithelium expression.** This is the cell-count-weighted mean of
the S6 CPM over the organ's epithelial cell types.

- **Expressed:** CPM ≥ 10 in that organ epithelium. The threshold is fixed
  a priori. Detection fractions are not used because the sci-RNA-seq3 data
  are shallow (EPCAM is detected in 28% of intestinal epithelial cells).
- **n_expr:** the number of the 8 organ epithelia in which the gene is
  expressed.

**Fetal breadth class:**

| Class | Definition |
|---|---|
| not detected in fetal intestinal epithelium | intestinal epithelium CPM < 10 |
| pan-fetal epithelial | intestine expressed and n_expr ≥ 5 of 8 |
| multi-organ fetal epithelial | intestine expressed and n_expr 3–4 |
| intestine-biased fetal | intestine expressed and n_expr ≤ 2 (intestine plus at most one other organ) |

**Annotations:**
- **Fetal-intestine compartment flag:** log2((intestinal epithelial CPM + 1)
  / (max intestinal non-epithelial cell-type CPM + 1)) < −3, the same logic
  as Gate E. This marks genes whose fetal intestinal signal is
  predominantly non-epithelial (for example myeloid).
- **Cells-only breadth:** n_expr over the 4 whole-cell organs (intestine,
  stomach, pancreas, kidney).
- **τ specificity index** across the 8 organ epithelia, computed on
  log2(CPM + 1).

**Caveats (disclosed):**
- The tables aggregate across fetuses, so there is no donor-level
  replication.
- The assay differs by organ (whole cells vs nuclei), so cross-organ
  comparisons are partly confounded; this is why the cells-only breadth is
  reported.
- Small stomach, pancreas and thymus epithelial groups come from 3–4
  fetuses.

## B. Cancer reactivation breadth (Toil TCGA + GTEx)

**Data.** `TcgaTargetGtex_rsem_gene_tpm`, stored as log2(TPM + 0.001) and
converted to log2(TPM + 1). The phenotype file is `TcgaTargetGTEX_phenotype`.

**Cancer types: carcinomas only (prespecified).** Each TCGA primary tumour
is compared with a normal reference made of GTEx matched tissue plus TCGA
solid-tissue normals of the same type:

| TCGA | GTEx reference |
|---|---|
| COAD, READ | Colon - Transverse |
| STAD | Stomach |
| ESCA | Esophagus - Mucosa |
| PAAD | Pancreas |
| LIHC, CHOL | Liver |
| LUAD, LUSC | Lung |
| BRCA | Breast - Mammary Tissue |
| KIRC, KIRP, KICH | Kidney - Cortex |
| BLCA | Bladder |
| PRAD | Prostate |
| THCA | Thyroid |
| UCEC | Uterus |
| CESC | Cervix (Ecto- + Endocervix) |
| OV | Ovary |
| ACC | Adrenal Gland |
| HNSC | none (TCGA normals only) |

- **Eligibility:** ≥ 10 tumours and ≥ 10 normal-reference samples.
- **Excluded:** non-carcinoma types (glioma, melanoma, sarcoma,
  haematological, germ-cell, thymoma, mesothelioma, paraganglioma, uterine
  carcinosarcoma).
- **Known reference mismatches** (reported, not corrected): OV against
  ovary (stromal tissue), and CHOL against liver.

**Statistic, per gene and cancer type:**
- Δ = median log2(TPM + 1) in tumours − median in the normal reference.
- Two-sided Mann–Whitney test.
- BH correction within each cancer type, across the annotated gene universe
  (338 ∪ Literature-31 ∪ calibrators).
- **Up:** Δ ≥ 1 and FDR < 0.05. **Down:** Δ ≤ −1 and FDR < 0.05.

**Cancer breadth class.** "CRC-up" means up in COAD or READ.

| Class | Definition |
|---|---|
| pan-cancer reactivated | up in ≥ 50% of eligible non-CRC carcinomas |
| multi-cancer reactivated | up in 25% to < 50% of eligible non-CRC carcinomas |
| CRC-biased reactivated | CRC-up and up in < 25% of eligible non-CRC carcinomas |
| not CRC-up in bulk | not up in COAD or READ, and up in < 25% of non-CRC carcinomas |

**Caveats (disclosed):**
- Bulk data reflect cell composition, as Gate E showed; tumour-vs-normal
  differences include stromal and immune shifts.
- GTEx and TCGA normals differ in batch.
- The reference tissue is not always the cell of origin.

## C. Combined classification (two dimensions)

- **Fetal axis:** pan-fetal / multi-organ / intestine-biased / not detected
  in fetal intestinal epithelium.
- **Cancer axis:** pan-cancer / multi-cancer / CRC-biased / not CRC-up.
- **Wording:**
  - "pan-fetal epithelial" and "intestine-biased fetal";
  - "pan-cancer reactivated" and "CRC-biased reactivated";
  - "pan-tissue pan-cancer oncofetal marker" only when a gene is **both**
    pan-fetal epithelial and pan-cancer reactivated.
  - Neither axis alone establishes "pan-oncofetal". Cao answers fetal
    tissue breadth, and TCGA/GTEx answers cancer versus adult normal.

## D. Calibrators (reported, never used to tune)

| Group | Genes |
|---|---|
| Expected pan-epithelial | EPCAM, KRT8, CDH1 |
| Expected intestine-biased | CDX2, VIL1, CDH17 |
| Non-epithelial | PTPRC, COL1A2 |
| Expected pan-cancer up | MKI67, TOP2A |

## E. Outputs

- **Restricted** (membership labels; git-ignored; mirrored to DATA):
  `results/Breadth_annotation.xlsx` with these sheets:
  - CIOC_8, ECOS_38, Extended_222, Candidates_338, Literature_31;
  - Cancer_types, Calibrators, Legend;
  - `results/Breadth_CIOC8_heatmap.pdf` (the main-text candidate figure),
    plus the 2D classification plot.
- **Public:** the code and this plan.
- **Script:** `core_oncofetal/scripts/breadth_annotation.py`.

## Application (2026-10-05; applied once)

Produced by `scripts/breadth_annotation.py`. Per-gene results are in the
restricted `results/Breadth_annotation.xlsx`; the figures are
`Breadth_CIOC8_heatmap.pdf` and `Breadth_2D_classification.pdf`, also
restricted.
- **Cancer types:** all 21 prespecified carcinomas were eligible (19 non-CRC).
- **Fetal organs:** 8 organ epithelia.

**Fetal tissue breadth:**

| Set | Pan-fetal epithelial | Multi-organ | Intestine-biased | Not detected in fetal intestinal epithelium | Not measured |
|---|---|---|---|---|---|
| CIOC (8) | 3 | 1 | 1 | 3 | 0 |
| ECOS-38 | 26 | 5 | 1 | 6 | 0 |
| Extended CIOC (222) | 159 | 10 | 7 | 43 | 3 |
| Cross-species (338) | 219 | 19 | 11 | 84 | 5 |

The fetal-intestine compartment flag is set for 96 of the 338 genes.

**Cancer reactivation breadth:**

| Set | Pan-cancer | Multi-cancer | CRC-biased | Not CRC-up in bulk | Not measured |
|---|---|---|---|---|---|
| CIOC (8) | 2 | 1 | 1 | 4 | 0 |
| ECOS-38 | 3 | 11 | 2 | 22 | 0 |
| Extended CIOC (222) | 12 | 59 | 33 | 115 | 3 |
| Cross-species (338) | 14 | 82 | 39 | 198 | 5 |

**Combined.** "Pan-tissue pan-cancer oncofetal" genes:

| Set | Genes |
|---|---|
| CIOC | 0 |
| ECOS-38 | 2 |
| Extended CIOC | 10 |
| Cross-species 338 | 12 |

**Calibrators (reported, not used for tuning):**
- **Behaved as expected:**
  - EPCAM, KRT8 and CDH1 are pan-fetal epithelial and pan-cancer.
  - CDX2 and CDH17 are intestine-biased fetal.
  - MKI67 and TOP2A are pan-cancer.
  - PTPRC is not detected in fetal intestinal epithelium.
- **COL1A2** (stromal) classifies as "pan-fetal epithelial", because
  stroma-level expression leaks into Cao epithelial pseudobulks at ≥ 10 CPM.
  The compartment flag catches it (−6.9).
  - **The fetal class must therefore be read together with the compartment
    flag.**
- **MKI67** is also flagged (erythroblasts).
- **Cao intestinal "Chromaffin cells"** express CDX2 and EPCAM (probably
  enteroendocrine cells) but are counted as non-epithelial by the frozen
  plan. This makes the flag slightly conservative against epithelial genes.

**Interpretation limits:**
- **The cancer axis measures tumour vs the patient's own tissue type.**
  Genes that define normal intestinal identity (for example CDX2) are not
  "CRC-up".
- **Bulk composition affects the cancer axis.** A gene's pan-cancer gain can
  reflect immune or stromal infiltration rather than tumour-cell expression.

## Revision v2: cancer breadth only (frozen 2026-10-05, before recomputation)

This revision was made after the v1 results were seen, by review decision.
It changes scope and interpretation labels only. The cancer-axis data,
cancer types, statistic and thresholds (Δ ≥ 1, FDR < 0.05; 25% and 50%
cut-offs) are unchanged.

**Rationale.** The CIOC defines oncofetal **identity**: fetal intestinal
epithelium > adult intestinal epithelium, conserved in mouse, and CRC
epithelium > normal epithelium. Fetal *tissue* specificity (fetal intestine
vs other fetal organs) and *cancer* specificity are different questions,
and neither is part of the definition.
- A gene expressed in several fetal organs, or reactivated in several
  carcinomas, is still a valid intestinal oncofetal marker.
- Organ-of-origin labels for shared developmental programs lead to an
  ontology problem that does not need solving.

**Scope changes:**
1. **The fetal tissue breadth axis (Cao 2020) is removed from all
   deliverables.**
   - The v1 code (`core_oncofetal/exploratory/breadth_fetal_cancer_v1.py`)
     and its outputs are kept as internal exploration only. The outputs are
     in DATA restricted `exploratory_breadth_v1/`.
   - It is not presented in the manuscript or the supplement.
   - The terms "pan-fetal", "pan-tissue" and "pan-oncofetal" are not used.
2. **Title:** "Pan-cancer breadth of CIOC reactivation in bulk tumours".
   It is an annotation, not a validation, and creates no hierarchy over the
   epithelial-resolved CIOC evidence.
3. **Bulk cannot overrule epithelial evidence.** Bulk expression is the sum
   over cell types of cell fraction × expression. A gene without a bulk CRC
   gain gets no cancer-breadth class.

**Per-gene annotation** (genes: the 338 cross-species candidates, with CIOC,
Extended CIOC and ECOS-38 labels):

| Column | Definition |
|---|---|
| CRC epithelial evidence | Gate C (Joanito + Pelka epithelial pseudobulk). Positive for all 338 by construction. |
| CRC bulk | **up** = Δ ≥ 1 and FDR < 0.05 in COAD or READ; **down** = Δ ≤ −1 and FDR < 0.05 in COAD or READ (and up in neither); else **n.s.** |
| Other carcinomas up | n of 19 eligible non-CRC carcinomas |
| Composition-sensitive | Gate E E1 FAIL in the Khaliq/Che atlases (strong non-epithelial attribution); this is an existing annotation |

**Interpretation:**

| CRC bulk | Other carcinomas up | Interpretation |
|---|---|---|
| up | < 25% | **CRC-biased reactivation** |
| up | 25% to < 50% | **multi-cancer reactivation** |
| up | ≥ 50% | **broad reactivation across carcinomas** |
| down | any | **epithelial reactivation not captured by bulk (CRC bulk down)**; bulk cancer breadth not interpretable |
| n.s. | any | **bulk cancer breadth not interpretable (CRC bulk n.s.)** |

- The suffix ", composition-sensitive" is added when E1 fails.
- The count of other carcinomas is always reported, but it is descriptive
  when CRC bulk is not up.

**Outputs:**
- **Restricted:** `results/Cancer_breadth_annotation.xlsx` (sheets CIOC_8,
  ECOS_38, Extended_222, Candidates_338, Cancer_types, Legend) and
  `results/Cancer_breadth_CIOC8.pdf` (CIOC heatmap, carcinoma panel only).
- **Script:** `core_oncofetal/scripts/cancer_breadth_annotation.py`.
- The v1 outputs (`Breadth_annotation.xlsx`, `Breadth_CIOC8_heatmap.pdf`,
  `Breadth_2D_classification.pdf`) are removed from `results/`.

### Application of v2 (2026-10-05; applied once)

Produced by `scripts/cancer_breadth_annotation.py`. The cancer statistic is
imported unchanged from the v1 code. Per-gene results are in the restricted
workbook `results/Cancer_breadth_annotation.xlsx`; the figure is
`Cancer_breadth_CIOC8.pdf`.

**Interpretation counts** (n; of which composition-sensitive):

| Set | CRC-biased | Multi-cancer | Broad | Epithelial reactivation not captured by bulk (CRC bulk down) | Bulk breadth not interpretable (CRC bulk n.s.) | Not measured |
|---|---|---|---|---|---|---|
| CIOC (8) | 1; 0 | 1; 0 | 2; 1 | 1; 1 | 3; 2 | 0 |
| ECOS-38 | 2; 0 | 8; 0 | 3; 0 | 4; 0 | 21; 0 | 0 |
| Extended CIOC (222) | 34; 0 | 37; 2 | 11; 2 | 17; 10 | 120; 19 | 3 |
| Cross-species (338) | 40; 0 | 43; 4 | 12; 3 | 39; 19 | 199; 43 | 5 |

**Reading.** Most candidates have no interpretable bulk breadth, because
CRC bulk does not show their epithelial gain. This is consistent with
cell-intrinsic epithelial reactivation being diluted or opposed by
composition in bulk tissue. It does not count against the
epithelial-resolved evidence.
