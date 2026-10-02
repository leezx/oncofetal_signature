# Conserved Intestinal Oncofetal Core — Methods (validated gate design)

- Version: Core method v1.0, 2026-10-02.
- Status: **frozen** before any CRC-gate value of the developmental-gate
  genes was read. The two developmental gates and their result were reviewed
  first (benchmark v2 addendum v1.4).
- Inputs are existing, version-controlled contrasts. **No statistic is
  recomputed for gene selection.**

## 1. Principle

The Core is the subset of literature-nominated intestinal oncofetal genes
that passes four sequential evidence gates:

- a human in-vivo developmental gate;
- a mouse in-vivo developmental gate;
- a CRC malignant-epithelium gate;
- an independent CRC replication gate.

Every gate is a **same-study, epithelium-resolved, replicate-level**
contrast. Supportive datasets are reported beside the Core but **never take
part in selection**. Thresholds are fixed and are not relaxed to enlarge the
Core: a Core of three to five genes is acceptable, because the Core is not
a signature.

## 2. Candidate universe

- **Source:** the 31 Step 1 literature-curated intestinal oncofetal markers
  (`step2_fetal/config/literature_candidates_31.tsv`).
- **Symbols:** current HGNC symbols, with legacy source symbols resolved
  (CCN1 = CYR61, CCN2 = CTGF; `step3_cancer/config/symbol_aliases.tsv`).
- **Mouse-defined markers:** LY6A and REG3B have no human one-to-one
  orthologue, so they cannot pass the human gates.
- **SPRR1A** has no mouse orthologue (Ensembl 116), so it cannot pass the
  mouse gate.
- **The pool is deliberately heterogeneous.** It mixes direct fetal genes,
  regeneration/revival genes, YAP-responsive genes, injury genes and
  oncofetal cancer genes. The gates supply the selection pressure.

## 3. Gate design rules (applied before gate choice)

- **Design eligibility comes first, literature recovery second.** A gate
  dataset must meet all of E1–E5 (benchmark v2 addendum v1.4):
  - E1: fetal and adult (or tumour and normal) come from the same study;
  - E2: primary in-vivo tissue;
  - E3: epithelium-resolved;
  - E4: ≥ 3 biological replicates per arm;
  - E5: a replicate-level model with per-gene FDR.
- **Among eligible datasets, literature recovery is the comparison metric.**
  It is measured as % of the 31 markers that are fetal-positive, or
  fetal-positive with FDR < 0.05 and log2FC ≥ 0.5.
- **TNFRSF12A is not an admission criterion** for any gate.
- **Developmental stage.** Human fetal donors < 9 weeks are excluded. The
  9-week breakpoint was derived from fetal-only data (Gao and Fawkner manual
  state analysis) before any adult comparison was re-opened. It was
  independently supported by HGCA: mid/late vs early fetal TNFRSF12A
  +1.60, P = 0.009.

## 4. Gates

A gene passes a gate only if it is **measured** (it passes the contrast's
expression filter), **log2FC ≥ 0.5** in the stated direction, and **FDR < 0.05**
(Benjamini–Hochberg within the contrast). An unmeasured gene fails.

| Gate | Dataset | Contrast | Unit, model | Direction |
|---|---|---|---|---|
| G1 Human developmental | Human Gut Cell Atlas (Elmentaite 2021) | fetal epithelium ≥ 9 PCW (10 donors, 9.2–17 PCW) vs healthy adult epithelium (7 donors) | donor raw-count pseudobulk of author-annotated epithelial cells (doublets and lymph node excluded, ≥ 50 cells); edgeR QL `~ stage`, `filterByExpr`, TMM | fetal > adult |
| G2 Mouse developmental | Pikkupeura 2023, GSE230581 (in-vivo arm; paper Fig 1D–G) | freshly isolated E16.5 proximal SI epithelium (3) vs adult proximal SI crypt epithelium (3) | biological replicate; edgeR QL; Ensembl 116 one-to-one orthologue (CXADR: high-confidence one2many) | fetal > adult |
| G3 CRC malignant | Joanito 2022 (Synapse syn26844071) | malignant epithelium vs normal epithelium (S1) | patient pseudobulk; edgeR QL `~ cohort + group`; frozen label map | malignant > normal |
| G4 CRC replication | Pelka 2021 (GSE178341) | tumour epithelium vs normal epithelium (S2) | patient pseudobulk; edgeR QL | tumour > normal |

**Core = genes that pass G1 ∧ G2 ∧ G3 ∧ G4.**

The **Developmental Core** (G1 ∧ G2) is reported as an intermediate set.
Under the frozen thresholds it is, from benchmark v2 addendum v1.4:
**GJA1, CLU, ANXA6, SPP1, RBP1.**

## 5. Supportive evidence (reported, never used for selection)

| Evidence | Role | Contrast |
|---|---|---|
| H-new2: Fawkner 2021 fetal ≥ 9 PCW vs Burclaff 2022 adult | human independent support (cross-study) | donor/sample pseudobulk, edgeR |
| Gao 2018 fetal LI ≥ 9 W vs GSE103154 adult LI | human support | effect-only (cross-platform, adult n = 2) |
| HGCA all fetal vs adult | human, stage-unresolved reference | edgeR |
| GSE44433: Hemmerling 2014, WT E17.5 vs 8-week LCM ileal epithelium | mouse independent in-vivo support | microarray, limma |
| Pikkupeura 2021 cultures (GSE160449), laminin and collagen | developmental/culture support | fetal enterospheres vs adult organoids |
| Joanito sensitivity: cohorts with both groups | CRC robustness | edgeR |
| Pelka progenitor check P: tumour vs normal stem/TA | CRC biology (not a gate) | edgeR |
| TCGA COAD/READ: tumour vs normal (T1), paired, vs GTEx | bulk CRC support | edgeR |

A supportive dataset that disagrees with a Core gene is reported and
discussed; it does not remove the gene.

## 6. Known design caveats (disclosed)

- **G1:** HGCA fetal (HDBR) and adult (transplant-donor) cohorts were
  collected and processed within one atlas framework but as separate
  cohorts. Residual cohort effects cannot be excluded.
  - Fetal tissue is small and large intestine; adult tissue spans the
    duodenum-to-rectum regions sampled.
- **G2:** the fetal arm is whole proximal SI epithelium and the adult arm is
  crypt-only epithelium, so the compartments are asymmetric.
- **G3/G4:** tumour epithelium is malignant cells and normal is
  normal-adjacent or healthy epithelium; patient pseudobulk.
  - Joanito is the primary gate; Pelka is the replication gate.
  - TCGA is bulk and does not veto.
- **Human-panel recovery in G1 is moderate:** 15/25 measured
  fetal-positive, 7/25 strict. Several literature candidates that are
  injury/revival genes are adult-high in G1 (ANXA1, AREG, CD44, EMP1, IL33,
  IL1RN, EPS8L1). This is treated as selection pressure, not dataset
  failure.

## 7. Expected and reported non-members

TNFRSF12A is literature-nominated but fails both developmental gates:

| Gate | log2FC | P / FDR |
|---|---|---|
| G1 | −0.01 | P = 0.99 |
| G2 | +0.46 | FDR not significant |

It is therefore not a Core member. Statement for the manuscript:
*TNFRSF12A is mechanistically linked to the oncofetal program but is not
itself a member of the stringent conserved oncofetal Core.* It was not
removed manually.

## 8. Outputs and data-use restriction

- **Gate statistics:** `core_oncofetal/results/` holds the evidence workbook,
  with one column per gate (log2FC, P, FDR, pass).
- **Joanito restriction:** the Joanito data-use terms forbid disclosure of
  derived material. The **full** workbook (G3 values and final Core
  membership) is written locally and mirrored to
  `DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/core_oncofetal/`.
  It is git-ignored.
- **Public workbook:** the public version omits the G3 column and the final
  Core label. It keeps the candidate universe, G1, G2, G4 and the supportive
  evidence.

## 9. Change control

Any change to a gate dataset, contrast, threshold or direction is a new
method version, recorded here before it is applied.
