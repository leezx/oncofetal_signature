# Step 2 human benchmark v2 — report

Date: 2026-10-02. Rules: [`QUALIFICATION_PLAN.md`](QUALIFICATION_PLAN.md)
(frozen in `5c3bcef` before any new dataset was opened). No signature was built.

## Verdicts (frozen rule: TNFRSF12A fetal-high with P < 0.05 **and** panel binomial P < 0.05)

| Contrast | TNFRSF12A log2FC (P) | Markers fetal-positive | Binomial P | Verdict |
|---|---|---|---|---|
| HGCA (original H1) | −0.46 (0.38) | 14/25 | 0.35 | fails |
| Gao-original, recomputed (fetal LI vs GSE103154) | −0.08 (0.45) | 14/28 | 0.57 | fails |
| Senger enterospheres | −0.22 (0.65) | 14/29 | 0.64 | fails |
| **H-new1** Fawkner Visium, fetal vs adult | −1.71 (0.001) | 14/24 | 0.27 | fails |
| H-new1 colon only | −3.01 (<0.001) | 13/24 | 0.42 | fails |
| **H-new2** Fawkner fetal scRNA vs Burclaff adult | **+1.83 (<0.001)** | **18/24** | **0.011** | **qualifies** |
| **H-new3** Gao fetal SI+LI vs Wang adult | +0.31 (0.81) | 17/24 | 0.032 | fails (rule 1) |
| H-new3-debug Gao fetal LI vs Wang colon+rectum | +0.70 (0.41) | 17/24 | 0.032 | fails (rule 1) |
| Pikkupeura LN (mouse culture) | +1.79 (FDR < 0.001) | 24/27 | < 0.001 | qualifies |
| Pikkupeura collagen (mouse culture) | +2.13 (FDR < 0.001) | 26/27 | < 0.001 | qualifies |

Priority controls in H-new2: TACSTD2 +3.97, CLU +2.25, TNFRSF12A +1.83 (all
P < 0.001); ANXA1 −0.95 (P = 0.16). Full values:
`results/tables/Benchmark_v2_31_marker_matrix.csv`; figure:
`results/figures/Benchmark_v2_31_marker_heatmap.{pdf,png}`.

## What each new contrast shows

**H-new1 (Visium) is dominated by tissue composition.** Fetal "epithelial-rich"
spots carry far less epithelial signal than adult spots (epithelial-marker UMI
fraction 0.4–0.8% vs 2.1%; stromal/immune fraction 0.4–0.8% vs 0.02%), and
EPCAM itself is adult-high (log2FC −2.6) while COL1A1/VIM are fetal-high. At
55 µm resolution, thin 12–19 PCW epithelium is diluted by mesenchyme, so this
contrast measures epithelial content more than developmental state. The adult
arm is also only two colon sections (possibly one donor). Not usable as
primary developmental evidence in its current form.

**H-new2 is the first human contrast that qualifies.** Hashtag-demultiplexed
fetal EPCAM+ samples (22 eligible, 8–22 PCW) vs Burclaff healthy adult
epithelium (3 donors, all six regions) recover TNFRSF12A, TACSTD2 and CLU as
strongly fetal-high. Caveats: study/platform is completely confounded with
stage; some hashtag samples are different regions of the same fetus (pseudo-
replication on the fetal side); adult donors span duodenum to distal colon
whereas fetal samples are terminal ileum and colon.

**H-new3 and the Gao debug pair localise the Gao problem to the fetal arm.**
Both adult references give similar TNFRSF12A levels (GSE103154 P1/P2: 26.2 /
18.5; Wang colon+rectum: 9.1–18.7 per million). Swapping the adult arm moves
the effect from −0.08 to +0.70 without statistical support, so by the frozen
reading table the problem lies in the **Gao fetal compartment / stage /
composition**, not the adult reference. Descriptively, Gao fetal TNFRSF12A
varies 34-fold across embryos (2.6–89.8), with the lowest values in 6–9-week
embryos. Panel bias (17/24) is significant in both H-new3 contrasts.

## Gene coverage notes

- Burclaff's published h5ad keeps 23,170 genes (SPP1, ANKRD1, SOX17 removed by
  the authors' gene filter); counts for the same annotated cells are therefore
  taken from the per-donor full-gene matrices (32,732 genes shared with
  Fawkner).
- Fawkner GEO titles mis-pair expression and hashtag libraries; pairs were
  fixed by ≥ 99% cell-barcode overlap (EPI2–HTO3, EPI3–HTO5, EPI pool 4 –
  "HTO_stromal_4").
- TACSTD2 is filtered as low-expression in Visium; LY6A/REG3B have no human
  one-to-one orthologue.

## Implication (for review; nothing is frozen from this)

- The literature-marker failure of Step 2 is **dataset-specific, not
  universal**: a cleanly epithelial fetal scRNA vs adult epithelial scRNA
  contrast recovers the expected fetal direction for TNFRSF12A, TACSTD2 and
  CLU.
- HGCA, Gao and Visium fail for identifiable reasons (adult reference /
  composition in HGCA is still unresolved; Gao fetal heterogeneity and stage;
  Visium spot composition).
- A replacement human primary developmental contrast should be built around
  H-new2-like data, with a second independent fetal epithelial scRNA source to
  break the study-stage confound. TNFRSF12A by gestational age × region ×
  epithelial subtype (Gao, Fawkner) is the natural manual follow-up.

## Review decisions (2026-10-02)

1. **H-new2 becomes the new human primary developmental contrast**, labelled a
   **cross-study developmental contrast** (Fawkner fetal vs Burclaff adult;
   study/platform confounded with stage). It requires independent validation
   and is not to be treated as final proof on its own.
2. **Visium (H-new1, H-new1-colon) is dropped from quantitative Step 2.** Its
   outputs are retained only as a failed dataset/QC record; no further work is
   planned on it.

Decisions are recorded per contrast in `config/review_decisions.tsv` and appear
as `review_decision` in `results/tables/Benchmark_v2_qualification.csv`; they
do not alter the frozen verdicts. Not yet decided: the manual TNFRSF12A
follow-up (gestational age × region × epithelial subtype).

## Addendum v1.1 — primary-tissue bulk (H-bulk1, H-bulk2)

Frozen in `cd39210` before data were opened. Source tracing showed that the
Finkbeiner 2015 fetal/adult comparison and the Senger 2018 "primary tissue"
comparison use the **same public samples**: 6 Roadmap Epigenomics fetal
small-intestine total-RNA libraries (SRP001371; six donors, 91–120 days) and
HPA E-MTAB-1733 adult small intestine/duodenum (ERP003613). H-bulk2 (duodenum
only) is therefore a subset of H-bulk1, not an independent dataset. Both arms
were taken from recount3 as raw counts (the authors share FPKM only) and
analysed with edgeR `~ stage`.

| Contrast | TNFRSF12A log2FC (P) | Fetal-positive | Binomial P | OLFM4 log2FC (sanity) | Verdict |
|---|---|---|---|---|---|
| H-bulk1 fetal SI (6) vs adult SI + duodenum (6) | −0.50 (0.20) | 13/26 | 0.58 | −6.87 | fails |
| H-bulk2 fetal SI (6) vs adult duodenum (2) | −0.90 (0.12) | 16/27 | 0.22 | −6.83 | fails |

Priority markers in H-bulk1: TACSTD2 +2.49 (P = 1e-4), CLU +0.14 (n.s.),
ANXA1 −2.66 (P < 1e-5), TNFRSF12A −0.50 (n.s.); GJA1 +0.55 (n.s.), RBP1 +1.50
(P < 1e-4), BASP1 +1.10 (P = 0.03).

**Reading.** Whole-tissue developmental contrasts recover established
intestinal maturation signals (OLFM4 −6.87 / −6.83; TACSTD2, RBP1, BASP1
fetal-high) but do not reproduce the epithelial TNFRSF12A fetal-high pattern,
potentially reflecting tissue-composition and library-preparation differences
(whole tissue in both arms; fetal total RNA vs adult poly(A)). This is an
interpretation, not a demonstrated cause: composition, platform, region,
developmental stage, and TNFRSF12A being a specific epithelial-state marker
cannot yet be distinguished.

**Evidence counting.** H-bulk1 and H-bulk2 share the same six fetal samples and
the same HPA source; they count as **1 independent dataset, 2 related
contrasts**, never as two independent bulk validations
(`config/review_decisions.tsv`, column `independent_source`).

**Evidence tiers** (`config/review_decisions.tsv`, column `evidence_tier`):
Tier A epithelial scRNA — H-new2 (qualifies); Tier B primary-tissue bulk —
H-bulk1, H-bulk2 (both fail); Tier C problematic references — HGCA, Gao,
Senger enterospheres (fail); Visium dropped.

This matches the stop condition stated in the review: with the bulk contrasts
also showing TNFRSF12A fetal ≤ adult, the automated benchmark stops.

## Review decisions — round 3 (2026-10-02)

1. The PR is kept; **the bulk benchmark is closed and no further human
   fetal/adult datasets are sought.**
2. H-bulk1/H-bulk2 = 1 independent dataset, 2 related contrasts.
3. The bulk failure is reported with the wording above (recovered maturation
   signal, TNFRSF12A not reproduced, cause not established).
4. The human evidence is read as **compartment-dependent**:

   | Compartment | TNFRSF12A fetal/oncofetal-high? |
   |---|---|
   | Whole tissue (H-bulk1/2, Visium) | not supported |
   | Mixed/heterogeneous epithelial references (HGCA, Gao) | unstable |
   | Clean epithelial scRNA (H-new2) | fetal-high |
   | Fetal epithelial culture (Pikkupeura, mouse) | fetal-high |
   | CRC malignant epithelium (Step 3) | strongly up |

5. **H-new2 is not used on its own for genome-wide signature discovery** yet
   (dataset-selection concern: it is the only human contrast passing the
   frozen TNFRSF12A rule). The final human discovery contrast is decided after
   the manual developmental-state analysis of Gao and Fawkner
   ([`step2_fetal_state/`](../../step2_fetal_state/)), which asks *when, where
   and in which epithelial state* TNFRSF12A is high in human fetal intestine.

## Addendum v1.2 — stage-resolved re-analysis (Gao, HGCA; Senger not splittable)

The plan was frozen in `618356e` before any stage-split contrast was
computed.

- **Breakpoint:** 9 weeks, taken from the fetal-only analysis in
  `step2_fetal_state/`.
- **Groups:** early fetal < 9 W and mid/late fetal ≥ 9 W, applied to the
  age as each source reports it.
- **Methods:** units, filters and statistics are unchanged, and the frozen
  rule is applied to every contrast. The recomputed HGCA all-fetal contrast
  reproduces the original (TNFRSF12A −0.46, P = 0.38).

| Dataset | Contrast | TNFRSF12A log2FC (P) | Panel fetal-positive | Binomial P | Verdict |
|---|---|---|---|---|---|
| HGCA | all fetal vs adult (16 vs 7) | −0.46 (0.38) | 14/25 (56%) | 0.35 | fails |
| | early (6) vs adult | −1.58 (0.006) | 13/23 (57%) | 0.34 | fails |
| | mid/late (10) vs adult | −0.01 (0.99) | 15/25 (60%) | 0.21 | fails |
| | mid/late vs early | **+1.60 (0.009)** | 14/24 | 0.27 | (within-fetal) |
| Gao-original (LI vs GSE103154) | all fetal (12) vs adult | −0.08 (0.45) | 14/28 (50%) | 0.57 | fails |
| | early (5) vs adult | −1.44 (0.38) | 10/28 (36%) | 0.96 | fails |
| | mid/late (7) vs adult | −0.01 (0.95) | 14/28 (50%) | 0.57 | fails |
| | mid/late vs early | +1.43 (0.36) | 19/28 | 0.044 | (within-fetal) |
| H-new3 (Gao SI+LI vs Wang) | all fetal (15) vs adult | +0.31 (0.81) | 17/24 (71%) | 0.032 | fails (rule 1) |
| | early (6) vs adult | −1.18 (0.20) | 12/24 (50%) | 0.58 | fails |
| | mid/late (9) vs adult | +0.69 (0.089) | **20/24 (83%)** | **0.0008** | fails (rule 1) |
| | mid/late vs early | +1.87 (0.070) | 20/28 | 0.018 | (within-fetal) |
| Senger | — | — | — | — | not splittable (all fetal 11–22.5 W GA) |

The answers to the frozen questions follow.

1. **Does TNFRSF12A become fetal-positive with support?** No.
   - Removing the early stratum moves TNFRSF12A from negative to about zero
     (HGCA −0.46 → −0.01; Gao-original −0.08 → −0.01).
   - It reaches +0.69 (P = 0.089) only in H-new3.
   - Every early-vs-adult contrast is negative (HGCA significantly so).
2. **Does the panel fraction rise with it?** Only in H-new3 (71% → 83%). For
   HGCA it barely moves (56% → 60%) and for Gao-original not at all (50% →
   50%).

**Promotion criterion (frozen): not met for any dataset.** No mid/late
contrast qualifies.

**Reading.**
- **Developmental-window mixing is real but explains only part of the
  failure.**
  - Early fetal epithelium is TNFRSF12A-low relative to adult.
  - Mid/late is higher than early in both datasets. This is significant in
    HGCA, an independent check of the 9-week breakpoint, which was derived
    from Gao/Fawkner fetal data only.
  - Pooling therefore dilutes the signal.
- **Removing the dilution does not produce a fetal-high TNFRSF12A in HGCA or
  Gao-original:** mid/late ≈ adult.
- **The only stage-resolved contrast that approaches the H-new2 pattern is
  H-new3-late.** It has a strong panel bias (83%) but TNFRSF12A is not
  significant.

Failure types:
- **Type I (developmental-window mixing):** contributes in HGCA and Gao.
- **Type II (composition/design):** remains the explanation for Roadmap/HPA
  bulk and Visium.
- **Neither type** accounts for HGCA/Gao mid/late ≈ adult. The adult
  reference and cross-study platform effects remain candidates.

Tables:
- `results/tables/Benchmark_v2_qualification.csv` (columns
  `stage_resolution`, `parent_contrast`), the matrix, and `*_DE.csv` for
  each stage contrast.
- `HGCA_stage_resolved_units.csv` and `*_units.csv` (Gao).

## Addendum v1.3 — stage-resolved H-new2 (Fawkner fetal vs Burclaff adult)

The plan was frozen in `9571e29`. Units and model are the same as H-new2;
fetal PCW comes from the authors' sample key.

| Contrast | TNFRSF12A log2FC (P) | Panel fetal-positive | Binomial P | Verdict |
|---|---|---|---|---|
| all fetal (22) vs adult (3) | +1.83 (< 0.001) | 18/24 (75%) | 0.011 | qualifies |
| early, 8 PCW (2) vs adult | +0.86 (0.015) | 18/25 (72%) | 0.022 | qualifies |
| mid/late ≥ 9 PCW (20) vs adult | +1.94 (< 0.001) | 16/22 (73%) | 0.026 | qualifies |
| mid/late vs early | **+1.08 (0.008)** | 13/23 | 0.34 | (within-fetal) |

**Reading.**
- **The rise across the breakpoint replicates.** Mid/late fetal is higher than
  early fetal in all three fetal scRNA sources:

  | Source | Mid/late vs early, TNFRSF12A log2FC (P) |
  |---|---|
  | HGCA | +1.60 (0.009) |
  | Fawkner | +1.08 (0.008) |
  | Gao SI+LI | +1.87 (0.07) |

- **What differs is where the adult reference falls.**
  - Fawkner/Burclaff: adult sits below even the early fetal epithelium.
  - HGCA and Gao: adult is about level with mid/late fetal.
- **Caveat:** the early Fawkner arm is 2 TI units (8 PCW). ABZ2 is flagged
  (8.4% sex-discordant cells), and ABF1 and AAU2 are kept as in the parent
  contrast. See `results/tables/Hnew2_stage_resolved_units.csv`.

**Combined table.** `marker_summary/results/Fetal_stage_decomposition_summary.csv`
gives every human fetal-vs-adult parent contrast at every stage resolution,
including Senger and H-bulk1/2, which are marked not splittable.
