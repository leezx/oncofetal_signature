# Fetal epithelial state analysis (Gao + Fawkner) — report

Date: 2026-10-02.

- **Plan:** [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md), frozen in `4747efe`.
- **Addenda:** v1.1 (hashtag QC, `7b898f7`) and v1.2 (label confidence, `c360264`). Both were committed before any test gene was tabulated.
- **Scope:** descriptive. No rule, no signature, no adult arm.

## Data after QC

| Dataset | Units (donor × region) | Donors | Age | Cells |
|---|---|---|---|---|
| Fawkner 2021, pools 1–3 | 20 | 11 | 8–20 PCW | 19,295 |
| Gao 2018 (≥ 20 cells per unit) | 27 | 15 embryos | 6–25 W (as reported) | 1,717 |

**Fawkner hashtag QC.** Hashtags were mapped to samples with the authors'
Mendeley sample overview. Sex agreed with the key in every unit with ≥ 50
cells.

- **Excluded:** ABF1 (25.6% sex-discordant cells) and AAU2 (16.4%). AAU2 also
  holds 98% of pool 4's cells, so pool 4 is unusable.
- **Flagged:** ABZ2 (8.4%). All Fawkner results are repeated without it.

The tables are `results/tables/Fawkner_hashtag_sample_checks.csv` and
`Unit_state_composition.csv`.

**State labels.** Labels come from the authors' marker lists, with the
analysed genes removed.

- **Fawkner (partly validated):** 58% of cells are resolved. TA has MKI67 in
  80% of cells and LGR5 in 32%; Stem has LGR5 in 21% and ASCL2 in 22%; OLFM4
  is absent from all fetal states.
- **Gao (not validated):** 54% of cells are unresolved and no Stem state is
  found. Gao states are therefore exploratory.

## Findings

### 1. Gao's TNFRSF12A instability comes from the ≤ 8-week embryos

| Gao units | n | TNFRSF12A log2(mean TPM + 1), median (range) |
|---|---|---|
| 6–8 W | 11 | 2.03 (0.00–6.50) |
| 9–25 W | 16 | 4.52 (2.46–5.51) |

- **Range:** 91-fold across all units, but only 8-fold within ≥ 9 W.
- **No age trend after 9 W:** Spearman ρ with age is 0.51 overall and 0.12
  within ≥ 9 W.
- **The early units differ in other ways too:**
  - lower EPCAM (ρ with TNFRSF12A = 0.59);
  - higher MKI67 (ρ = −0.51);
  - higher CCN1/CCN2.

  They are an earlier, less mature (or lower-quality) epithelial population.
  Composition by the exploratory Gao labels explains little of the variance
  (fixed-expression counterfactual variance 0.14 vs observed 3.67, log2).
- **Caveat:** whether Gao's "W" means post-conception or gestational weeks is
  not resolved here.

### 2. Fawkner: TNFRSF12A is flat across 8–20 PCW

- **Range:** units span log2 CPM 5.1–7.0 (3.7-fold).
- **No age trend:** ρ with age = 0.25, or 0.13 without ABZ2.
- **The two 8-PCW units are only modestly lower:** median 5.30 vs 6.15.
- **The developmental axis itself is clearly captured** by the controls over
  the same units:
  - LGR5 ρ = 0.92;
  - ASCL2 ρ = 0.66;
  - MKI67 ρ = −0.30;
  - the YAP target ANKRD1 ρ = −0.73, high only at 8 PCW.

### 3. No single state carries TNFRSF12A in Fawkner

Values are medians over unit × state pseudobulks.

| State | log2 CPM | % positive cells |
|---|---|---|
| Stem | 6.28 | 61 |
| BEST4/OTOP2 | 6.24 | 29 |
| Mature enterocyte | 5.96 | 20 |
| TA | 5.88 | 50 |
| Fetal progenitor | 5.65 | 26 |
| EEC/secretory progenitor | 5.46 | 28 |
| Goblet | 5.42 | 40 |

- **TNFRSF12A is expressed across fetal epithelial states.** Stem and TA have
  the highest detection, but % positive is partly a depth effect.
- **Stem shows a mild decline with age:** ρ = −0.61 (P = 0.046); without ABZ2,
  −0.57 (P = 0.09).
- **No 12–20-PCW progenitor state is selectively TNFRSF12A-high.**
- **Composition explains little.** Unit-to-unit variance is mostly within
  state: fixed-composition variance 0.137 vs observed 0.152; fixed-expression
  variance 0.041.

### 4. Region

- **Fawkner:** colon minus TI within donor has median +0.14 (range −1.1 to +0.9, 8 pairs), so there is no consistent region effect.
- **Gao:** LI minus SI has median +1.19 (9 of 12 positive).

### 5. The other test genes are more state-structured

- **ANXA1 marks the Fawkner Fetal progenitor state:** median 6.55 against ≤ 1
  in the other resolved states. It rises with age in Gao (ρ = 0.65).
- **TACSTD2** is highest in BEST4/OTOP2 and Goblet, and about 0.3 in mature
  enterocytes.
- **CLU** is highest in Stem and TA (8.8).

## Interpretation (for review)

- **Answer to the question.** Within 8–25 weeks, TNFRSF12A is a broadly
  expressed fetal-epithelial gene. Its level is stable from about 9 weeks, it
  is detected most often in stem/TA cells, and it is low or erratic only in
  6–8-week Gao embryos. It does not mark a discrete fetal progenitor
  sub-state.
- **Why the human contrasts disagree.** Gao's fetal arm mixes ≤ 8-week embryos
  with later ones, which plausibly explains the Gao contrasts' failure.
  - Gao-original's fetal LI arm includes 6/7/8-week embryos.
  - Fawkner (8–20 PCW, mostly ≥ 12) has no such early stratum, and H-new2 is
    built on it.
  - This is consistent with the HGCA/Gao vs H-new2 discordance, but it is not
    tested against adult here.
- **What a fetal-state program would look like.** For TNFRSF12A, a "fetal
  epithelial progenitor/state program" would be defined by developmental
  window (≥ 9–12 weeks) and stem/TA compartment, not by a unique sub-state.
  ANXA1 is the clearer state-specific marker (Fetal progenitor).
- **Open choices (not run):**
  - an adult comparison by matched state, e.g. Fawkner Stem/TA vs Burclaff
    adult Stem/TA;
  - a Gao ≥ 9-week re-analysis against adult.
  Either would be a new contrast and needs its own frozen rule.

## Outputs

- **`results/tables/`:**
  - Unit level: `Unit_gene_values.csv`, `Unit_state_gene_values.csv`.
  - Age and composition: `Age_correlations.csv`, `Early_vs_later_units.csv`,
    `Composition_decomposition{,_summary}.csv`.
  - Region: `Within_donor_region_differences.csv`.
  - Labels and QC: `Unit_state_composition.csv`, plus the label and QC
    tables.
- **`results/figures/`:**
  - `Fig1_unit_gene_heatmap` (unit × gene, with age, region and dataset).
  - `Fig2_TNFRSF12A_state_by_age` (unit × state points with % positive).
  - `Fig3_{TACSTD2,CLU,ANXA1}_state_by_age_Fawkner`.
- **Source data:** in `results/source_data/`. Cell-level objects and labels
  live in `DATA/.../2026-10-02_step2_fetal_state_v0.1/cells/`.
