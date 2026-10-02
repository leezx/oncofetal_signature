# Fetal epithelial state analysis (Gao + Fawkner) — plan

- Version: plan v1, 2026-10-02
- Status: **frozen before any expression value of the analysed genes was read**
  (only sample metadata, author annotations and author marker lists were read).
- Requested by review round 3 of benchmark v2. Descriptive: there is **no
  pass/fail rule and no signature**.

## Question

Not "is TNFRSF12A fetal-high?" but: **at what time, in which region, and in
which epithelial state is TNFRSF12A high in human fetal intestine?**
Specifically, is the 34-fold embryo-to-embryo TNFRSF12A variation in Gao driven
by early (6–9 week) embryos, by region, or by epithelial-state composition?

## Data (all fetal; no adult arm)

| Dataset | Cells | Unit | Age | Region |
|---|---|---|---|---|
| Fawkner-Corbett 2021 (GSE158702), 4 EPCAM+ 10x pools | QC ≥ 500 UMI, ≥ 200 genes; hashtag call as in benchmark v2 | sample = donor × region | 8–22 PCW | TI, proximal colon, distal colon, hindgut |
| Gao 2018 (GSE95630) | author-labelled epithelial cells, Tissue SI or LI | embryo × region | 6–25 W (as reported) | SI, LI |

**Fawkner hashtag → sample** (`config/sample_key_fawkner.tsv`): authors'
Mendeley supplement (doi:10.17632/gncg57p5x9.2, sheet "1. Sample Overview"),
TotalSeq hashtag number. GEO pool = run: pool 1 = run 3, pool 2 = run 4,
pool 3 = run 5, pool 4 = run 2 (matched on the PCW list in each GEO record).
Before use, each hashtag unit is checked against the key, without any analysed gene:

- **Sex**: XIST vs Y genes (RPS4Y1, DDX3Y, UTY, KDM5D) must match the key.
- **Region**: SATB2 (colon) must be higher in colon/hindgut than in TI samples
  of the same donor. This also settles the run-3 discrepancy between the
  TotalSeq and in-house hashtag numbers for ACB1/ACB3.

A unit failing either check is excluded from the age/region layers and reported.

## Epithelial state (identical procedure in both datasets)

- Reference: authors' Fawkner sub-cluster marker lists (sheets 3 and 10),
  collapsed into 9 states (`config/state_map.tsv`): Stem, TA, Fetal progenitor,
  Early enterocyte, Mature enterocyte, BEST4/OTOP2, Goblet, Paneth,
  EEC/secretory progenitor. Proximal/distal sub-clusters are merged; region is
  a separate layer.
- Marker sets: for each state, the top 30 markers per author cluster by
  avg_logFC (p_val_adj < 0.05), pooled. **All 13 analysed genes and their
  aliases are removed**, so the controls (LGR5, OLFM4, ASCL2, MKI67) stay
  independent checks of the labels.
- Cells: log-normalised; per-gene z-score; state score = mean z of the
  state's markers. Leiden clusters on the dataset's own HVG/PCA graph
  (Fawkner resolution 1.0, Gao 1.5); each cluster gets the state with the
  highest cluster-mean score.
- Gao author `Group` is not used for labels (it partly tracks age and region).
- Labels are written to disk and frozen before the analysed genes are
  tabulated.

## Genes (`config/genes.tsv`)

- Test: TNFRSF12A, TACSTD2, CLU, ANXA1.
- Controls: EPCAM, OLFM4, LGR5, ASCL2, MKI67.
- YAP/fetal context: YAP1, CCN2 (CTGF), CCN1 (CYR61), ANKRD1.

## Outputs

1. **Unit-level heatmap**: rows are units (donor × region), columns are genes.
   Colour = per-gene z-score within each dataset (platforms differ). Rows are
   annotated with dataset, donor, age and region. Values:
   - Fawkner: pseudobulk log2(CPM + 1).
   - Gao: log2(mean TPM + 1).
2. **TNFRSF12A by state × age**: one point = one unit × state pseudobulk/mean
   (minimum cells: Fawkner 10, Gao 5), never a cell. A companion panel shows
   % TNFRSF12A-positive cells (count > 0) per unit × state. The same figure
   is made for TACSTD2, CLU and ANXA1.
3. **Composition decomposition** per dataset: unit TNFRSF12A = Σ_state
   (fraction × state mean). Compare the observed age trend with:
   - a fixed-composition counterfactual (dataset-mean state fractions);
   - a fixed-expression counterfactual (dataset-mean state expression).
4. Unit × state cell counts, state fractions, label-validation table (control
   genes by state), and hashtag sex/region checks.

## Descriptive statistics only

- Spearman correlation with age is computed across units, per dataset and per
  state.
- Regions are compared within donor where a donor has more than one region.
- Age bins for display: ≤ 9, 10–13, 14–17, ≥ 18 weeks.
- No P value is used as a gate.

## Interpretation (written before results)

- **Early-embryo effect**: if TNFRSF12A is low at ≤ 9 weeks and higher later
  within the same state, embryo-to-embryo variance is mainly developmental
  timing.
- **Composition effect**: if the fixed-expression counterfactual reproduces
  the trend but the fixed-composition one does not, the variance is
  state-composition driven.
- **State program**: a state that is TNFRSF12A-high in both datasets would
  define a candidate fetal epithelial state program. Comparison with adult is
  not part of this plan.

## Addendum v1.1 — hashtag QC (2026-10-02, after sex/region checks, before any analysed gene was read)

The unit-level checks (`results/tables/Fawkner_hashtag_sample_checks.csv`)
showed sex agreement for every unit with ≥ 50 cells. Two refinements follow.

1. **The region check is applied only among units with ≥ 50 cells.** SATB2
   means of 10–20-cell units are too noisy to arbitrate.
2. **Per-cell sex discordance is added** (cell typed XX if XIST > 0 and Y
   genes = 0, XY if the reverse; discordant if the type contradicts the key).
   - Background is 2–4%.
   - Units above 10% are excluded: ABF1 (25.6%) and AAU2 (16.4%).
   - AAU2 holds 98% of pool 4's called cells, so pool 4's hashtagging is
     unreliable and the pool is not used.
   - ABZ2 (8.4%, about 3× background; the only large ≤ 9-PCW unit) is kept
     and flagged. Every Fawkner result is repeated without it.
   - Sex-discordant cells are removed from all units.

## Addendum v1.2 — state-label confidence (2026-10-02, after the label check, before any analysed test gene was read)

The forced argmax labelling from plan v1 was checked against the control genes
(`results/tables/*_state_label_check_controls.csv`) and the cluster scores
(`*_cluster_state_scores.csv`).

**Fawkner — partly validated.**
- TA is MKI67-highest and OLFM4 is absent from every fetal state.
- 7 clusters (about 5,700 cells) score ≤ 0 for every state, so their
  labels are arbitrary (e.g. "Paneth").
- Most Stem clusters lead by < 0.1.

**Gao — not validated.**
- No cluster is labelled Stem.
- The "Paneth" and "TA" clusters have low EPCAM, suggesting low-quality cells.
- Most margins are < 0.1.

**Rule added (both datasets).** A cluster keeps its state only if the
winning score is > 0.2 and it leads the runner-up by ≥ 0.1. Otherwise it is
labelled `Unresolved`, and that label is reported but not interpreted as a
state.

**Gao state layer: exploratory only.** Gao contributes to the donor, age and
region layers. Its state-level panels are shown for completeness, but
conclusions about epithelial state rest on Fawkner.
