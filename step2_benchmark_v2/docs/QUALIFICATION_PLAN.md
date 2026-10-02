# Step 2 human benchmark v2 — qualification plan

- Version: step2_benchmark_v2 plan v1
- Date: 2026-10-02
- Status: **frozen before any new dataset was opened** (only GEO series and
  sample-level metadata were read).
- Scope: dataset qualification only. **No fetal-high signature, genome-wide
  intersection, or threshold for gene selection is produced.**

## Purpose

The original Step 2 human contrasts (HGCA, Gao + GSE103154) failed the
31-marker biological QA. This benchmark adds three human fetal-vs-adult
contrasts and asks, for every contrast, whether the 31 Step 1 positive markers
are recovered — and, as a debugging experiment, whether a failure is caused by
the adult reference or by the fetal compartment.

## Contrasts

Positive log2FC always means higher in fetal epithelium.

| ID | Fetal | Adult | Unit | Statistic |
|---|---|---|---|---|
| H-new1 | Fawkner-Corbett 2021 Visium (GSE158328): fetal colon 12/19 PCW + fetal small intestine 12 PCW sections | same study: adult colon sections | section | edgeR QL on raw-count epithelial-spot pseudobulk, `~ stage` |
| H-new1-colon | as H-new1, fetal colon sections only | adult colon sections | section | sensitivity, same model |
| H-new2 | Fawkner-Corbett 2021 EPCAM+ scRNA (GSE158702), hashtag-demultiplexed samples | Burclaff 2022 (GSE185224) healthy adult epithelium, all regions | fetal sample / adult donor | edgeR QL on raw counts, `~ stage` (study confounded with stage) |
| H-new3 | Gao 2018 fetal SI + LI epithelium (GSE95630) | Wang 2020 (GSE125970) adult ileum + colon + rectum epithelium | embryo / adult sample | effect-only, cross-platform (see below) |
| H-new3-debug | Gao fetal **LI** epithelium (identical to Gao-original) | Wang adult colon + rectum | embryo / adult sample | isolates the adult reference: only the adult arm differs from Gao-original |
| Gao-original (recomputed) | Gao fetal LI epithelium | GSE103154 adult LI (P1, P2) | embryo / adult donor | same effect-only machinery as H-new3 |

Existing contrasts carried into the matrix unchanged: HGCA (H1), Senger
(GSE101531), Pikkupeura LN and collagen arms (GSE160449, mouse).

### Epithelial selection (no candidate marker is used)

- **Visium spots (H-new1)**: epithelial score = fraction of spot UMIs from
  `EPCAM, CDH1, KRT8, KRT18, KRT19, KRT20, CLDN7, VIL1`. Within each section,
  spots in the top quartile of this score that also have `EPCAM` > 0 are
  epithelial-rich spots. Stromal carry-over is reported (fraction of UMIs from
  `VIM, COL1A1, COL3A1, DCN, ACTA2, PTPRC`) per section and is not corrected.
- **Fawkner scRNA (H-new2)**: all cells of the four EPCAM+ epithelial pools
  (EPI1–3, EPI pool 4); cells assigned to a single hashtag (highest HTO count
  ≥ 2× the second and ≥ 10 counts) form one fetal sample; negatives and
  doublets are dropped. QC: ≥ 500 UMIs and ≥ 200 genes per cell.
- **Burclaff (H-new2)**: all cells of the authors' annotated epithelial object;
  one pseudobulk per donor.
- **Wang (H-new3)**: all annotated epithelial cells; one pseudobulk per sample.
- **Gao fetal (H-new3)**: author-labelled fetal small- and large-intestinal
  epithelial cells (oesophagus/stomach excluded), aggregated per embryo, using
  the same author labels as the original Step 2 H2.
- Eligibility: ≥ 50 cells (scRNA) or ≥ 20 spots (Visium) per unit.

### Cross-platform effect-only statistic (H-new3, H-new3-debug, Gao-original)

Gao fetal values are TPM; adult values are UMI-based. Each unit is the mean
per-cell expression on a per-million scale (TPM for Gao; UMI CPM for Wang and
GSE103154). Effect: `log2((median fetal unit + 1) / (median adult unit + 1))`.
Statistical support: Welch t-test on `log2(unit + 1)` between fetal and adult
units. Platform is completely confounded with stage; the P value is descriptive.

## Frozen qualification rule (applied identically to every contrast)

1. **Biological positive control**: TNFRSF12A must be fetal-high with
   statistical support — log2FC > 0 **and** P < 0.05 (single prespecified
   gene, nominal P).
2. **Panel bias**: among the 31 markers measured in the contrast, the number
   with log2FC > 0 must exceed 50% by a one-sided exact binomial test,
   P < 0.05.

A contrast **qualifies** only if both hold. TACSTD2, CLU and ANXA1 are reported
at the top of every table but are not individually required. Markers without a
human orthologue (LY6A, REG3B) are not counted in human contrasts.

## Interpretation rules for the Gao debug pair

| Gao-original | H-new3-debug | Reading |
|---|---|---|
| TNFRSF12A ≤ 0 | TNFRSF12A > 0, supported | adult reference GSE103154 is the problem |
| TNFRSF12A ≤ 0 | TNFRSF12A ≤ 0 | fetal compartment/stage/composition is the problem |

If, after this benchmark, TNFRSF12A is still fetal ≤ adult in every human
contrast, the automated pipeline stops and TNFRSF12A is examined manually by
donor × gestational age × region × epithelial subtype.

## Outputs

- `results/tables/Benchmark_v2_31_marker_matrix.csv` — one row per marker
  (TACSTD2, CLU, ANXA1, TNFRSF12A first), log2FC / P / FDR per contrast.
- `results/tables/Benchmark_v2_qualification.csv` — rule 1, rule 2, verdict per contrast.
- `results/tables/*_units.csv` — per-unit inclusion (cells/spots, library size).
- `results/tables/*_marker_unit_values.csv` — per-unit values of the 31 markers
  (donor-level consistency).
- Raw and intermediate data stay under `DATA`; this directory holds code,
  tables, and figures only.

## Addendum v1.1 — primary-tissue bulk contrasts (frozen 2026-10-02, before data opened)

Requested after the v2 review. Source tracing (done before any expression value
was read) established that both proposed bulk contrasts draw on the **same
public samples**:

- Finkbeiner et al. 2015 (Stem Cell Reports) reprocessed 6 Roadmap
  Epigenomics fetal small-intestine total-RNA RNA-seq samples (GEO GSE18927 /
  SRA SRP001371: GSM1059486, GSM1059507, GSM1059508, GSM1059517, GSM1059519,
  GSM1059521; 91–120 days) and 6 adult small-intestine samples from the Human
  Protein Atlas (ArrayExpress E-MTAB-1733 / ENA ERP003613: 2 duodenum + 4 small
  intestine). Their repository shares Cufflinks FPKM only.
- Senger et al. 2018 Table 2 ("scraped mucosae") lists the same Roadmap fetal
  samples and the two E-MTAB-1733 duodenum samples (and repeats GSM1059508 for
  two ages — an accession inconsistency in the paper).

Therefore H-bulk2 is a **duodenum-only subset of H-bulk1**, not an independent
dataset. Both are independent of Fawkner/Burclaff (H-new2).

| ID | Fetal | Adult | Model |
|---|---|---|---|
| H-bulk1 | 6 Roadmap fetal SI (SRP001371) | all HPA duodenum + small-intestine samples (ERP003613) | edgeR QL, `~ stage`, recount3 raw counts |
| H-bulk2 | the same 6 fetal SI | HPA duodenum samples only | same; minimum group size 2 |

- Data: recount3 (Monorail, GENCODE v26) gene counts for both SRA projects,
  i.e. the same uniform pipeline as Step 3 TCGA/GTEx. Authors' FPKMs are not
  used for statistics.
- Known confounds, reported not corrected: different laboratories; fetal
  libraries are total RNA (ribo-depleted), HPA libraries are poly(A); whole
  tissue (not epithelium-only) in both arms.
- Author sanity control (reported only): OLFM4 fetal < adult.
- Same frozen qualification rule (TNFRSF12A fetal-high P < 0.05; panel
  binomial P < 0.05). Tier labels for reporting: Tier A epithelial scRNA
  (H-new2); Tier B primary-tissue bulk (H-bulk1, H-bulk2); Tier C problematic
  references (HGCA, Gao, Senger enterospheres); Visium dropped.

## Addendum v1.2 — stage-resolved re-analysis (frozen 2026-10-02, before any stage-split contrast was computed)

**Breakpoint.** 9 weeks. It comes from the fetal-only manual analysis in
`step2_fetal_state/` (Gao TNFRSF12A: ≤ 8 W low and erratic, plateau from
9 W). No adult data were used to choose it. The breakpoint is applied to the
age as each source reports it:

- **HGCA:** post-conception weeks (e.g. "8.4Wk").
- **Gao:** "W", convention not stated (possibly gestational).
- **Senger:** "weeks gestational age".

| Stage | Definition |
|---|---|
| Early fetal | reported age < 9 weeks |
| Mid/late fetal | reported age ≥ 9 weeks |
| Adult | unchanged |

**Contrasts, per dataset.** These use the same units, filters and statistics
as the original contrast; only the fetal arm is subset.

1. All fetal vs adult (the original, recomputed).
2. Early fetal vs adult.
3. Mid/late fetal vs adult.
4. Mid/late fetal vs early fetal (within-fetal, same platform).

**Datasets.**

| Dataset | Fetal age distribution | Early / late units | Statistic |
|---|---|---|---|
| HGCA (H1) | 6.1–17 PCW | 6 donors (6.1–8.4) / 10 donors (9.2–17); adult 7 | edgeR QL `~ stage`, donor pseudobulk |
| Gao-original (Gao fetal LI vs GSE103154) | 6–25 W | 5 / 7 embryos; adult 2 | effect-only (module 04) |
| H-new3 (Gao fetal SI + LI vs Wang) | 6–25 W | 6 / 9 embryos; adult Wang samples | effect-only (module 04) |
| Senger enterospheres | 11–22.5 W GA | 0 / 6 | **not splittable**: no early stratum, so mid/late vs adult is identical to the original |
| H-bulk1/2 (Roadmap) | 13–17 W | 0 / 6 | not splittable (already all ≥ 9 W) |

H-new3-debug is not split (it is a debug contrast only). For the within-fetal
contrast 4, Gao uses the same effect-only statistic (no platform
confounding); HGCA uses edgeR.

**Readouts.** Every stage contrast receives the unchanged frozen rule (rule 1
TNFRSF12A; rule 2 panel binomial). Two questions are reported:

1. Does TNFRSF12A become fetal-positive with statistical support?
2. Does the fetal-positive fraction of the whole 31-marker panel rise at the
   same time?

**Promotion criterion.** A dataset is a candidate for promotion back into
primary human developmental evidence only if its mid/late-vs-adult contrast
qualifies under the frozen rule **and** its panel fraction rises compared with
all-fetal. A TNFRSF12A-only reversal is not enough. Promotion remains a
review decision.
