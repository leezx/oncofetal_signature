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
