# Replacement dataset qualification report

Date: 2026-10-02

## Decision summary

The four proposed datasets do not support a single automatic replacement of
the invalidated Step 2 design.

- **Pikkupeura GSE160449 qualifies as a strong mouse culture-state benchmark.**
  It recovers positive fetal direction for 25 of 27 measured literature
  candidates and all five priority markers. It remains culture-based and cannot
  alone replace an in-vivo developmental criterion.
- **GSE44433 is informative but does not pass as a clean universal positive
  benchmark.** It recovers 10 of 23 measured candidates and four of five
  priority markers; `Tacstd2` is slightly negative. It is useful independent
  in-vivo evidence, not a sole discovery dataset.
- **Senger GSE101531 fails the proposed positive-marker direction benchmark.**
  Only 14 of 27 measured candidates and one of five priority markers are
  fetal-positive. Culture matching reduces composition confounding but does not
  make the contrast a faithful positive-control dataset for this marker panel.
- **Fawkner-Corbett GSE158702 is not eligible for fetal/adult FC.** It has no
  matched adult arm. It confirms fetal epithelial expression/detection only and
  should be used to localize fetal cell states after an annotated object is
  available.

No new hard gate or conserved signature is defined by this benchmark.

## Priority-marker results

| Gene | Senger human log2FC | Fawkner fetal detection | Pikkupeura LN log2FC | Pikkupeura collagen log2FC | GSE44433 log2FC |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | -1.054 | 0.061 | 9.440 | 12.052 | -0.168 |
| GJA1 | 2.510 | 0.074 | 10.088 | 9.456 | 0.691 |
| CLU | -0.499 | 0.537 | 8.495 | 8.304 | 0.479 |
| ANXA1 | -1.772 | 0.039 | 7.382 | 8.097 | 2.396 |
| TNFRSF12A | -0.217 | 0.301 | 1.785 | 2.126 | 0.309 |

Positive log2FC denotes fetal-high. Fawkner values are fractions of retained
fetal EpCAM-positive cells with non-zero expression and are not fold changes.

## Dataset-level recovery

| Dataset | Measured | Positive/detected | Priority measured | Priority positive/detected |
|---|---:|---:|---:|---:|
| Human Senger | 27 | 14 | 5 | 1 |
| Human Fawkner | 31 | 26 with median pool detection > 0 | 5 | 5 with median pool detection > 0 |
| Mouse Pikkupeura | 27 | 25 | 5 | 5 |
| Mouse GSE44433 | 23 | 10 | 5 | 4 |

Pikkupeura arm directions disagree for `Areg` and `Mmp7`; `L1cam` is negative
in both arms. GSE44433 includes five WT fetal and five WT adult samples only;
all TnfΔARE and maternal-genotype comparison groups were excluded.

## Quality assessment

### High-severity constraints

1. Fawkner lacks an adult reference, so any Fawkner fetal/adult FC would be a
   fabricated cross-study comparison.
2. Senger and Pikkupeura are culture-state experiments. Their results may
   capture developmental identity but also matrix, organoid, or YAP-associated
   culture effects.
3. GSE44433 is an older microarray with incomplete marker coverage and contains
   multiple genotype/maternal groups; only the WT-to-WT subset is valid here.

### Reproducibility safeguards

- Raw data and checksums remain under DATA, not Git.
- Candidate genes were prespecified before inspection.
- Human-to-mouse mapping uses Ensembl release 116 one-to-one orthologues plus
  the independent NCBI mouse `gene_info` table; it does not depend on HGCA
  tested-gene coverage.
- Fawkner raw-feature matrices use the same fixed cell filter in all four pools:
  at least 500 total UMI and 200 detected genes.
- No marker was used to tune a DE threshold.

## Review recommendation

Do not simply replace HGCA with Senger. Promote Pikkupeura only as a mouse
culture benchmark and retain GSE44433 as independent in-vivo replication. Human
discovery still requires a dataset with a genuine, metadata-compatible adult
epithelial reference. Fawkner should remain fetal-state annotation evidence
unless such an adult reference is prospectively specified and justified.
