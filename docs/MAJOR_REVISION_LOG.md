# Major revision log

## 2026-10-02 — Step 2 developmental dataset suitability failure

### Severity and status

- Severity: **major / analysis-invalidating**
- Affected step: Step 2, conserved developmental axis
- Current disposition: **do not freeze Step 2 and do not begin Step 3**
- Threshold status: unchanged; this is not a threshold-tuning decision
- Data status: all existing inputs and outputs are retained for audit, but the
  706-gene intersection must not be treated as a validated developmental set

### Trigger

The gate-by-gate audit of 31 literature-curated positive intestinal oncofetal
markers exposed a systematic conflict between prior biological evidence and the
selected dataset contrasts. Several expected positive markers were estimated as
adult-high rather than fetal-high, and many others were absent from the tested
HGCA gene universe. This pattern was not visible from global DEG counts or from
the final selected-gene intersection alone.

### Observed evidence

The frozen audit table is
`step2_fetal/results/tables/Literature_31_gate_failure_audit.csv`.

- Only 5 of 31 candidates passed the final H1 ∩ H2 ∩ M1 definition:
  `GJA1`, `CLU`, `ANXA6`, `RBP1`, and `BASP1`.
- HGCA tested 23 of the 31 candidates. Of those, 11 had negative fetal/adult
  log2FC: `ANXA1`, `IL33`, `TNFRSF12A`, `AREG`, `CD44`, `EDN1`, `IL1RN`,
  `MSLN`, `EMP1`, `LAMC2`, and `EPS8L1`.
- Eight candidates were not represented in the HGCA tested output:
  `LY6A`, `SPRR1A`, `CCN2`, `CCN1`, `SOX17`, `REG3B`, `MMP7`, and `L1CAM`.
- Fifteen candidates first failed at the HGCA primary gate; eight were not
  tested by HGCA; only three first failed at a later gate.
- Gao tested 27 candidates, of which 14 had negative log2FC. Mouse tested 19,
  of which 4 had negative log2FC.
- Illustrative contradictions include `ANXA1` (HGCA log2FC -2.045),
  `TNFRSF12A` (-0.460), `AREG` (-3.739), `CD44` (-4.215), `EMP1` (-5.636),
  and `LAMC2` (-1.468).
- `TACSTD2` was positive but did not pass HGCA primary (log2FC 0.368,
  FDR 0.572), despite passing both Gao H2 and mouse M1.

The genome-wide cross-species diagnostic also showed only modest effect-size
concordance (11,049 tested one-to-one genes; Spearman rho 0.222), although H1
and M1 calls were enriched for overlap (998 genes; Fisher OR 2.54,
P = 2.4e-79). Enrichment alone does not resolve the literature-marker conflict.

### Interpretation

The computational pipeline is reproducible and the frozen gates were applied as
specified. The newly identified problem is the **fitness of the selected data
and biological contrasts for the intended oncofetal definition**, not an
identified sign inversion or software failure.

The current comparisons are not biologically equivalent:

- HGCA combines heterogeneous fetal and healthy-adult intestinal epithelial
  compartments.
- Gao uses a different platform and large-intestinal sampling scheme.
- Mouse compares E16.5 proximal small-intestinal epithelium with adult proximal
  crypt epithelium.

Differences in anatomical region, developmental time, epithelial composition,
cell-state annotation, and tested-gene coverage can dominate the intended
developmental signal. At present, no single cause has been proven; dataset and
metadata suitability must be re-audited before any replacement analysis.

### Consequences

1. The 706-gene set is reclassified as a **provisional output from an unsuitable
   or insufficiently validated contrast**, not the final conserved fetal-high
   universe.
2. The volcano plots, DEG tables, orthology table, and diagnostic statistics are
   retained as reproducibility records, but must not be used for biological
   claims or CRC-high filtering.
3. PR #1 remains an audit trail rather than approval to proceed to Step 3.
4. No gate will be relaxed to rescue literature markers. The failure must be
   addressed at the dataset, metadata, annotation, or contrast-definition level.

### Required remediation before rerun

1. Re-audit HGCA metadata and tested-gene coverage, including the exact fetal
   and adult regions, epithelial subtypes, donor composition, gene identifiers,
   and reasons for the eight missing candidates.
2. Verify expression direction for a prespecified positive-control panel before
   genome-wide DE. At minimum include `TACSTD2`, `ANXA1`, `TNFRSF12A`, `GJA1`,
   `CLU`, `AREG`, `CD44`, `EMP1`, and `LAMC2`.
3. Decide whether the current whole-epithelium/crypt contrasts answer the same
   biological question. If not, select better-matched in-vivo datasets or
   redefine comparable epithelial compartments before running edgeR.
4. Predefine dataset-level acceptance criteria. A replacement dataset must show
   expected direction for the positive-control panel without using those genes
   to tune differential-expression thresholds.
5. Rerun Step 2 under a new versioned analysis directory. Do not overwrite the
   current files; preserve them as the failed-analysis checkpoint.

### Reopening criterion

Step 2 can be reconsidered only after the dataset/contrast audit is documented,
the positive-control discrepancies are explained or resolved, and a new
versioned run passes both computational QA and biological positive-control QA.
