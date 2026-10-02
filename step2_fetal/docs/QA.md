# Step 2 QA record

Checkpoint date: 2026-10-02.

## Source integrity

- HGCA local MD5 matched the upstream MD5: `2a149b8cf04567569707e9d1fab27209`.
- SHA-256 checks passed for HGCA, all Gao inputs, mouse counts, and Ensembl mapping.
- Raw files are confined to DATA and are not tracked in Git.

## Metadata and replicate checks

- HGCA: 16 fetal and 7 adult donors after the ≥50-cell eligibility rule.
- HGCA proliferative control: 7 adult donors.
- Gao: 849 author-labelled fetal large-intestinal epithelial cells aggregated
  into 12 embryo identifiers; adult donors P1 and P2 retained separately.
- Mouse: 3 fetal and 3 adult in-vivo biological replicates.
- Mouse non-empty gene symbols: 24,587; duplicated rows: 0.

## Result assertions

- HGCA primary fetal-high: 4,733 genes.
- Complete H1: 4,645 genes.
- Gao H2: 5,271 genes.
- Mouse M1: 3,140 genes.
- Strict H1 ∩ H2 ∩ M1: 706 genes.
- Every final row was asserted to satisfy all frozen numerical and directional gates.

## Figure QA

- Human volcano: 4,733 fetal-high and 5,238 adult-high genes.
- Mouse volcano: 3,140 fetal-high and 2,270 adult-high genes.
- Two mouse genes without symbols are retained using Entrez labels.
- PDF is vector output; PNG is 600 dpi.
- Direction, thresholds, labels, clipping, and combined-panel readability were
  visually inspected at final size.

## Reviewer sanity checks

- `Mouse_log2FC > 5`: 425 genes.
- Minimum fetal mean raw CPM among these genes: 1.522.
- Median fetal mean raw CPM among these genes: 9.630.
- Genes with fetal mean raw CPM < 1: 0.
- Conclusion: the extreme-right mouse volcano cluster is supported by measurable
  fetal expression and near-zero adult expression, not a retained near-zero/near-zero artifact.
- The 706-gene HGCA–mouse effect-size scatter has descriptive Spearman rho =
  0.485; this value is not a gate.
- Requested marker labels present in the 706: GJA1 and CLU. Requested but absent:
  TACSTD2, ANXA1, TNFRSF12A, EMP1, and LAMC2.

## Known limitation

Gao is effect-only replication because the released fetal and adult matrices are
processed TPM/UMI-normalized TPM from only two adult donors. It is not used as an
FDR gate and cannot rescue failure in HGCA or mouse.
