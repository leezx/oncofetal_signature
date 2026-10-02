# Dataset qualification plan

## Objective

Determine whether four proposed replacement datasets reproduce the direction
and fetal expression of 31 prespecified literature markers before selecting a
new developmental discovery dataset.

## Primary output

One row per literature gene with Senger human log2FC, Fawkner fetal-expression
metrics, Pikkupeura mouse log2FC, and GSE44433 mouse log2FC. The five priority
markers are `TACSTD2`, `GJA1`, `CLU`, `ANXA1`, and `TNFRSF12A`.

## Guardrails

1. No new fetal-high signature or genome-wide intersection is generated.
2. Marker directions are diagnostic; they are not used to tune statistical
   thresholds.
3. Fawkner-Corbett has no matched adult arm and therefore cannot yield an FC.
4. Organoid/culture comparisons are labelled as such and cannot by themselves
   establish an in-vivo developmental criterion.
5. GSE44433 uses only WT offspring of WT dams at both ages.
6. Dataset qualification requires metadata, gene-coverage, sample-grain, and
   positive-control checks; computational reproducibility alone is insufficient.

## Review decision

Results must be reviewed before any dataset is promoted to discovery, primary
validation, or hard-gate status.
