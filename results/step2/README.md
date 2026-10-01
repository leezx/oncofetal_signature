# Step 2 results — developmental axis

Frozen analysis version: `2026-10-01_step2_developmental_axis_v0.1`

## Outcome

The strict intersection contains **706 Conserved Fetal-High Candidates**.

| Gate | Tested/mapped genes | Passing genes |
|---|---:|---:|
| HGCA fetal vs adult primary | 17,779 | 4,733 |
| HGCA H1 including proliferative control | 17,779 | 4,645 |
| Gao H2 effect-only validation | 24,153 | 5,271 |
| Mouse M1 in vivo fetal vs adult | 13,314 | 3,140 |
| H1 ∩ H2 ∩ M1 after Ensembl one-to-one mapping | 13,397 mapped | 706 |

## Files

- `HGCA_DE.csv`: HGCA primary and proliferative-control evidence joined by gene.
- `HGCA_proliferative_control.csv`: the standalone negative-control contrast.
- `Gao_validation.csv`: fetal and both adult-donor expression values plus the H2 decision.
- `Mouse_DE.csv`: in vivo E16.5 fetal epithelium versus adult crypt epithelium.
- `Cross_species_evidence.csv`: all mapped genes and all three gate decisions.
- `Conserved_Fetal_High.csv`: strict H1 ∩ H2 ∩ M1 candidates only.

Raw and processed data are intentionally excluded from Git and are managed under
`/Volumes/Stelligen_SSD/Stelligen/DATA`. See `docs/STEP2_EXECUTION_REPORT.md` for
the full audit trail and exact data locations.
