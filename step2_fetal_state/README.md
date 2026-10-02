# step2_fetal_state — manual fetal epithelial state analysis (Gao + Fawkner)

This module asks when, where and in which epithelial state TNFRSF12A is high
in human fetal intestine. It was requested by benchmark v2 review round 3.
- Descriptive only: no signature and no adult arm.
- Plan: [`docs/ANALYSIS_PLAN.md`](docs/ANALYSIS_PLAN.md) (frozen; addenda v1.1–v1.2).
- Results: [`docs/REPORT.md`](docs/REPORT.md).
- Run: `bash scripts/run_fetal_state.sh`. Scripts 01–03 are Python and 04 is R.
- Inputs live under `DATA/scRNAseq/GSE158702_FawknerCorbett2021` (including
  the authors' Mendeley supplement) and `DATA/scRNAseq/GSE103239_Gao2018`.
