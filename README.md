# oncofetal_signature

Construction of an evidence-tiered **intestinal oncofetal reference signature**
(Core-OF / Extended-OF / Intestinal-Fetal / CRC-Oncofetal), built from first
principles (fetal vs adult intestine in mouse + human, CRC malignant vs normal
epithelium, tissue-specificity vs other fetal organs) and benchmarked against
published fetal / regenerative / revSC / revCSC signatures. Methodological
foundation for the first half of the TWEAKR paper.

**Status:** Step 2 is **invalidated pending dataset/contrast review**. The
previous 706-gene result is retained as a failed-analysis checkpoint and must
not be used for Step 3. CRC-high analysis is blocked.

- Founding task: [`docs/TASK_BRIEF.md`](docs/TASK_BRIEF.md)
- Current checkpoint: [`WORKLOG.md`](WORKLOG.md)
- Major revision record: [`docs/MAJOR_REVISION_LOG.md`](docs/MAJOR_REVISION_LOG.md)
- Step 2 analysis package: [`step2_fetal/`](step2_fetal/)
- Module inventory: [`docs/MODULE_INDEX.md`](docs/MODULE_INDEX.md)

## Data

No raw expression matrices or pseudobulk intermediates live in this repo. They
are managed under `/Volumes/Stelligen_SSD/Stelligen/DATA/`. Version-controlled
DE tables, figure source data, and figures live in `step2_fetal/results/`.

## Repository layout

```text
docs/                       project-level documentation
step2_fetal/
  config/                   frozen paths, thresholds, and Ensembl query
  docs/                     plan, runbook, execution report, QA report
  scripts/                  ordered acquisition/analysis/plotting code
  results/
    tables/                 DEG and cross-species evidence tables
    source_data/            plotting source data
    figures/                PDF and 600-dpi PNG figures
WORKLOG.md                  chronological project record
```
