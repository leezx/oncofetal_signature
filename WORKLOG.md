# Worklog — intestinal oncofetal signature

Chronological append-only record; newest entries are at the bottom.

## Progress tracker

| Step | Status | PR | Notes |
|---|---|---|---|
| 1. Literature candidate collection | completed | — | Curated evidence tables under `docs/step1_literature_candidates/` |
| 2. Conserved developmental axis | completed, PR open | [#1](https://github.com/leezx/oncofetal_signature/pull/1) | 706 strict H1 ∩ H2 ∩ M1 candidates |
| 3. CRC-high axis | not started | — | Next planned analysis |

## 2026-10-01 — Step 2 design frozen (PR #1, commit `3a6c329`)

**What**: Defined the human discovery, human replication, mouse in-vivo, adult
proliferative-control, orthology, and strict-intersection gates.

**Why**: The developmental axis had to be fixed before examining CRC expression,
preventing downstream results from changing the fetal-high definition.

**How**: HGCA uses donor-level raw-count pseudobulk and edgeR quasi-likelihood;
Gao is effect-only independent replication; GSE230581 supplies in-vivo mouse
evidence; final mapping is Ensembl one-to-one only.

**Review**: User-approved after simplifying the design to the minimum frozen
gates and explicitly excluding organoid data as primary evidence.

## 2026-10-01 — Step 2 executed (PR #1, commit `8363953`)

**What**: Downloaded and checked source data, generated pseudobulks, ran H1/H2/M1,
mapped Ensembl orthologues, and wrote the strict conserved set.

**How**: HGCA included 16 fetal and 7 adult donors; Gao retained 12 fetal embryos
and both adult donors; mouse included 3 fetal and 3 adult biological replicates.
All raw/processed data remained under `DATA`; only version-controlled results
were added to Git.

**Real findings**: 4,645 genes passed H1, 5,271 passed H2, 3,140 passed M1, and
706 passed the strict H1 ∩ H2 ∩ M1 intersection. Download and result checksums,
script syntax, and all frozen-gate assertions passed.

**Review**: The Gao corrected supplementary annotation resolved the only material
metadata concern. No integration, trajectory, GSEA, or organoid substitution was added.

## 2026-10-02 — DEG figures and reproducibility checkpoint (PR #1)

**What**: Added complete human/mouse DEG tables, individual and paired volcano
plots, plotting source data, a normalized Step 2 directory, executable runbook,
rendered report, and module inventory.

**How**: Volcano plots use the frozen thresholds (`|log2FC| >= 0.5`,
`FDR < 0.05`) and an R-only rendering workflow. PDF vector and 600-dpi PNG
outputs are retained. Two mouse rows lacking symbols remain represented by
Entrez IDs in plotting source data.

**Real findings**: Human HGCA contains 4,733 fetal-high and 5,238 adult-high DEGs;
mouse contains 3,140 fetal-high and 2,270 adult-high DEGs under the symmetric
volcano classification. Visual QA and source-data assertions passed.

**Review**: Repository visibility was changed from private to public at the
user's request. This checkpoint corrects stale README status and DATA checksum
pointers while preserving unrelated user files.

## 2026-10-02 — Step 2 reviewer sanity checks (PR #1)

**What**: Implemented only the three requested checks: mouse symbol duplication,
raw per-sample CPM for mouse genes with log2FC > 5, and a 706-gene human–mouse
effect-size scatter.

**How**: The mouse module now stops if a duplicated non-empty symbol is detected,
exports raw CPM without TMM-adjusted library sizes for the six biological
samples, and flags fetal mean CPM below 1 descriptively. The scatter uses the
frozen 706-gene universe and does not add a correlation gate.

**Real findings**: There were zero duplicate symbols. All 425 genes with mouse
log2FC > 5 had fetal mean raw CPM ≥ 1.522; none met the low-count flag. The
706-gene effect sizes showed descriptive Spearman rho = 0.485. Only GJA1 and CLU
of the seven requested labels belong to the final 706.

**Review**: The sanity checks support the existing mouse volcano and do not
justify changing thresholds or re-running Step 2 with a more complex model.
