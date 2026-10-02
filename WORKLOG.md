# Worklog — intestinal oncofetal signature

Chronological append-only record; newest entries are at the bottom.

## Progress tracker

| Step | Status | PR | Notes |
|---|---|---|---|
| 1. Literature candidate collection | completed | — | Curated evidence tables under `docs/step1_literature_candidates/` |
| 2. Conserved developmental axis | **invalidated; dataset review required** | [#1](https://github.com/leezx/oncofetal_signature/pull/1) | 706-gene output is provisional and must not be used downstream |
| 3. CRC-high axis | **in progress — T1 + S2 done; S1 (Joanito) blocked on Synapse terms acceptance** | branch `plan/step3-cancer-axis` | Run genome-wide and independently of Step 2 at user instruction; fetal ∩ CRC intersection deferred |

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

## 2026-10-02 — Cross-species conservation and candidate gate audit (PR #1)

**What**: Added the two requested non-gating diagnostics: a complete 31-candidate
gate-failure table and a genome-wide human–mouse comparison across every
one-to-one gene tested in both datasets.

**How**: Module 08 joins frozen H1 and M1 calls without using Gao to define the
tested universe, computes Spearman correlation, constructs the H1×M1 quadrant
table, runs a one-sided Fisher enrichment test, and exports exact figure source
data. No threshold or model was changed.

**Real findings**: Among 11,049 tested one-to-one genes, rho = 0.222. H1 and M1
overlapped at 998 genes (OR = 2.54, P = 2.4e-79); 37.2% of H1 genes passed M1
and 38.7% of M1 genes passed H1. TACSTD2, ANXA1, and TNFRSF12A all first failed
HGCA primary, not Gao. TACSTD2 still passed Gao and mouse, ANXA1 passed mouse,
and TNFRSF12A passed neither downstream gate.

**Review**: The former rho = 0.485 is now explicitly documented as conditioned
on the final 706-gene selection and is not presented as genome-wide conservation.

## 2026-10-02 — MAJOR REVISION: developmental dataset suitability failure

**Trigger**: The 31-candidate gate audit revealed a systematic contradiction
between the selected developmental contrasts and literature-curated positive
oncofetal markers. Only five candidates passed the final intersection. In HGCA,
11 of 23 tested candidates had negative fetal/adult log2FC and eight candidates
were absent from the tested output; 15 first failed HGCA primary.

**Assessment**: The workflow is computationally reproducible, but reproducibility
does not establish biological validity. The selected datasets and contrasts are
not sufficiently fit for the intended conserved fetal-high definition. The
precise cause—region, developmental stage, epithelial composition, annotation,
gene coverage, or a combination—has not yet been isolated.

**Decision**: Step 2 is no longer complete. The 706-gene result is provisional,
must not be used for biological claims or Step 3, and will not be rescued by
changing thresholds. Existing outputs remain frozen as an audit trail.

**Remediation**: Re-audit dataset metadata, gene coverage, epithelial subsets,
and a prespecified positive-control panel before choosing or rerunning the
developmental contrasts. Any replacement analysis must use a new versioned
directory and preserve the current failed-analysis checkpoint.

**Full record**: See `docs/MAJOR_REVISION_LOG.md`.

## 2026-10-02 — Replacement dataset qualification benchmark

**Scope**: Evaluated the 31 prespecified literature markers in four proposed
replacement resources before permitting any new hard developmental gate.
No genome-wide signature was constructed.

**Contrasts**: GSE101531 used six fetal versus three adult human epithelial
enterospheres. GSE158702 was restricted to fetal EpCAM-positive expression
validation because it has no adult arm. GSE160449 evaluated fetal versus adult
mouse cultures separately in LN and collagen arms. GSE44433 used only five WT
E17.5 and five WT eight-week laser-microdissected ileal epithelial samples.

**Findings**: Pikkupeura recovered positive fetal direction for 25/27 measured
candidates and all five priority markers. GSE44433 recovered 10/23 and four of
five priority markers, with Tacstd2 slightly negative. Senger recovered only
14/27 and one of five priority markers. Fawkner confirmed fetal expression but
cannot provide a fetal/adult fold change.

**Decision for review**: Pikkupeura qualifies as a strong mouse culture
benchmark; GSE44433 is useful independent in-vivo replication but not a sole
discovery dataset; Senger fails this marker-direction qualification; Fawkner is
fetal-state annotation evidence only. A suitable human discovery dataset with a
matched adult epithelial reference remains unresolved. Step 3 remains blocked.

**Audit**: See `step2_dataset_benchmark/docs/BENCHMARK_REPORT.md` and
`step2_dataset_benchmark/results/tables/`.

## 2026-10-02 — Step 3 started at user instruction (independent of Step 2)

**Decision**: The user instructed that Step 3 begin now although Step 2 is
invalidated. Step 3 is therefore run genome-wide and does not use the
provisional 706-gene output in any way. The fetal ∩ CRC intersection is
deferred until a replacement Step 2 passes biological QA, and Step 3 results
must not be used to choose the replacement Step 2 datasets or thresholds.

**Plan**: `step3_cancer/docs/ANALYSIS_PLAN.md`, committed (`a2db19f`) before
any tumour/normal result was examined. CRC-high = T1 (TCGA COAD+READ primary
tumour vs normal, log2FC ≥ 0.5, FDR < 0.05, `~ tissue + project`) ∩ S1
(Joanito malignant vs normal epithelium, same thresholds, plus a malignant >
normal stem/TA proliferation control) ∩ S2 (Pelka tumour vs normal epithelium,
log2FC > 0, FDR < 0.05). New relative to Step 2: an oncofetal-independent
dataset admission QA (8 canonical tumour-up and 10 tumour-down CRC markers;
≥ 80% in the expected direction per gated contrast).

**Executed**:
- recount3 TCGA COAD/READ + GTEx COLON (uniform Monorail processing). After
  excluding FFPE, non-01/11 sample types, and duplicate aliquots: 624 primary
  tumours (458 COAD, 166 READ) and 51 normals. T1: 28,204 genes tested,
  10,266 tumour-high, 6,011 normal-high (all GENCODE gene types). Paired
  sensitivity: 50 pairs. GTEx sensitivity: 433 transverse colon samples
  (direction only; source fully confounded).
- Pelka GSE178341: 62 tumour and 35 eligible normal patient pseudobulks
  (one normal pseudobulk < 50 cells excluded); 14,696 genes tested.
- Admission QA: T1 and S2 recovered 18/18 panel genes in the expected
  direction. The Pelka proliferation-control contrast is not a gated contrast
  and is reported as not applicable.
- Bulk vs epithelial effect sizes: Spearman rho = 0.702 over 14,460 jointly
  tested genes.

**Descriptive candidate notes (not gates)**: TACSTD2, TNFRSF12A, SPP1, MMP7,
RBP1, IL1RN and LAMC2 pass both T1 and S2. CLU and ANXA6 are normal-high in
bulk but tumour-high in Pelka epithelium; ANXA1 is flat in bulk but tumour-high
in epithelium — candidates for the `Epithelial_only` label once S1 exists.
EMP1 is normal-high in both bulk and Pelka whole-epithelium, yet higher in
tumour than normal stem/TA cells. CCN1/CCN2 were matched through their
GENCODE v26 symbols CYR61/CTGF (`config/symbol_aliases.tsv`); LY6A and REG3B
have no human annotation; GJA1, SPRR1A, ANKRD1, SOX17 and L1CAM were removed
by Pelka expression filtering.

**Blocker**: Joanito Synapse syn26844071 has a self-sign access requirement
(9606933) that the account owner must accept in the Synapse web UI. No
`CRC_high` call is made until S1 runs; current labels are
`T1_pass_S1_pending` / `Not_T1`.
