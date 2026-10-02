# Worklog — intestinal oncofetal signature

Chronological append-only record; newest entries are at the bottom.

## Progress tracker

| Step | Status | PR | Notes |
|---|---|---|---|
| 1. Literature candidate collection | completed | — | Curated evidence tables under `docs/step1_literature_candidates/` |
| 2. Conserved developmental axis | **invalidated; dataset review required** | [#1](https://github.com/leezx/oncofetal_signature/pull/1) | 706-gene output is provisional and must not be used downstream |
| 2b. Human benchmark v2 | **completed; in review** | [#3](https://github.com/leezx/oncofetal_signature/pull/3) | H-new2 adopted as cross-study human primary (needs independent validation); Visium dropped; no signature built |
| 3. CRC-high axis | **executed (plan v3); in review** | [#2](https://github.com/leezx/oncofetal_signature/pull/2) | CRC_high = Joanito S1 ∩ Pelka S2 passing P; Joanito-derived results restricted (not in git) |

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

## 2026-10-02 — Step 3 plan v2 after user review

**Review verdict**: "Approve with one conceptual modification: do not let TCGA
bulk veto a malignant-epithelial CRC-high gene."

**Change**: `CRC_high` = S1 ∩ S2 (Joanito primary + Pelka replication). TCGA T1
is reported per gene as `bulk_support` and no longer gates the label. No
threshold, model, dataset, or QA rule changed; Joanito had not been examined.

**Also applied**: final tables now carry `hgnc_symbol` (current HGNC, e.g.
CCN1/CCN2) plus `source_symbol` (CYR61/CTGF in GENCODE v26/Pelka). LY6A and
REG3B are recorded as having no direct human one-to-one orthologue (Ensembl
116: Ly6a → LY6S and Reg3b → REG1B, both `ortholog_one2many`, confidence 0) in
`config/nonhuman_candidate_orthology.tsv`, replacing the earlier "no human
gene" wording.

**Confirmed unchanged**: Pelka log2FC > 0 + FDR < 0.05; direction-only stem/TA
control; patient-level pseudobulk only; admission-QA markers never used to tune
thresholds; EMP1's pattern (normal-high in whole epithelium, higher than normal
stem/TA) is reported as-is.

**Gate before S1**: the Joanito label mapping must be shown to the user and
frozen before any Joanito DE. Interim labels: 6,296 `S2_pass_S1_pending`,
22,147 `S2_fail_S1_pending`.

## 2026-10-02 — Joanito downloaded; label map frozen (plan v3)

**Data**: Synapse terms accepted by the account owner; epithelial count matrix,
epithelial metadata, and clinical table downloaded with MD5 verification.
Terms: non-commercial use; no transfer or disclosure of data or derived
material. Because this repository is public, **no Joanito-derived numbers,
tables, figures, or labels are committed**; they are kept in the git-ignored
`step3_cancer/restricted/` record and under `DATA` (see below).

**Label map** (`step3_cancer/config/joanito_label_map.tsv`, rules only): Malignant =
iCMS2/iCMS3 cells in `Tumor`/`Tumor-2` samples (pooled per patient); Normal =
iCMS `Normal` cells in `Normal` samples; excluded = normal-like cells in tumour
samples, iCMS2/3-labelled cells in normal samples, and lymph-node samples.
Joanito has no normal stem/TA labels and none were constructed. The map was
reviewed and frozen before any Joanito pseudobulk or DE (plan v3).

## 2026-10-02 — Step 3 plan v3 executed (S1 Joanito)

**Review (v3)**: S1 is Joanito malignant vs normal epithelium only; the Pelka
tumour vs normal stem/TA comparison is an independent progenitor-specificity
check (P), not a substitute Joanito gate. A cohort-confounding check
(cohort × group patient table) and a sensitivity S1 restricted to cohorts with
both groups were added.

**Outcome (qualitative; numbers restricted)**: the full and sensitivity S1
models were highly concordant, so Joanito is retained; S1 passed the dataset
admission QA; `CRC_high` (S1 ∩ S2 passing P) and all other labels were
produced. Detailed counts, statistics, and gene-level observations are in the
restricted record.

**Open review items (no rule changed)**: whether to add an ambient
plasma-cell (immunoglobulin) flag; the progenitor check P removes few genes,
so separation of stem/WNT programmes from oncofetal genes is left to the fetal
intersection rather than to retuning P.

**Data-terms remediation**: commit `59e345a` had pushed Joanito-derived label
counts and per-patient cell counts. The branch was rebuilt from `c9916f8`
without any restricted file and force-pushed; restricted outputs are
git-ignored and mirrored to
`DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/`.

## 2026-10-02 — Cross-step 31-marker summary statistics; symbol-lookup corrections

**What**: `marker_summary/` tabulates log2FC, P value and FDR for the 31 Step 1
positive markers across 16 dataset contrasts (Step 2 primary, Step 2
benchmarks, Step 3), with a per-dataset measurement-status grid. Descriptive
only; no gate computed.

**Corrections found while building it**:
1. Step 2 audit (`Literature_31_gate_failure_audit.csv`, MAJOR_REVISION_LOG)
   reported CCN1/CCN2 as "not tested in HGCA". HGCA and Gao did test them under
   `CYR61`/`CTGF` (HGCA log2FC CCN1 0.73, CCN2 1.66). The other six "absent"
   genes are a mix: SPRR1A, SOX17, MMP7 and L1CAM are in the HGCA annotation but
   were removed by `filterByExpr`; LY6A and REG3B have no human one-to-one
   orthologue. The invalidated Step 2 outputs are left frozen; this note
   corrects the record.
2. Step 2 benchmark: Senger and Fawkner lookups now resolve CYR61/CTGF, and
   Fawkner absent genes are NA instead of 0 (see BENCHMARK_REPORT correction).
   Senger 14/29 (was 14/27); Fawkner 28/29 detected (was 26/31). No decision
   changes.

**Data terms**: the full summary contains Joanito-derived values and is git-ignored; a public version without Joanito columns or Step 3 labels is committed.

## 2026-10-02 — Step 2 human benchmark v2 (three new human contrasts)

**Plan**: `step2_benchmark_v2/docs/QUALIFICATION_PLAN.md`, frozen (`5c3bcef`)
before any new dataset was opened. Rule: TNFRSF12A fetal-high with P < 0.05 and
the 31-marker panel significantly fetal-biased (one-sided binomial P < 0.05).
No signature is built.

**Data added**: GSE158328 (Fawkner Visium), GSE158702 hashtag libraries
(Fawkner fetal scRNA), GSE185224 (Burclaff adult epithelium), GSE125970 (Wang
adult epithelium). DATA `link.md` files and `dataset.index.md` updated.

**Results**: only H-new2 (Fawkner fetal EPCAM+ scRNA, 22 hashtag samples, vs
Burclaff adult, 3 donors) qualifies among human contrasts: TNFRSF12A +1.83,
TACSTD2 +3.97, CLU +2.25 (all P < 0.001), ANXA1 −0.95 (n.s.); 18/24 markers
fetal-positive (binomial P = 0.011). H-new1 Visium fails and is
composition-confounded (fetal epithelial-rich spots carry ~4× less epithelial
UMI fraction; EPCAM adult-high). H-new3 passes panel bias (17/24) but not
TNFRSF12A. The Gao debug pair (Gao-original −0.08 → debug +0.70, n.s.) places
the Gao failure in the fetal arm (TNFRSF12A varies 34-fold across embryos,
lowest at 6–9 weeks), not the adult reference. Pikkupeura mouse cultures
qualify. Details: `step2_benchmark_v2/docs/BENCHMARK_V2_REPORT.md`.

**Technical fixes during the run** (documented in code): GEO mis-pairs the
Fawkner GEX/HTO libraries (re-paired by ≥ 99% barcode overlap); Burclaff's h5ad
keeps 23,170 genes, so counts for the same annotated cells come from the
full-gene per-donor matrices; the shared edgeR script gained an optional
minimum-units argument (H-new1 adult arm has two sections).

**Not done**: no replacement Step 2 signature; no change to Step 3; the
proposed manual TNFRSF12A follow-up (donor × gestational age × region ×
subtype) awaits review.

## 2026-10-02 — Review decisions on benchmark v2

1. **H-new2 = new human primary** developmental contrast, explicitly labelled a
   cross-study developmental contrast (Fawkner fetal vs Burclaff adult). It
   needs independent validation and is not final proof on its own.
2. **Visium dropped** from quantitative Step 2; kept only as a failed
   dataset/QC record.

Recorded in `step2_benchmark_v2/config/review_decisions.tsv` (surfaced as
`review_decision` in the qualification table), the v2 report, READMEs, the
dataset registry, and the 31-marker summary notes. Frozen verdicts unchanged.
Open: manual TNFRSF12A follow-up; choice of an independent validation source
for H-new2.
