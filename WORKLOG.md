# Worklog — intestinal oncofetal signature

Chronological append-only record; newest entries are at the bottom.

## Progress tracker

| Step | Status | PR | Notes |
|---|---|---|---|
| 1. Literature candidate collection | completed | — | Curated evidence tables under `docs/step1_literature_candidates/` |
| 2. Conserved developmental axis | **invalidated; dataset review required** | [#1](https://github.com/leezx/oncofetal_signature/pull/1) | 706-gene output is provisional and must not be used downstream |
| 2b. Human benchmark v2 | **closed (review round 3); stage-resolved addendum v1.2 run** | [#3](https://github.com/leezx/oncofetal_signature/pull/3) | H-new2 = tier A (cross-study, not used alone for discovery); bulk H-bulk1/2 = 1 dataset, fail; Visium dropped; no more datasets |
| 2c. Fetal epithelial state (manual) | **executed; in review** | (this PR) | Gao early-embryo effect; TNFRSF12A pan-fetal-epithelial, stem/TA-enriched detection; no adult arm |
| Core. Conserved Intestinal Oncofetal Core | **CIOC method v4.0 frozen (last revision)** | [#6](https://github.com/leezx/oncofetal_signature/pull/6) | Literature31 (candidate definition) → H(replicated) ∧ M ∧ C(replicated); v1.0 intersection = stringent sensitivity set; membership restricted (Joanito) |
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

## 2026-10-02 — Benchmark v2 addendum: primary-tissue bulk (H-bulk1, H-bulk2)

**Plan**: addendum v1.1 frozen (`cd39210`) before data were opened.
**Source tracing**: Finkbeiner 2015 and Senger 2018 "primary tissue" both use 6
Roadmap fetal SI total-RNA samples (SRP001371, six donors, 91–120 days) and HPA
E-MTAB-1733 adult SI/duodenum (ERP003613). H-bulk2 (duodenum only) is a subset
of H-bulk1. Senger Table 2 lists GSM1059508 for two ages (paper inconsistency).
**Data**: recount3 raw counts for both studies (new DATA dataset
`bulkRNAseq/recount3_Roadmap_HPA_intestine`); HPA technical runs summed per
sample (6 adult samples).
**Results**: H-bulk1 TNFRSF12A −0.50 (P = 0.20), 13/26 markers fetal-positive;
H-bulk2 −0.90 (P = 0.12), 16/27; both fail. OLFM4 sanity control strongly
adult-high (−6.9). TACSTD2 +2.49 and RBP1 +1.50 fetal-high; ANXA1 −2.66
adult-high; CLU flat.
**Tiers recorded**: A = H-new2 (qualifies); B = H-bulk1/2 (fail); C = HGCA, Gao,
Senger (fail); Visium dropped.
**Stop condition reached**: bulk contrasts also give TNFRSF12A fetal ≤ adult;
the manual single-cell TNFRSF12A examination is the proposed next step (not
started).

## 2026-10-02 — Benchmark v2 review round 3: bulk benchmark closed

**Decisions**: PR #3 kept; bulk benchmark closed; no further human fetal/adult
datasets. H-bulk1/H-bulk2 count as 1 independent dataset, 2 related contrasts
(`independent_source` column added to `review_decisions.tsv`). Bulk failure
worded as "maturation signal recovered, TNFRSF12A not reproduced, cause not
established". Human evidence read as compartment-dependent. H-new2 is not to be
used alone for genome-wide discovery (dataset-selection concern).
**Next**: manual developmental-state analysis of Gao + Fawkner
(`step2_fetal_state/`).

## 2026-10-02 — Fetal epithelial state analysis (Gao + Fawkner), `step2_fetal_state/`

**Plan**: frozen (`4747efe`) before expression was read; addenda v1.1 (hashtag
QC: per-cell sex discordance; ABF1/AAU2 excluded, ABZ2 flagged) and v1.2
(state-label confidence; Gao states exploratory) committed before test genes
were tabulated. Fawkner hashtags mapped to samples via the authors' Mendeley
supplement (new file in DATA raw).
**Results**: Gao TNFRSF12A instability is driven by ≤ 8 W embryos (median
log2 2.0 vs 4.5 at ≥ 9 W; 91-fold range overall, 8-fold within ≥ 9 W, no age
trend after 9 W). Fawkner TNFRSF12A is flat across 8–20 PCW while LGR5/ASCL2/
MKI67/ANKRD1 track age. TNFRSF12A is expressed across all fetal epithelial
states (highest detection Stem 61%, TA 50%); no discrete TNFRSF12A-high
progenitor sub-state; variance mostly within state. ANXA1 marks the Fawkner
Fetal-progenitor state. No consistent region effect in Fawkner; Gao LI > SI.
**Open**: state-matched fetal vs adult comparison, or Gao ≥ 9 W vs adult (new
frozen contrast either way); final human discovery contrast still undecided.

## 2026-10-02 — Benchmark v2 addendum v1.2: stage-resolved re-analysis

**Plan**: frozen (`618356e`) before computation; breakpoint 9 weeks from the
fetal-only state analysis; early < 9 W, mid/late ≥ 9 W, reported ages.
**Scope**: HGCA (6 early / 10 mid-late / 7 adult donors), Gao-original and
H-new3 (Gao 5–6 early / 7–9 mid-late embryos); Senger not splittable (all
fetal ≥ 11 W GA); Roadmap bulk already all ≥ 9 W.
**Results**: early vs adult TNFRSF12A negative everywhere (HGCA −1.58,
P = 0.006); mid/late vs adult ≈ 0 (HGCA −0.01, Gao-original −0.01) or +0.69
(H-new3, P = 0.089, panel 20/24); mid/late vs early positive (HGCA +1.60,
P = 0.009 — independent support for the breakpoint). No contrast qualifies;
promotion criterion not met; panel fraction rises only in H-new3 (71 → 83%).
**Updated**: qualification + matrix + heatmap (21 contrasts), stage DE tables,
review_decisions, marker summary (public committed; full mirrored to DATA
restricted), dataset registry notes.

## 2026-10-02 — Benchmark v2 addendum v1.3: stage-resolved H-new2; results audit

**Audit**: every contrast run so far is present in `marker_summary/results`
(public: all except Joanito; full: restricted mirror in DATA). Gap found:
H-new2 (Fawkner) had no stage split. Frozen as v1.3 (`9571e29`), same
breakpoint and rules.
**Results**: early (2 units, 8 PCW) vs adult TNFRSF12A +0.86 (P = 0.015);
mid/late vs adult +1.94; mid/late vs early +1.08 (P = 0.008). Within-fetal
rise across 9 weeks now seen in HGCA, Fawkner (significant) and Gao (P = 0.07).
**Added**: `marker_summary/results/Fetal_stage_decomposition_summary.csv`
(all human fetal-vs-adult parents × stage resolution); marker summary now 36
contrasts; matrix/qualification/heatmap 24 contrasts.

## 2026-10-02 — Benchmark v2 addendum v1.4: Core-gate dataset selection benchmark

**Plan**: frozen (`da9a1cb`): design eligibility E1–E5 (same study, in vivo,
epithelium-resolved, ≥ 3 replicates/arm, replicate-level FDR) first, then
31-marker recovery (% positive; + FDR < 0.05; + log2FC ≥ 0.5). TNFRSF12A not
an admission criterion. No gate selected.
**Correction**: Pikkupeura pure in-vivo contrast = GSE230581 (Step 2 M1), not a
second dataset; independent mouse in-vivo replicate = GSE44433.
**Results**: eligible — HGCA ≥ 9 PCW 15/25 positive, 7 FDR-supported (7
adult-high); Pikkupeura in vivo 18/25, 15 FDR-supported; GSE44433 10/23, 8.
Ineligible H-new2 ≥ 9 PCW 16/22, 14 FDR-supported. HGCA ∩ Pikkupeura in vivo
FDR-supported fetal-high: GJA1, CLU, ANXA6, SPP1, RBP1.
**Outputs**: `Core_gate_dataset_recovery.csv`, `Core_gate_marker_by_dataset.csv`
(step2_benchmark_v2/results/tables and marker_summary/results).

## 2026-10-02 — Core oncofetal method v1.0 and gate workbook (`core_oncofetal/`)

**Method**: frozen (`9e150dd`) before CRC-gate values of the developmental
genes were read. Gates: G1 HGCA ≥ 9 PCW fetal vs adult epithelium; G2
GSE230581 E16.5 vs adult crypt (mouse in vivo); G3 Joanito malignant vs normal;
G4 Pelka tumour vs normal; pass = measured, log2FC ≥ 0.5, FDR < 0.05.
Supportive evidence (H-new2, Gao ≥ 9 W, GSE44433, Pikkupeura cultures, TCGA,
Joanito sensitivity, Pelka P) reported only.
**Results**: Developmental Core (G1 ∧ G2) = GJA1, CLU, ANXA6, SPP1, RBP1.
TNFRSF12A fails G1 and G2. Final Core membership (needs G3) is recorded only
in the restricted workbook (DATA restricted_joanito/core_oncofetal/).
**Outputs**: public `core_oncofetal/results/Core_oncofetal_gate_statistics_public.{xlsx,csv}`.

## 2026-10-02 — Core method v2.0: literature-anchored evidence framework

**Why**: review — a four-dataset intersection makes one imperfect human dataset
(HGCA, strict literature recovery 7/25) an absolute veto, and treated "not
measured" as fail.
**Method**: frozen (`ef58297`) before labels were computed. Axes L (curator code
contains A, or B with ≥ 2 named studies), H (HGCA/H-new2/Gao ≥ 9 W: ≥ 2/3
fetal-positive and ≥ 1 with log2FC ≥ 0.5 & FDR < 0.05), M (GSE230581), C (Joanito
and Pelka). Calls pass/fail/not evaluable. v1.0 retained as stringent
intersection sensitivity set.
**Public results**: L passes 16/31 (RBP1 fails: its curated provenance has no
named study). H passes TACSTD2 (via H-new2/Gao), not SPP1 (Gao −0.01, H-new2 not
measured). TNFRSF12A fails H and M. Final Core labels are restricted (Joanito).
**Outputs**: 31-gene evidence matrix workbook (public / restricted).

## 2026-10-02 — Core method v3.0 (final framework): literature audit, missing-data rule

**Literature audit**: all 31 genes searched (symbols, mouse symbols, aliases)
in 60 full texts (55 Step 1 papers + Mustata 2013, Pikkupeura 2023, Elmentaite
2021, Fernandez-Vallone 2016, Karo-Atar 2022, Vaquero-Siguero 2026); A/B/C
classified from each study's own data (`config/literature_audit_31.tsv`).
L passes 20/31 (v2.0: 16): AREG, EREG, EDN1, IL1RN gained a second primary B
study. RBP1 has no source in the corpus or a Europe PMC search (fails L as
unsourced). Fumagalli 2025 not located.
**Mouse axis**: GSE44433 replication rules evaluated (coverage 23/30); both
direction-based rules would exclude genes on n.s. differences (identities
restricted; wording corrected 2026-10-05); review
kept GSE230581 alone (decision disclosed in Methods); GSE44433 contradiction
flag added to the matrix.
**H missing data**: NA neither lowers the denominator nor supports; 2
evaluable → 2/2 required; < 2 → not evaluable.
**Frozen**: `05f8812` (before labels). Final Core labels are restricted
(Joanito); public workbook updated.

## 2026-10-02 — CIOC v3.0 permanently frozen; presentation finalised

**Decision (review)**: v3.0 frozen permanently; no further membership
optimisation. Name: Conserved Intestinal Oncofetal Core (CIOC); claim =
developmental evidence + malignant epithelial reacquisition, not fetal-specific
markers.
**Presentation-only changes (membership unchanged)**: H wording ("concordant
across ≥ 2 independent comparisons, statistical support in ≥ 1"); H-new2 = statistical-
support dataset; matrix shows all three human log2FC plus the support source;
discordance flags (opposite human contrast, significant GSE44433 adult-high,
significant TCGA tumour-low) in red; v1.0 renamed "sensitivity analysis using
single-dataset hard intersections". Methods §0/§11/§12 updated; the public
GJA1/Joanito statement removed. Restricted manuscript text (definition, Methods,
Results, legends, reviewer answers) in `results/CIOC_manuscript_text.md`
(git-ignored, mirrored to DATA restricted).
**Next (not started)**: program coherence of the CIOC in independent scRNA/spatial
data; relation to revCSC/proCSC, TWEAKR and YAP.

## 2026-10-02 — CIOC method v4.0: literature nominates, data validate

**Why (review)**: the v2.0–v3.0 literature gate re-filtered the same literature
source that defined the 31 candidates, with a threshold ("A, or B with ≥ 2
studies") that was never part of the candidate definition.
**Change**: all 31 = Literature candidate YES; A/B/C classes, primary studies
and resolved/unresolved status are annotation; Core = H ∧ M ∧ C. Applied to all 31.
The affected gene was confirmed in the original frozen candidate list
(`3b08761`, unchanged). **Disclosed**: which gene the change would affect was
known before v4.0 was approved (identity restricted). [Wording corrected
2026-10-05 to remove a restricted membership inference.]
**Outputs**: Methods v4.0; provenance as annotation; public evidence matrix;
restricted gate funnel workbook (`build_gate_funnel.py`) and manuscript text.
Final membership restricted (Joanito). Open curation item: unresolved
nomination sources (see provenance table).

## 2026-10-02 — Data-QC fix: panel genes exempt from genome-wide expression filter

**Bug (review)**: "not measured" mixed orthology gaps, expression-filter removal,
mapping and zero expression. **Audit** (`feature_audit_31.py`): no mapping failures;
every non-orthology NA in the gate contrasts came from edgeR `filterByExpr` (e.g.
GJA1 present in both CRC matrices) or true zero (mouse Sox17, SPRR1A in H-new2);
Gao MIF zero in every unit (implausible; not_available_in_source).
**Fix** (frozen `ae871ec` before recomputation): 31 panel genes exempt from the
filter in all edgeR gate contrasts (any counts kept); zero_expression evaluable
(no enrichment); fixed NA vocabulary; low-count display flag. Original contrasts
reproduced (|Δlog2FC| ≤ 0.0005). Gate rules unchanged.
**Effect (public)**: SPP1 now passes H (H-new2 +7.30); MSLN passes M (+7.18);
GJA1 evaluable in Pelka (+3.89); LY6A and IL1RN now measured mouse fails. Final
membership restricted (Joanito).

### 2026-10-02 — CIOC final freeze (wording only)
**FDR universe confirmed**: rescued panel genes are fitted jointly with the
genome-wide filtered genes; BH is computed once over that enlarged tested set
(4–6 extra hypotheses per contrast), not a separate 31-gene BH. Max |ΔFDR| vs
original ≤ 3.7e-3 (all genes), ≤ 8.6e-4 (panel genes); no panel gate call
changes. Methods now state the exemption rationale and FDR universe, member
wording, no-ranking and CRC-compartment statements; NE wording no longer says
"not measured". No numerical or membership change; membership frozen
(restricted).

### 2026-10-02 — CIOC-extended (exploratory)
Frozen gates H, M, C applied genome-wide without the Literature-31
restriction (`core_oncofetal/scripts/build_extended_core.py`). Universe and
counts are in the restricted funnel workbook; the CIOC is fully recovered
inside the extended set. Outputs restricted (Joanito), git-ignored, mirrored
to DATA.

### 2026-10-02 — Level 2 renaming; Extended CIOC plan
Review: the genome-wide H∧M∧C set includes strong, concordant non-epithelial
programs (collagen, smooth-muscle, endothelial, lymphoid genes), so it is
not an epithelial oncofetal signature. Renamed: H∧M∧C = cross-species
fetal–CRC candidates; H∧C = human conserved fetal–CRC candidates (Level 2,
discovery universe, not for scoring). Script renamed to
`build_genomewide_candidates.py`; old CIOC_extended outputs removed.
Internal validation: genome-wide search recovers all CIOC genes and no other
Literature-31 gene. Extended CIOC (Level 3) = Level 2 ∩ Gate E (epithelial
compatibility) ∩ CIOC coherence, rules to be frozen before computation
(`core_oncofetal/docs/EXTENDED_CIOC_PLAN.md`). Open decision: independent
all-compartment CRC atlas (Tabula Sapiens LI is adult-normal only).

### 2026-10-05 — Extended CIOC (construction complete)
Rules v1.0 frozen before computation (`4f2d9bd`):
- **Gate E:** E1 fails if the log2 ratio is < −3 in both atlases or < −5 in
  either; E2 requires detection in ≥ 5% of tumour-derived epithelial cells
  in ≥ 1 atlas.
- **Coherence:** within-stratum metacells, leave-one-out CIOC score,
  cell-cycle adjusted, against 1,000 expression-matched null sets; pass if
  positive in both atlases and empirical P < 0.05 in ≥ 1.

Applied once: 338 → 189 (Gate E) → 38 Extended CIOC, which includes 2 of
the 8 CIOC genes. Four CIOC genes fail Gate E because they are dominated by
non-epithelial compartments in tumour tissue; two fail coherence.
Signature construction ends here. Outputs are restricted (git-ignored,
mirrored to DATA).

### 2026-10-05 — ECOS-38 rename; Extended CIOC v2.0 (epithelial-only)
**Rule freeze.** The rule was frozen before computation (`72a98b6`). The
plan separates state specificity from lineage specificity.

**Rename (membership unchanged).** The v1.0 output "Extended CIOC" (38
genes) is now **ECOS-38**, the Epithelial-Compatible Oncofetal Signature
for mixed-cell data.
- Script renamed: `ext_04_extended_cioc.py` → `ext_04_ecos38.py`.
- Superseded files are kept unmodified in DATA
  `restricted_joanito/core_oncofetal/superseded_2026-10-05_v1.0_named_Extended_CIOC/`:

  | File | md5 |
  |---|---|
  | `Extended_CIOC_calls.csv` | `553e042ba848f12de3b8b88ced84a479` |
  | `Extended_CIOC.gmt` | `8a0a56c552239bf0c4dc94fea39b329e` |
  | `Extended_CIOC.xlsx` | `5ebf7a9f62c961f9a3ce2fbff9d7e196` |

- The ECOS_38 files were produced from these by relabelling only. The 38
  genes were asserted identical.

**Extended CIOC v2.0** (`ext_05_extended_cioc.py`, applied once).
- **Rule:** epithelial detectability (≥ 5% in ≥ 1 atlas) and T > 0 in both
  atlases. Gate E E1 is an annotation only.
- **Result:** 338 → 222 detectable → 222 Extended CIOC. As recorded before
  computation, coherence positivity excluded no gene.
- **Annotations:** 56 genes are high-confidence CIOC-coherent (P < 0.05 in
  ≥ 1 atlas); 10 reach P < 0.05 in both atlases.
- **Overlap:** all 8 CIOC genes and all 38 ECOS-38 genes are included.
- **Reproduction:** the 189 genes tested in v1.0 reproduce exactly.

Signature construction is complete.

### 2026-10-05 — Methods.md check; generated references
`core_oncofetal/Methods.md` was checked against the scripts, unit tables and
provenance files and then corrected. The file is restricted (it contains
memberships): it is git-ignored and mirrored to DATA. The original is backed
up as `Methods_original_2026-10-05.md` (md5
`05122d57eee6ee3aa4abe550cc2f7d80`).

**Corrections:**
- **Gao comparison:** now stated explicitly as Gao fetal LI ≥ 9 W (7
  embryos) vs GSE103154 adult LI (P1, P2), cross-platform, effect-only, with
  a descriptive Welch P.
- **Unit counts:** given for every gate dataset.
- **H-new2:** hashtag-sample pseudoreplication and the two retained
  QC-failed units are disclosed.
- **9-week breakpoint:** its origin is stated (Gao fetal-only TNFRSF12A).
- **Statistical support:** disclosure of which dataset supplied it (H-new2
  alone for 4 CIOC genes and 48 of the 338; Gao alone for 2 of the 338).
- **Orthology exceptions:** CXADR, LY6A, REG3B and SPRR1A.
- **FDR universe** for the prespecified panel.
- **Unresolved provenance:** RBP1, SPRR1A and EPS8L1.
- **Extended CIOC:** stated that coherence excluded no gene (membership =
  H/M/C + detectability).
- **Gate E:** compartments now include mast cells, and the timing of the −5
  arm is stated.
- **Che:** chemotherapy disclosed.
- **Revision history:** disclosed.

**References.** `core_oncofetal/scripts/build_references.py` generates
`results/Methods_references.{md,tsv}`. Sources:
- `literature_audit_31.tsv` (A/B/C primary studies);
- `paper_inventory.tsv` DOIs, plus `config/reference_doi_supplement.tsv`
  for 6 DOIs resolved by Crossref search with titles verified;
- `DATASET_REGISTRY.csv`;
- `config/reference_methods.tsv`.

Metadata comes from Crossref, cached in `config/crossref_cache.json`. The
dataset registry gained the Khaliq, Che and Tabula Sapiens rows, and the
final gate roles were added for HGCA, Gao, GSE230581 and GSE44433.

### 2026-10-05 — Supplementary tissue / cancer breadth annotation
The plan was frozen before computation (`9c958ee`) and applied once by
`core_oncofetal/scripts/breadth_annotation.py`.

**Inputs:**
- Cao 2020 fetal atlas (GSE156793 aggregated tables): 8 organ epithelia.
- UCSC Toil TCGA + GTEx: 21 carcinomas vs matched normal.

**Results.** All 338 candidates were annotated, with no membership change.
- Pan-tissue pan-cancer oncofetal genes: 0 of the CIOC, 2 of ECOS-38, 10
  of the Extended CIOC and 12 of the 338.
- Per-gene classes are in the restricted workbook.

**Calibrators.** EPCAM, CDX2 and MKI67 behaved as expected. COL1A2 leaks
into the Cao epithelial pseudobulks and is caught only by the compartment
flag, which is disclosed.

**Records.** New DATA datasets were registered with `link.md` files. Cao
and Toil were added to the dataset registry and the references (14 dataset
references). Per-gene outputs are restricted.

### 2026-10-05 — Restricted-membership wording scrub (public docs)
A leak check found that public docs named genes in ways that implied CIOC,
338 or Extended CIOC membership, or stated Joanito-derived results.

Affected files:
- `docs/HANDOFF.md`;
- `core_oncofetal/README.md`;
- `EXTENDED_CIOC_PLAN.md` and `BREADTH_ANNOTATION_PLAN.md`;
- the RBP1 disclosure and GSE44433-consequence passages in
  `Method_core_oncofetal_validated.md`;
- WORKLOG lines.

The names were replaced with restricted pointers. The text is clean from
this commit onward. **Earlier pushed commits still contain these mentions in
history.** Removing them requires a history rewrite and force-push of
`analysis/core-oncofetal`; that is the user's decision and has not been done.
