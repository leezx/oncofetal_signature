# Conserved Intestinal Oncofetal Core (CIOC) — Methods

- **Version:** Core method v4.0, 2026-10-02. It supersedes v3.0, v2.0 and
  v1.0. v1.0 is retained as a sensitivity analysis (§7).
- **What v4.0 changes.** It is a logic correction, not a threshold change:
  **candidate nomination is separated from validation.**
  - The 31 literature-curated genes *define* the candidate universe.
  - Literature provenance becomes an annotation (§2) and is no longer a gate.
  - In v2.0–v3.0 a second literature filter ("A, or B with ≥ 2 studies") was
    applied to genes that were already literature candidates. That filtered
    the same source twice, with a threshold that was never part of the
    candidate definition.
- **Status: frozen.**
  - The v3.0 statement "permanently frozen" is superseded by this correction,
    which review approved as the last framework revision.
  - Disclosure: when v4.0 was approved, it was already known that RBP1 was
    the gene blocked only by the literature gate (it passed H, M and C under
    v3.0).
  - The correction applies identically to all 31 genes.
  - RBP1 was confirmed as an original member of the frozen candidate list
    (present since `3b08761`; the file is unchanged).
  - No further framework revisions will be made.
- **Earlier disclosure (v2.0–v3.0):** the per-gene human and mouse
  developmental values were visible when the rules were set. The rules were
  specified by review, not tuned on the outcome. TACSTD2 is a pre-declared
  sanity check (§8), not a target.
- **Inputs** are existing, version-controlled contrasts. No statistic is
  recomputed for gene selection.

## 0. Nomenclature and claim

- **Name:** *Conserved Intestinal Oncofetal Core (CIOC)*. It is not called
  "the intestinal oncofetal signature".
- **Operational definition:** we first assembled 31 candidate genes
  implicated in intestinal fetal, regenerative/revival or oncofetal biology
  from the literature. We then defined the CIOC by sequentially requiring:
  - replicated human developmental evidence;
  - in-vivo mouse fetal enrichment;
  - malignant epithelial reactivation in two independent CRC cohorts.
- **Claim scope:** the CIOC is a set of literature-nominated genes that show
  **developmental evidence** together with **malignant epithelial
  reacquisition**. It is **not** a set of universally fetal-specific
  markers. Statements about members use "developmental evidence" or
  "developmental enrichment", never "fetal-specific".
- **Disclosure:** membership is recorded in the restricted workbook (Joanito
  data-use terms, §10).

- **Wording for the members.** Write "N of the 31 literature-curated
  candidates satisfied all transcriptomic validation criteria." Do not write that all members have established literature
  evidence for intestinal fetal biology: RBP1's provenance is unresolved and
  stays recorded as such in the provenance table.
- **No ranking.** Members are not ranked or tiered in the manuscript. The
  evidence matrix shows every value as computed, including values at the
  threshold margin; thresholds are not adjusted in response.
- **Compartment.** CRC validation was performed specifically in epithelial /
  malignant epithelial pseudobulks. Downstream spatial or bulk readouts need
  cell-type resolution before a member-high signal (for example SPP1, which is
  also macrophage-expressed) is read as CIOC-high cancer cells.

## 1. Principle

**Literature nominates; data validate.** The CIOC is not an intersection of
single-dataset significance calls. Starting from the 31 literature
candidates, a Core gene must show:

- reproducible human fetal-associated expression (Gate H);
- in-vivo mouse fetal-high expression (Gate M);
- independently replicated CRC malignant-epithelial gain (Gate C).

**No single imperfect human developmental dataset holds a veto.** The rule is
fully deterministic: no gene is added or removed by judgement.

## 2. Candidate universe

- **Source:** the 31 Step 1 literature-curated candidates
  (`step2_fetal/config/literature_candidates_31.tsv`).
- **Provenance source:** the merged v1+v2 curated list
  (`step1_Literature_curated_intestinal_oncofetal_marker_candidates_merged_v1_v2.csv`).
- **Every one of the 31 is a literature candidate** ("Literature candidate =
  YES"). This means membership of the curated candidate universe. It does
  **not** mean equally strong published evidence that the gene is a bona fide
  fetal intestinal marker; the strength and nature of the evidence differ by
  gene and are given by the provenance annotation below.
- **Curation errors:** a gene found to have been added to the 31 in error
  would be removed from the candidate universe as a curation correction, not
  failed at a gate. None was found.
- **Symbols:** current HGNC, with legacy aliases resolved (CCN1 = CYR61,
  CCN2 = CTGF).
- **Orthology:**
  - LY6A and REG3B have no human one-to-one orthologue, so the human axes
    are not evaluable for them.
  - SPRR1A has no mouse orthologue, so the mouse axis is not evaluable for
    it.

## 3. Gate rules

Each gate returns **pass**, **fail** or **not evaluable (NE)**. NE means the
required dataset has no evaluable value for the gene (no one-to-one
orthologue, or `not_available_in_source`; after the data-QC amendment no gate
NA arises from expression filtering). **NE is never reported as a biological
fail.**

| Gate | Status | Datasets |
|---|---|---|
| (candidate definition) | not a gate | 31 literature candidates; provenance = annotation |
| H | mandatory, replicated | HGCA ≥ 9 PCW, H-new2 ≥ 9 PCW, Gao ≥ 9 W |
| M | mandatory | GSE230581 |
| C | mandatory, replicated | Joanito, Pelka |

Unless stated otherwise, **support** below means log2FC ≥ 0.5 and FDR < 0.05.

### Literature provenance (annotation; full-text audit; not a gate)

Evidence comes from a full-text audit of all 31 candidates
(`config/literature_audit_31.tsv`).

**Audit corpus (60 full texts):**
- the 55 Step 1 papers;
- Mustata 2013, Pikkupeura 2023 and Elmentaite 2021 (local PDFs);
- Fernandez-Vallone 2016, Karo-Atar 2022 and Vaquero-Siguero 2026 (Europe
  PMC open access).
- Not located: Fumagalli 2025, cited only for ANXA1.

**Search.** Every gene was searched in every text by human symbol, mouse
symbol and protein alias (e.g. TROP2, CX43, SCA1, CTGF, CYR61, FN14,
osteopontin), not only in the papers originally cited for it.

**Classification.** Each gene × study record was classified **from that
study's own data**:

| Class | Definition |
|---|---|
| **A** | fetal intestine: in-vivo fetal epithelium, or fetal-derived spheroids/organoids compared with adult |
| **B** | the gene is induced in, or marks, an injury-, infection- or YAP-driven regenerative fetal-like / revival epithelial state |
| **C** | CRC oncofetal / fetal-like tumour state |
| not counted | a mention (e.g. a gene used as a staining or proliferation marker), restating another paper, review statements, other tissues, or not found in the text |

**Annotation output** (`config/literature_provenance_31.tsv`, generated by
`scripts/literature_provenance.py`), per gene:
- the audited evidence classes (A/B/C);
- the supporting primary studies;
- a provenance status: *resolved*, or *unresolved* when no primary source was
  identified.

Unresolved provenance (RBP1, SPRR1A, EPS8L1) is reported as a curation gap,
not as biological evidence. Under v2.0–v3.0 these classes fed a literature
gate (20/31 passed after the audit); v4.0 uses them only as annotation.

**Data corrections from the audit:** AREG, EREG, EDN1 and IL1RN each gained a
second primary B study. RBP1 has no source in the corpus or in a Europe PMC
search.

### H — Human developmental, replicated (mandatory)

Three human fetal (≥ 9 weeks) vs adult epithelial contrasts:

| Contrast | Role |
|---|---|
| HGCA ≥ 9 PCW vs adult | primary: same atlas, donor pseudobulk edgeR |
| H-new2 ≥ 9 PCW: Fawkner vs Burclaff | cross-study, donor/sample pseudobulk edgeR |
| Gao fetal LI ≥ 9 W vs GSE103154 adult LI | effect-only, Welch |

- **Pass** if both hold:
  1. at least 2 of the 3 contrasts are measured and fetal-positive
     (log2FC > 0);
  2. at least 1 contrast shows support.
- **Wording for reporting:** human developmental evidence required
  concordant fetal enrichment across at least two independent comparisons,
  with statistical support in at least one dataset.
- **Role of H-new2:** H-new2 (Fawkner fetal × Burclaff adult) is
  cross-study, so stage is collinear with study/platform.
  - It serves as a **statistical-support dataset**, not a discovery dataset.
  - Direction is evaluated independently in the within-study (HGCA) and
    reference (Gao) datasets.
  - For every gene the matrix shows all three log2FC values and names the
    dataset(s) supplying statistical support. The H call is never shown
    alone.
- **Missing data:** NA neither lowers the denominator nor provides
  support.
  - With 2 evaluable contrasts, both must be fetal-positive (2/2).
  - With fewer than 2 evaluable contrasts, the call is **NE**.
- **Fail** otherwise.
- HGCA is the primary evidence but has no veto.

### M — Mouse in-vivo developmental (mandatory)

GSE230581 (Pikkupeura 2023 in vivo): E16.5 proximal SI epithelium vs adult
crypt epithelium, 3 vs 3, edgeR. The human gene is mapped through its Ensembl
one-to-one orthologue (CXADR: high-confidence one2many).

- **Pass** if the gene shows support.
- **NE** if there is no one-to-one orthologue or no evaluable value (NA).
- **Fail** otherwise.
- GSE44433 is replication evidence only (not selecting). It is reported
  beside every gene.
- **Decision record (v3.0).** A mouse replication requirement was
  evaluated. Candidate rules:
  - GSE230581 support plus GSE44433 fetal-positive where measured;
  - both positive with at least one supported;
  - GSE230581 support plus no significant GSE44433 contradiction.
- **GSE44433 coverage:** 23 of 30 mouse-applicable markers (older
  microarray; missing probes include CCN1, MIF, BASP1, CXADR, EPS8L1 and
  REG3B).
- **Per-gene consequences were shown before the decision.** The two
  direction-based rules would exclude TACSTD2 (GSE44433 −0.17, FDR 0.18) and
  SPP1 (−0.33, FDR 0.057) on non-significant differences. All three rules
  would exclude AREG, EREG and EMP1, which are significantly adult-high in
  GSE44433.
- **Review chose GSE230581 alone (unchanged from v2.0).** This choice was
  made with those consequences visible, and is disclosed as such.
  Significant GSE44433 contradictions are flagged in the evidence matrix.

### C — CRC malignant-epithelial, replicated (mandatory)

Two contrasts:

| Contrast | Model |
|---|---|
| Joanito 2022: malignant vs normal epithelium | patient pseudobulk, edgeR `~ cohort + group` |
| Pelka 2021: tumour vs normal epithelium | patient pseudobulk, edgeR |

- **Pass** if **both** contrasts show support.
- **Fail** if either measured contrast lacks support.
- **NE** if neither evaluable contrast fails but at least one has no evaluable value (NA).
- TCGA bulk is support only, with no veto.

### Panel-gene evaluation and NA provenance (data-QC amendment, frozen before recomputation)

**Problem.** A 31-gene feature audit (`scripts/feature_audit_31.py`) found
that almost every non-orthology NA in the gate contrasts came from the
genome-wide edgeR `filterByExpr` step, not from the gene being absent.
- **Present in the raw matrices but filtered:** for example, GJA1 in both CRC
  epithelial matrices (Pelka 3,017 counts; Joanito count in the restricted audit).
- **True zeros:** for example, mouse Sox17.
- **Symbol mapping:** no failures were found. Symbols, historical symbols
  and Ensembl IDs were all checked.

**Rule.** Genome-wide DE uses standard expression filtering
(`filterByExpr`), whose purpose is to limit the multiple-testing burden of
low-count genes. The 31 CIOC candidates were prespecified independently of
the expression data, so candidate-level validation retains any candidate
with nonzero counts and applies the same statistical model, with low-count
candidates explicitly flagged (item 4). The genome-wide filter therefore does
not decide whether a candidate can be evaluated.
1. **Filter exemption.** In every edgeR gate contrast, a panel gene is kept
   if it has any counts in the contrast's samples:
   `keep = filterByExpr OR (panel gene AND total count > 0)`.
   - The model, normalisation and dispersion estimation are otherwise
     identical to the original contrast.
   - **FDR universe.** The rescued candidates are fitted jointly with the
     genome-wide filtered genes in one edgeR fit, and Benjamini–Hochberg FDR
     is computed once over that enlarged tested set (genome-wide filtered
     genes + rescued candidates; 4–6 extra hypotheses per contrast). There is
     no separate 31-gene BH, so rescued and non-rescued candidates share one
     hypothesis family.
   - This applies to HGCA ≥ 9 PCW, H-new2 ≥ 9 PCW, GSE230581, Joanito and
     Pelka. Gao has no expression filter beyond all-zero.
2. **Zero expression.** A panel gene with zero counts in all relevant
   samples is called `zero_expression`.
   - It is evaluable and shows no enrichment: it counts as not positive in
     H and fails M or C.
   - An implausible zero in a source matrix (a ubiquitously expressed gene
     at zero in every unit) is flagged as a likely quantification artefact
     and treated as `not_available_in_source`.
3. **NA vocabulary.** NA reasons are restricted to:
   - `not_in_annotation`;
   - `no_1to1_orthologue`;
   - `zero_expression`;
   - `filtered_low_expression`;
   - `mapping_failure`;
   - `not_available_in_source`.

   The label "not measured" is no longer used. After the exemption, a remaining NA
   can only be one of these.
4. **Low-count display flag.** Panel genes whose mean CPM is below 1 in both
   arms of a contrast are flagged as low expression. The flag is shown in
   the matrix and is not selecting.

The gate rules (H, M, C thresholds) are unchanged. CIOC membership is
recomputed only after this amendment is applied to all 31 genes and all gate
contrasts.

**Implementation and verification.**
- **Scripts:** `scripts/panel_edger.R` and `scripts/panel_contrasts.py`.
  - They reproduce each original contrast exactly: HGCA ≥ 9 PCW and
    GSE230581 use `filterByExpr(design)`, `calcNormFactors` and robust QL;
    H-new2, Pelka and Joanito use `filterByExpr(group)`, `normLibSizes` and
    QL, with Joanito `~ cohort + group`.
  - The only change is the panel exemption.
- **Reproduction check** (genes matched by gene ID to the original DE
  tables):

  | Contrast | Tested genes (original → panel) | max \|ΔFDR\|, all genes | max \|ΔFDR\|, panel genes |
  |---|---|---|---|
  | HGCA ≥ 9 PCW | 17,804 → 17,808 | 1.1e-3 | 4.6e-4 |
  | H-new2 ≥ 9 PCW | 14,914 → 14,920 | 7.3e-4 | 2.9e-4 |
  | GSE230581 | 13,314 → 13,318 | 1.6e-3 | 8.6e-4 |
  | Pelka | 14,696 → 14,701 | 5.2e-4 | 7.0e-5 |

  - Joanito passes the same checks (values restricted).
  - log2FC reproduces to within 0.0005 in every contrast.
  - FDR in the rescued run equals BH over all tested genes (checked
    directly).
  - One genome-wide, non-panel gene crosses 0.05 (H-new2 SERPINB8,
    0.04998 → 0.05001). No panel-gene gate call changes.
- **Values:** `results/Panel_gate_values_public.csv`; Joanito values are
  restricted.
- **Remaining NA in the gate contrasts:**
  - `no_1to1_orthologue`: LY6A and REG3B in the human contrasts; SPRR1A in
    mouse.
  - `zero_expression`: SPRR1A in H-new2; Sox17 in GSE230581.
  - `not_available_in_source`: Gao MIF, zero in every Gao and GSE103154
    unit, which is implausible for a ubiquitous gene.
  - No `mapping_failure`.
- **Supportive datasets** (GSE44433, TCGA, cultures) were not re-run. Their
  NAs are labelled "supportive dataset; absent or filtered, not re-audited".

## 4. Core definition

**Conserved Intestinal Oncofetal Core (CIOC) = literature candidates passing Gates H, M and C.**
Every other gene is labelled "not Core" and lists, for each axis, whether it
failed or was not evaluable.

## 5. Evidence matrix (primary display)

The matrix has one row per literature candidate (31). Its columns are:

- Literature candidate (YES for all 31), evidence class annotation (or "provenance unresolved") and primary studies;
- HGCA ≥ 9 PCW, H-new2 ≥ 9 PCW, Gao ≥ 9 W, then the H call;
- mouse in vivo GSE230581 (M call) and mouse replication GSE44433;
- Joanito and Pelka, then the C call;
- TCGA (support);
- the deterministic Core label.

Each cell gives log2FC with significance, or NA.

### Discordance display (reported, never selecting)

Each gene carries a discordance flag listing:
- any human developmental contrast in the opposite direction;
- a significantly adult-high GSE44433 result;
- a significantly tumour-low TCGA bulk result.

Flags are displayed prominently (red) next to the calls, so heterogeneity is
visible rather than summarised away. Bulk discrepancies are described only as
*consistent with* compartment differences; no cause is claimed.

## 6. Supportive evidence (reported, never selecting)

- GSE44433 (mouse in-vivo replication).
- Pikkupeura 2021 cultures (laminin, collagen).
- HGCA all-fetal (stage-unresolved).
- Joanito sensitivity (cohorts with both groups).
- Pelka progenitor check P.
- TCGA T1 / paired / vs GTEx.

## 7. Sensitivity analysis using single-dataset hard intersections (former method v1.0)

- **Rule:** literature candidates ∩ G1 HGCA ≥ 9 PCW ∩ G2 GSE230581 ∩ G3 Joanito ∩ G4
  Pelka.
- **Pass in each gate:** log2FC ≥ 0.5 and FDR < 0.05.
- **NA** (no evaluable value) means not evaluable, never "fail".
- **Developmental intersection (G1 ∧ G2):** GJA1, CLU, ANXA6, SPP1, RBP1.
- **Purpose:** the hard-intersection result is a **sensitivity analysis**,
  not an alternative Core. It shows that the conclusions do not depend on the
  evidence-integration framework. Its membership is in the restricted
  workbook.
- **GJA1** is evaluated in CRC under the panel-gene exemption. In the public
  Pelka contrast it is +3.89 (FDR < 1e-4, low count). Its Joanito status is
  in the restricted workbook.

## 8. Pre-declared sanity checks

- **TACSTD2 (TROP2):** if v2.0 excludes it, the axis that excludes it and the
  cross-dataset evidence for that exclusion are reported explicitly.
- **TNFRSF12A:** a literature candidate (audited classes A and B). Its status is reported whatever
  it is. Statement if excluded: *TNFRSF12A is mechanistically linked to the
  oncofetal program but is not itself a member of the stringent conserved
  oncofetal Core.*
- **Neither gene's status may motivate a rule change.**

## 9. Known design caveats (disclosed)

- **HGCA:** fetal (HDBR) and adult (transplant-donor) cohorts were collected
  separately within one atlas framework.
- **H-new2:** stage is collinear with study/platform.
- **Gao:** cross-platform, with adult n = 2.
- **GSE230581:** fetal whole epithelium vs adult crypt-only epithelium.
- **Literature provenance:** this comes from the full-text audit.
  - Fumagalli 2025 could not be located.
  - Gene-specific statements may also exist only in figures or supplements
    that were not text-searchable.
  - RBP1, SPRR1A and EPS8L1 have no identified primary source (provenance
    unresolved). They remain candidates because they are in the frozen
    candidate list. If a later check shows a curation error, the gene is
    removed from the universe, not failed at a gate.

## 10. Outputs and data-use restriction

- **Restricted (Joanito terms forbid disclosing derived material):** the
  Joanito columns, the C call and the Core labels are written only to the
  full workbook (git-ignored, mirrored to
  `DATA/.../restricted_joanito/core_oncofetal/`).
- **Public workbook** keeps the literature annotation, H, M, Pelka, TCGA and
  the supportive evidence.
- **Gate funnel workbook** (`results/CIOC_gate_funnel_31.xlsx`, restricted;
  `scripts/build_gate_funnel.py`): one row per candidate, with Literature
  candidate = YES, then Gates H, M and C, the reasons for filtering, and a
  final Core column (YES or blank).

## 11. Change control

v4.0 is frozen and is the last framework revision. Any future change to a
gate, dataset, threshold or rule defines a different gene set with a new
name; it is not a revision of the CIOC.

## 12. Change log

| Version | Change |
|---|---|
| v1.0 | Four-dataset intersection (now the stringent sensitivity set). |
| v2.0 | Evidence framework (L, H replicated, M, C replicated); NE introduced. |
| v3.0 (final) | L evidence from the full-text audit (rule unchanged); explicit H missing-data rule; mouse replication evaluated and not adopted (GSE44433 supportive, decision disclosed). No further framework revisions. |
| v3.0, presentation only (after freeze) | CIOC name and claim scope (§0); H reporting wording and H-new2 role; per-gene H support source; discordance flags; v1.0 renamed "sensitivity analysis using single-dataset hard intersections". **Membership unchanged.** |
| v4.0 | Logic correction: literature defines the candidate universe (all 31 = YES); provenance (A/B/C, studies, resolved/unresolved) is annotation, not a gate. Core = H ∧ M ∧ C. Applied identically to all 31 genes. Last framework revision. |
| v4.0 data-QC amendment | Panel genes exempt from the genome-wide expression filter; explicit NA provenance vocabulary; zero expression treated as evaluable; low-count flag. Gate rules unchanged. Frozen before recomputation. |
| v4.0 final freeze (wording only) | Rationale and FDR universe of the panel exemption stated (one joint fit, BH over genome-wide filtered + rescued genes; FDR reproduction table); NE wording no longer says "not measured"; §0 member wording, no-ranking and compartment statements. No numerical or membership change. **CIOC membership frozen** (membership in the restricted workbook, §10). |
