# Session handoff: oncofetal signature construction

Last updated: 2026-10-05.
- Branch: `analysis/core-oncofetal`, last commit `bc07628`, pushed.
- PR: #6, chained on PRs #2 → #3 → #4 → #5, all open for review.

**Read this first, then `core_oncofetal/HANDOFF_RESTRICTED.md`.** That file
is local only and git-ignored; it holds memberships and key values. A copy
is in DATA `restricted_joanito/core_oncofetal/`.

## 1. Where the project stands

**Signature construction is complete and frozen.** Do not change gates,
thresholds or membership unless a documented review decision explicitly
asks for it, and then follow the freeze protocol in §5.

There are four fixed objects:

| Object | Size | Definition | Intended use | Built by |
|---|---|---|---|---|
| **CIOC Core** (Conserved Intestinal Oncofetal Core) | 8 | Literature-31 candidates passing Gates H ∧ M ∧ C (method v4.0 + data-QC amendment) | Stringent biological anchor | `core_oncofetal/scripts/build_core_gates.py`, `build_gate_funnel.py` |
| **Cross-species fetal–CRC candidates** | 338 | The same frozen H/M/C gates applied genome-wide | Discovery universe, **not a scoring signature** | `build_genomewide_candidates.py` |
| **Extended CIOC** (v2.0, epithelial-only) | 222 | 338 ∩ epithelial detectability (≥ 5% of tumour-derived epithelial cells in ≥ 1 atlas) ∩ positive CIOC coherence in both Khaliq and Che | Epithelial-resolved (cell-level) data | `ext_05_extended_cioc.py` |
| **ECOS-38** (Epithelial-Compatible Oncofetal Signature) | 38 | 338 ∩ Gate E (compartment compatibility) ∩ stringent coherence (positive in both atlases, empirical P < 0.05 in ≥ 1) | Bulk RNA-seq and mixed-cell or unresolved spatial data | `ext_04_ecos38.py` |

Also kept, never used for scoring:
- "Human conserved fetal–CRC candidates" (H ∧ C, 1,111 genes);
- "Core human" (H ∧ C among the 31).

**Agreed facts to state honestly:**
- **Extended CIOC coherence excluded no gene.** All 222 detectable genes
  are positive in both atlases, so membership = H/M/C + detectability. The
  56 "high-confidence CIOC-coherent" genes (empirical P < 0.05 in ≥ 1
  atlas) are an annotation only.
- **ECOS-38 sits inside the Extended CIOC numerically.** ECOS-38 equals the
  high-confidence genes that pass Gate E, and all 38 are Extended CIOC
  members. Conceptually, present it as a separate branch from the 338.
- **Internal validation.** The genome-wide search recovers all 8 CIOC genes,
  and no other Literature-31 gene passes.
- **State vs lineage specificity.** SPP1, CCN2, CLU and RBP1 are induced in
  epithelium but are dominated by myeloid, stromal, endothelial or mast
  compartments in tumour tissue.
  - They fail Gate E, so they are absent from ECOS-38.
  - They remain in the CIOC and the Extended CIOC.
  - An 8-gene CIOC score in bulk data is microenvironment-dominated.
- **Terminology:**
  - Tumour-tissue epithelial cells in Khaliq and Che are "tumour-derived
    epithelial cells". They are never called malignant (no CNV inference).
  - ECOS-38 is "epithelial-compatible", never "epithelial-specific".
  - Khaliq and Che are "independent replication within the refinement
    stage", not external validation.

## 2. Gate definitions (frozen)

The full text is in `core_oncofetal/docs/Method_core_oncofetal_validated.md`
(CIOC) and `core_oncofetal/docs/EXTENDED_CIOC_PLAN.md` (Level 2, Gate E,
coherence, ECOS-38, Extended v2.0).

"Support" = log2FC ≥ 0.5 and FDR < 0.05.

**Gate H: human developmental** (fetal ≥ 9 weeks vs adult)
- **Datasets:**
  - **HGCA:** 10 fetal donors (9.2–17 PCW) vs 7 adults; edgeR QL robust,
    donor pseudobulk.
  - **H-new2:** Fawkner-Corbett fetal (20 hashtag samples from 10 fetuses,
    9–20 PCW) vs Burclaff adult (3 donors); edgeR QL.
  - **Gao-original late:** Gao fetal LI, 7 embryos at 9–25 W (GSE95630),
    vs GSE103154 adult LI (P1, P2). It is cross-platform and effect-only,
    with a median-based log2FC and a **descriptive** Welch P.
- **Pass:** ≥ 2 of the 3 contrasts fetal-positive, and ≥ 1 supported.
  - NA neither lowers the denominator nor supports.
  - With only 2 evaluable contrasts, both must be fetal-positive.
  - With < 2 evaluable contrasts, the call is NE.
- **The 9-week breakpoint** comes from fetal-only Gao TNFRSF12A.

**Gate M: mouse in vivo**
- GSE230581: E16.5 vs adult proximal small-intestinal crypt epithelium,
  3 vs 3, edgeR QL robust.
- Mapping uses Ensembl 116 one-to-one orthologues, with these exceptions:
  - CXADR uses a one-to-many orthologue;
  - LY6A and REG3B are mouse-only markers;
  - SPRR1A has no mouse orthologue.
- GSE44433 is a supportive flag only.

**Gate C: CRC epithelial**
- Joanito (61 malignant vs 24 normal patient pseudobulks, `~ cohort + group`)
  **and** Pelka (62 tumour-specimen vs 35 normal-specimen epithelium) must
  both show support.
- TCGA is supportive only.

**Data-QC amendment**
- The 31 prespecified genes are exempt from `filterByExpr`.
- They are fitted jointly, and BH is computed over the enlarged tested set.
- Zero expression is evaluable and counts as no enrichment.
- NA reasons use a fixed vocabulary.
- A low-count flag is shown when mean CPM < 1 in both arms.

**Gate E (ECOS-38 only)**
- **Ratio:** log2((tumour-derived epithelial CPM + 1) / (top non-epithelial
  compartment CPM + 1)), using patient-median pseudobulks.
- **E1 fails** if the ratio is < −3 in both atlases, or < −5 in either.
- **E2 detectability:** ≥ 5% of tumour-derived epithelial cells in ≥ 1
  atlas.

**Coherence**
- **Units:** within each sample (≥ 200 epithelial cells), k-means metacells
  of about 20 cells, with ≥ 10 metacells.
- **Score and statistic:**
  - The CIOC score leaves the gene out when the gene is itself a CIOC
    member.
  - Values are residualised on S/G2M.
  - Spearman correlation, then T = mean Fisher z over ≥ 3 samples.
- **Null:** 1,000 expression-bin-matched random gene sets give an empirical
  one-sided P.

## 3. Code: run order and locations

All commands run from the repo root.

**Upstream steps** (frozen; reruns only reproduce):
- `step2_fetal/`, `step2_benchmark_v2/` (scripts 01–11) and
  `step2_fetal_state/`;
- `step3_cancer/` (Pelka and Joanito pseudobulks; `04_pseudobulk_edger.R`);
- `marker_summary/scripts/build_marker_summary.py`.

**CIOC module** (`core_oncofetal/scripts/`), in order:

| Order | Script | What it does |
|---|---|---|
| 1 | `literature_provenance.py` | A/B/C provenance annotation from `config/literature_audit_31.tsv` |
| 2 | `feature_audit_31.py` | NA-provenance audit for the 31 genes |
| 3 | `panel_contrasts.py` + `panel_edger.R` | Re-runs the gate contrasts with the panel exemption |
| 4 | `build_core_gates.py` | Evidence matrix, CIOC labels, public and restricted tables |
| 5 | `build_gate_funnel.py` | 31-gene funnel workbook (restricted) |
| 6 | `build_genomewide_candidates.py` | Genome-wide H/M/C: 338 and 1,111 genes |
| 7 | `ext_01_prepare_atlases.py` | Khaliq and Che raw AnnData |
| 8 | `ext_02_annotate_compartments.py` | Cluster and marker-detection compartment labels |
| 9 | `ext_03_gate_e_profile.py` | Blinded compartment profile and calibrators |
| 10 | `ext_04_ecos38.py` | ECOS-38 (rules v1.0) |
| 11 | `ext_05_extended_cioc.py` | Extended CIOC v2.0; imports the ext_04 functions; reproduction asserted |
| 12 | `build_references.py` | `results/Methods_references.{md,tsv}` from provenance files, metadata via Crossref (cached) |

**DATA locations** (never committed):
- **Root:** `/Volumes/Stelligen_SSD/Stelligen/DATA`.
- **Panel re-runs:**
  `2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_core_oncofetal_panel_v0.1/`.
- **Atlas work:**
  `2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_extended_cioc_v0.1/`.
  This holds compartment labels, the Gate E profile and coherence strata.
- **New atlases:** `scRNAseq/GSE200997_Khaliq2022/` and
  `scRNAseq/GSE178318_Che2021/`, each with a `link.md` and
  `processed/v0.1/*_raw.h5ad`.
- **Tabula Sapiens LI:** `1.Databases/TabulaSapiens/raw/TS_Large_Intestine.h5ad.zip`.
  The script unzips it to a scratch path; this needs editing in a new
  session (see `TS_LOCAL` in `ext_04_ecos38.py`).
- **Restricted mirror:**
  `2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/core_oncofetal/`.

## 4. Restricted-data rule (Joanito, Synapse syn26844071)

Never commit Joanito-derived material, including:
- C calls and Joanito values;
- **any membership list** (CIOC, 338, Extended CIOC, ECOS-38), because all
  depend on Gate C;
- `Methods.md`, the manuscript text and the funnel workbooks.

These files are git-ignored (see `.gitignore`) and mirrored to the DATA
restricted folder. Public docs give counts and rules only. Before
committing, check that diffs contain no Joanito values or gene lists. All
other work is committed normally, and PRs are used for review.

## 5. Working protocol (user preferences)

**Review loop.** The user pastes reviewer-style feedback, often in Chinese
and CNS-reviewer framed. Implement it exactly, and report back concisely
with numbers.

**Freeze protocol:**
1. Write the rule into the docs and commit it ("freeze(...)") **before**
   computing.
2. Apply it once.
3. Report without tuning, rescue or manual removal.

If a predictable outcome is visible before applying, record it in the
freeze text rather than changing the rule.

**Other rules:**
- Never let downstream CRC, TWEAKR, survival or spatial results feed back
  into membership.
- Report outcomes honestly, including inconvenient ones (for example, coherence
  had no selective power, and CIOC genes fail Gate E).
- Keep `WORKLOG.md` updated at each substantial checkpoint. Commit messages
  end with the Co-Authored-By line.

## 6. Open items / next steps

1. **Methods.md** (`core_oncofetal/Methods.md`; restricted, git-ignored) has
   been corrected. It still needs:
   - journal-style reference numbering (the list is in
     `results/Methods_references.md`);
   - "Supplementary Table X" placeholders replaced;
   - user confirmation that "Gregorieff/TGFB1" means Chen et al. 2023 (Cell
     Stem Cell, "TGFB1 induces fetal reprogramming…"), which is currently
     included.
2. **`core_oncofetal/results/CIOC_manuscript_text.md`** (restricted)
   predates the Extended CIOC and ECOS-38. Sync it with Methods.md if it is
   still used.
3. **RBP1, SPRR1A and EPS8L1** have unresolved literature provenance
   (disclosed; RBP1 is retained in the CIOC by design).
4. **Supportive datasets** (GSE44433, TCGA, organoid cultures) still carry
   old NA labels ("not re-audited"). Run the same feature audit before
   making final figures.
5. **Next phase: validation and application only.** External CRC atlases,
   spatial data (Xenium/Visium), survival, therapy response and
   TWEAKR/YAP. None of these may alter the signatures.
6. **PR chain** #2 → #6 is awaiting review.
