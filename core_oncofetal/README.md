# core_oncofetal — Conserved Intestinal Oncofetal Core (CIOC)

- **Method (v4.0, frozen; last framework revision):** [`docs/Method_core_oncofetal_validated.md`](docs/Method_core_oncofetal_validated.md).
- **Literature nominates; data validate.**
  - The 31 literature-curated genes are the candidate universe (Literature
    candidate = YES).
  - Provenance from a full-text audit (`config/literature_audit_31.tsv` →
    `config/literature_provenance_31.tsv`) is annotation only: A/B/C classes,
    primary studies, resolved or unresolved.
- **Gates** (Core = H ∧ M ∧ C):
  - **H** human developmental, replicated: HGCA ≥ 9 PCW, H-new2 ≥ 9 PCW and
    Gao ≥ 9 W, with ≥ 2/3 fetal-positive and ≥ 1 at log2FC ≥ 0.5 and
    FDR < 0.05;
  - **M** mouse in vivo (GSE230581);
  - **C** CRC replicated (Joanito and Pelka).
- **Panel-gene evaluation (data-QC amendment):** gate contrasts are re-run
  with the 31 panel genes exempt from the genome-wide expression filter
  (`scripts/feature_audit_31.py`, `scripts/panel_edger.R`,
  `scripts/panel_contrasts.py`). NA reasons use a fixed vocabulary; "not
  measured" is no longer used.
- **Calls** are pass, fail or not evaluable. A missing value (NA) is never treated
  as a fail.
- **Former v1.0** (a four-dataset intersection) is kept as a *sensitivity
  analysis using single-dataset hard intersections*, not as an alternative Core.
- **Evidence matrix:** each gene shows all three human log2FC values, the
  dataset(s) giving statistical support, and discordance flags in red.
- **Build:**
  `python3 core_oncofetal/scripts/literature_provenance.py && python3 core_oncofetal/scripts/feature_audit_31.py && python3 core_oncofetal/scripts/panel_contrasts.py && python3 core_oncofetal/scripts/build_core_gates.py && python3 core_oncofetal/scripts/build_gate_funnel.py`.
  Nothing is recomputed.
- **Public outputs:** `results/Core_oncofetal_gate_statistics_public.{xlsx,csv}`
  (31-gene evidence matrix, supportive evidence, literature provenance, long
  table, rules).
- **Restricted (Joanito data-use terms):** the Joanito column, the C call,
  Core labels, the sensitivity set, the 31-gene gate funnel
  (`results/CIOC_gate_funnel_31.xlsx`) and the manuscript text
  (`results/CIOC_manuscript_text.md`) are git-ignored and
  mirrored to `DATA/.../restricted_joanito/core_oncofetal/`.

## Genome-wide fetal–CRC candidates (Level 2) and Extended CIOC (Level 3)

`scripts/build_genomewide_candidates.py` applies the frozen gates H, M and C
genome-wide, removing only the Literature-31 restriction.
- **Outputs:**
  - cross-species fetal–CRC candidates (H ∧ M ∧ C);
  - human conserved fetal–CRC candidates (H ∧ C).
- **Status:** a discovery universe, **not a signature and not for scoring**.
- **Internal validation:** the CIOC is recovered in full, and no other
  Literature-31 gene passes (assertion in the script).
- **Restricted outputs** (Joanito; git-ignored, mirrored to DATA):
  - `Genomewide_fetal_CRC_candidates_calls.csv`
  - `Genomewide_fetal_CRC_candidates.xlsx`
  - `Genomewide_fetal_CRC_candidates.gmt`
- **Final objects** (restricted outputs; see `docs/EXTENDED_CIOC_PLAN.md`).
  Two branches come from the 338 candidates; neither is defined as a
  subset of the other.
  - **Extended CIOC v2.0** (`scripts/ext_05_extended_cioc.py`, rule frozen
    at `72a98b6`):
    - **Rule:** epithelial detectability, plus positive CIOC coherence
      within tumour-derived epithelial cells in both Khaliq and Che.
    - **Use:** cell-resolved epithelial data.
    - **Outputs:** `Extended_CIOC.{xlsx,gmt}`, `Extended_CIOC_calls.csv`.
  - **ECOS-38**, the Epithelial-Compatible Oncofetal Signature
    (`scripts/ext_04_ecos38.py`, rules v1.0 frozen at `4f2d9bd`):
    - **Rule:** compartment compatibility (Gate E) plus stringent
      coherence.
    - **Use:** bulk RNA-seq and unresolved or mixed-cell spatial data.
    - **Description:** epithelial-compatible, not epithelial-specific.
    - **Outputs:** `ECOS_38.{xlsx,gmt}`, `ECOS_38_calls.csv`.
- **State vs lineage specificity:** non-epithelial expression (for example
  SPP1, CCN2, CLU, RBP1) limits mixed-cell scoring. It is not evidence
  against epithelial oncofetal-state membership.
