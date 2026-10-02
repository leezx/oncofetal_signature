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
- **Calls** are pass, fail or not evaluable. Not measured is never treated
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
