# core_oncofetal — Conserved Intestinal Oncofetal Core

- **Method (frozen v3.0, final framework):** [`docs/Method_core_oncofetal_validated.md`](docs/Method_core_oncofetal_validated.md).
- **Rule:** a literature-anchored evidence framework; a Core gene must pass
  all four axes:
  - **L** literature provenance from a full-text audit of all 31 genes
    (`config/literature_audit_31.tsv` → `config/literature_provenance_31.tsv`);
  - **H** human developmental, replicated: HGCA ≥ 9 PCW, H-new2 ≥ 9 PCW and
    Gao ≥ 9 W, with ≥ 2/3 fetal-positive and ≥ 1 at log2FC ≥ 0.5 and
    FDR < 0.05;
  - **M** mouse in vivo (GSE230581);
  - **C** CRC replicated (Joanito and Pelka).
- **Calls** are pass, fail or not evaluable. Not measured is never treated
  as a fail.
- **Former v1.0** (a four-dataset intersection) is kept as the *stringent
  intersection set* sensitivity analysis.
- **Build:**
  `python3 core_oncofetal/scripts/literature_provenance.py && python3 core_oncofetal/scripts/build_core_gates.py`.
  Nothing is recomputed.
- **Public outputs:** `results/Core_oncofetal_gate_statistics_public.{xlsx,csv}`
  (31-gene evidence matrix, supportive evidence, literature provenance, long
  table, rules).
- **Restricted (Joanito data-use terms):** the Joanito column, the C call,
  Core labels and the stringent set are in the git-ignored full workbook,
  mirrored to `DATA/.../restricted_joanito/core_oncofetal/`.
