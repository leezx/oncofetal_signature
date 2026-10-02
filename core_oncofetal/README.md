# core_oncofetal — Conserved Intestinal Oncofetal Core

- **Method (frozen v1.0):** [`docs/Method_core_oncofetal_validated.md`](docs/Method_core_oncofetal_validated.md).
- **Gate chain:** Literature 31 → G1 HGCA ≥ 9 PCW (human in vivo) → G2
  GSE230581 (mouse in vivo) → G3 Joanito (CRC) → G4 Pelka (CRC replication).
- **Pass rule:** measured, log2FC ≥ 0.5, FDR < 0.05.
- **Supportive evidence** is reported and never used for selection.
- **Build:** `python3 core_oncofetal/scripts/build_core_gates.py`. It reads
  the marker summary and recomputes nothing.
- **Public outputs:** `results/Core_oncofetal_gate_statistics_public.{xlsx,csv}`
  and `results/Core_oncofetal_gate_long_public.csv`.
  - The workbook has one column per gate (log2FC; P; FDR | call), with
    supportive evidence on separate sheets.
  - It omits the Joanito G3 column and the final Core label.
- **Restricted outputs (Joanito data-use terms):** the full workbook with G3
  and Core membership is git-ignored. It is mirrored to
  `DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/core_oncofetal/`.
- **Developmental Core (G1 ∧ G2, public):** GJA1, CLU, ANXA6, SPP1, RBP1.
