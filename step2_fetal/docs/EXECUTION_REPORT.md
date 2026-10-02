# Step 2 execution report

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: execution
- Origin Date: 2026-10-01
- Verification Status: VERIFIED (download checks, metadata audit, scripted analysis)
- Version Label: step2_developmental_axis_v0.1

## Frozen implementation

The approved gates in `ANALYSIS_PLAN.md` were applied without
adding integration, trajectory, region-specific differential expression, GSEA,
or additional validation cohorts.

### H1: HGCA discovery

- Source object: `epi_raw_counts02_v2.h5ad` from the Human Gut Cell Atlas portal.
- Upstream and local MD5: `2a149b8cf04567569707e9d1fab27209`.
- Object dimensions: 142,113 cells × 33,538 genes.
- Inclusion: author-labelled epithelial cells from fetal or healthy adult
  intestinal regions; lymph node and paediatric samples excluded.
- Biological replicate: `Sample name` (donor), with all eligible regions and
  libraries combined.
- Eligibility: at least 50 retained cells per donor.
- Primary pseudobulk: 16 fetal donors and 7 adult donors.
- Adult proliferative pseudobulk: the same fetal donors versus 7 adult donors,
  using author labels `Stem cells`, `Proximal progenitor`, `Distal progenitor`,
  `TA`, and `Progenitor (NEUROG3+)`.
- Model: edgeR quasi-likelihood, `~ stage`, TMM normalization and
  `filterByExpr`; fetal is the positive coefficient.
- Result: 4,733 genes passed log2FC ≥ 0.5 and FDR < 0.05; 4,645 also had
  fetal/adult-proliferative log2FC > 0 and passed H1.

### H2: Gao 2018 independent validation

- Fetal source: GSE95630 TPM matrix plus the corrected Supplementary Table 2.
- Adult source: GSE103154 UMI-normalized TPM matrix and GEO library metadata.
- Inclusion: author-labelled fetal large-intestinal epithelial cells (849
  cells) and all adult large-intestinal crypt/villus epithelial cells.
- Replicates: 12 fetal embryos after combining additional libraries from the
  same age/embryo identifier; adult donors P1 and P2.
- Expression: mean TPM/normalized TPM within each donor. `Gao_log2FC` is
  `log2((median fetal donor expression + 1) / (median adult donor expression + 1))`.
- Direction check: the fetal donor median must exceed P1 and P2 separately.
- This is explicitly effect-only validation; no significance gate was applied.
- Result: 5,271 genes passed log2FC ≥ 0.5 and both donor-direction checks.

### M1: mouse in vivo evidence

- Source: GSE230581 `in_vivo_counts`, already present in the project archive
  and copied into the managed bulk-RNA-seq dataset directory.
- Comparison: three E16.5 proximal small-intestinal epithelial replicates
  versus three adult proximal small-intestinal crypt epithelial replicates.
- Model: edgeR quasi-likelihood, `~ stage`, TMM normalization and
  `filterByExpr`; fetal is the positive coefficient.
- Result: 3,140 of 13,314 tested genes passed log2FC ≥ 0.5 and FDR < 0.05.

### Cross-species intersection

- Mapping: Ensembl release 116 BioMart, `ortholog_one2one` only.
- The exact XML query and raw response are frozen under the DATA hierarchy.
- HGCA was joined to Ensembl by stable human gene ID; Gao was joined by the
  HGCA gene symbol; mouse evidence was joined by Ensembl mouse gene symbol.
- Result: 13,397 HGCA-tested genes had a one-to-one mapping; 706 passed H1,
  H2, and M1 and are written to `results/tables/Conserved_Fetal_High.csv`.

## DATA locations

- HGCA: `/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/HGCA_Elmentaite2021`
- Gao: `/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq/GSE103239_Gao2018`
- Mouse: `/Volumes/Stelligen_SSD/Stelligen/DATA/bulkRNAseq/GSE230581`
- Ensembl mapping: `/Volumes/Stelligen_SSD/Stelligen/DATA/1.Databases/Ensembl_orthologues/release_116`
- Versioned working results: `/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-01_step2_developmental_axis_v0.1`

## Reproduction order

1. `scripts/00_download_inputs.sh`
2. `scripts/01_hgca_pseudobulk.py`
3. `scripts/02_hgca_fetal_high.R`
4. `scripts/03_gao_validation.py`
5. `scripts/04_mouse_fetal_high.R`
6. `scripts/05_cross_species_intersection.py`
7. `scripts/06_plot_volcano.R`

The complete command sequence is encoded in `scripts/run_step2.sh` and explained
in `docs/RUNBOOK.md`.
