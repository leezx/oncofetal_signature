# Step 2 — Conserved Fetal-High Intestinal Genes

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan revision
- Origin Date: 2026-10-01
- Verification Status: PLAN REVIEW REQUIRED; no analysis has been run
- Version Label: step2_plan_v2

## Objective

Identify genes that show reproducible fetal-high expression in intestinal
epithelium across two independent human datasets and mouse in vivo intestine.

This step is deliberately narrow. It does not attempt to reconstruct the full
trajectory of intestinal development or model every source of technical
variation. It tests a simple evidence chain:

```text
HGCA human discovery
        ∩
Gao independent human replication
        ∩
mouse in vivo developmental evidence
        =
Conserved Fetal-High Candidates
```

The analysis stops when this set is frozen. It must not be optimized using CRC
expression, TNFRSF12A/TWEAKR association, or downstream manuscript results.

## Scope and non-goals

The primary analysis:

- includes epithelial cells only;
- uses donors or animals, never individual cells, as biological replicates;
- estimates fetal-versus-adult effects independently within each dataset;
- uses an adult proliferative epithelial compartment as a negative control;
- uses an in vivo mouse dataset for the required cross-species evidence;
- uses Ensembl one-to-one orthologues for the conserved set.

The primary analysis does not require:

- cross-dataset batch correction or joint integration;
- region-specific statistical gates;
- leave-one-donor-out analysis;
- developmental-stage regression;
- human–mouse developmental-time conversion;
- effect-size meta-analysis or rank aggregation;
- a second mouse validation cohort;
- rescue of one-to-many orthologues.

Organoid and fetal spheroid datasets may be shown as optional validation but
cannot establish primary developmental evidence.

## Dataset roles

| Role | Dataset | Primary comparison | Evidence rule |
|---|---|---|---|
| Human discovery | Elmentaite et al. 2021 Human Gut Cell Atlas (HGCA) | fetal vs adult intestinal epithelium | H1 discovery plus adult proliferative control |
| Independent human validation | Gao et al. 2018; fetal GSE95630 and adult GSE103154 | fetal vs adult large-intestinal epithelium | H2 direction and effect size; no FDR gate |
| Mouse primary | one qualifying in vivo fetal/adult intestinal epithelial dataset | fetal vs adult intestinal epithelium | M1 effect size and FDR |
| Culture validation | fetal/adult organoid or spheroid data | dataset-specific | optional annotation only; never a gate |

HGCA is the discovery dataset because it contains fetal, paediatric, and adult
human gut cells from multiple anatomical regions. Gao is analyzed separately as
an external replication cohort; its adult arm contains only two individuals,
so it is not treated as a second high-powered discovery analysis.

## Preflight requirements

Before differential expression, create and freeze a dataset manifest containing:

- accession, source URL, assay, genome build, and matrix provenance;
- donor/animal ID, library ID, developmental stage, age, region, chemistry, and
  batch when available;
- original and harmonized epithelial labels;
- number of retained epithelial cells per biological replicate;
- exclusions and reasons;
- availability of raw counts;
- source-file checksums.

Cell inclusion and label harmonization must be completed without reference to
the candidate-gene results. Paediatric samples are excluded from the primary
fetal-versus-adult contrast.

### Epithelial inclusion

Restrict all primary comparisons to high-confidence epithelial cells after
dataset-appropriate QC and doublet removal. Confirm epithelial identity using
multiple epithelial markers and exclude immune, endothelial, and mesenchymal
cells using multi-gene evidence.

Do not use TACSTD2, CLU, ANXA1, LY6A, TNFRSF12A, or other proposed oncofetal
markers as required inclusion markers.

Each human donor or mouse animal must contribute at least 50 retained epithelial
cells to a single-cell pseudobulk. This is a pragmatic eligibility threshold,
not a claim that 50 cells guarantees statistical power.

## H1 — human discovery in HGCA

### Pseudobulk unit

Aggregate raw counts to:

```text
donor × developmental stage
```

for all eligible epithelial cells. If one donor contributes multiple intestinal
regions or technical libraries, combine them within that donor for the primary
analysis. A donor × region observation must not be treated as an independent
biological replicate.

Region-specific results may be reported descriptively as a sensitivity analysis
but do not determine H1.

### Primary discovery comparison

Compare:

```text
fetal intestinal epithelium vs adult intestinal epithelium
```

using edgeR quasi-likelihood differential expression on donor-level raw-count
pseudobulks with:

```text
~ developmental_stage
```

as the default design. Add a technical covariate only when it is not perfectly
confounded with developmental stage and the available donor count supports its
estimation. Do not use normalized, integrated, batch-corrected, or imputed
expression values for differential expression. Apply `edgeR::filterByExpr`
using the donor-level developmental-stage groups before model fitting.

The primary HGCA fetal-high criterion is:

- log2FC ≥ 0.5; and
- Benjamini–Hochberg FDR < 0.05.

### Adult proliferative epithelial negative control

Construct a single adult proliferative epithelial compartment by combining
annotated adult:

- stem cells;
- progenitor cells;
- transit-amplifying cells;
- cycling epithelial cells.

Aggregate its raw counts at the donor level and compare fetal epithelium with
adult proliferative epithelium. This contrast asks whether a gene is fetal-high
rather than merely high in proliferating or crypt-progenitor cells.

The negative-control criterion is directional:

- fetal/adult-proliferative log2FC > 0.

It is not required to meet a separate FDR threshold. Separate adult stem, TA,
progenitor, or cycling comparisons may be shown descriptively but are not hard
gates.

A gene passes **H1** when it passes both the primary HGCA fetal-high criterion
and the adult proliferative directional control.

## H2 — independent human validation in Gao et al. 2018

Restrict validation to:

- fetal large-intestinal epithelial cells; and
- adult large-intestinal epithelial cells.

Exclude oesophagus and stomach. Aggregate at the biological donor/sample level
where the metadata permit. Analyze Gao independently from HGCA; do not batch
correct or integrate the two datasets.

Because the adult reference contains only two individuals, statistical
significance is not a hard replication requirement. A gene passes **H2** when:

- the overall fetal/adult log2FC is ≥ 0.5; and
- fetal expression is directionally higher than each adult donor separately.

The operational donor-consistency check must be frozen before testing. By
default, the median normalized fetal pseudobulk expression must exceed the
normalized expression of each adult donor. Report any nominal P value or FDR,
but do not use it to rescue or reject H2.

If donor-level raw-count pseudobulk cannot be reconstructed, label this result
`effect-only validation` and document the expression scale used.

## M1 — mouse in vivo developmental evidence

Select one mouse dataset that provides:

- freshly isolated in vivo fetal and adult intestinal epithelium;
- identifiable biological replicates in both stages;
- raw gene-level counts;
- no injury, tumour, treatment, transgene activation, organoid culture, or
  spheroid culture in the primary comparison.

Purified epithelial bulk RNA-seq and epithelial scRNA-seq with animal-level
pseudobulk are both acceptable. A well-replicated purified epithelial bulk
dataset is not downgraded merely because it is not single-cell.

For scRNA-seq, combine eligible regions and technical libraries within each
animal for the primary pseudobulk. For sorted bulk RNA-seq, use the animal-level
libraries directly after confirming that technical replicates are not counted
as independent animals.

Compare:

```text
mouse fetal vs adult in vivo intestinal epithelium
```

A one-to-one mouse orthologue passes **M1** when:

- log2FC ≥ 0.5; and
- Benjamini–Hochberg FDR < 0.05.

If no existing dataset meets the admission requirements, stop at dataset
selection or acquisition. Do not promote an organoid or spheroid dataset to
primary evidence.

## Cross-species mapping and final definition

Map human and mouse genes using Ensembl one-to-one orthologues from a frozen,
documented Ensembl release.

Genes without a one-to-one orthologue are excluded from the primary conserved
set but may be retained in a separate species-specific table. No one-to-many
rescue, cross-species meta-analysis, weighted score, or rank aggregation is
required for this step.

A gene is a **Conserved Fetal-High Candidate** only when:

1. it passes HGCA discovery and the adult proliferative control (**H1**);
2. it passes Gao directional/effect-size replication (**H2**); and
3. its one-to-one mouse orthologue passes in vivo mouse evidence (**M1**).

Therefore:

```text
Conserved Fetal-High = H1 ∩ H2 ∩ M1
```

Additional labels are retained for auditability:

| Label | Definition |
|---|---|
| Conserved Fetal-High Candidate | passes H1, H2, and M1 |
| Conserved discovery holdout | passes H1 and M1 but not H2; excluded from the strict set |
| Human fetal-high only | passes H1 and H2 but not M1 or lacks a one-to-one orthologue |
| Proliferation-associated reject | passes primary HGCA fetal/adult criterion but has fetal/adult-proliferative log2FC ≤ 0 |
| Species-specific developmental candidate | passes within species but has no one-to-one orthologue |
| Culture-supported | optional annotation added after primary classification; never changes pass/fail status |

## Minimal diagnostics and hard stops

The following diagnostics are required:

- donor/animal-level expression plots for final candidates and prespecified
  markers;
- pseudobulk library sizes and retained cell counts for single-cell datasets;
- PCA/MDS of pseudobulks colored by stage and known technical variables;
- the number of donors/animals contributing to every contrast;
- region composition per donor in HGCA;
- enrichment of canonical S/G2M genes among the proliferation-associated
  rejects;
- explicit display of TACSTD2, CLU, ANXA1, LY6A, and TNFRSF12A regardless of
  whether they pass.

Stop and revise the affected contrast when:

- developmental stage is perfectly confounded with an unmodelled technical
  batch;
- raw counts or biological replicate identifiers are unavailable for HGCA or
  the mouse primary dataset;
- fewer than three usable biological replicates per stage are available for
  HGCA or mouse primary evidence;
- epithelial labels cannot be defended without using the target genes;
- Gao cannot be separated into fetal and adult large-intestinal epithelial
  samples or its adult donors cannot be distinguished.

Do not replace a stopped pseudobulk analysis with cell-level Wilcoxon testing.

## Multiple testing and reporting

- Apply BH correction over all genes tested in each primary discovery dataset.
- Analyze HGCA, Gao, and mouse independently; do not pool their P values.
- Treat H1 ∩ H2 ∩ M1 as a deterministic evidence rule.
- Report effect size even when a gene fails a gate.
- Define all thresholds and labels before examining CRC or TWEAKR-associated
  results.

## Required output

Generate one auditable gene-level table containing at least:

| Column | Meaning |
|---|---|
| `gene` | human gene symbol/stable identifier |
| `HGCA_log2FC` | fetal vs adult epithelial effect |
| `HGCA_FDR` | BH-adjusted discovery P value |
| `HGCA_proliferative_control_log2FC` | fetal vs adult proliferative effect |
| `H1_pass` | HGCA discovery plus negative-control result |
| `Gao_log2FC` | fetal vs adult large-intestinal effect |
| `Gao_direction_consistent` | fetal exceeds each adult donor |
| `H2_pass` | independent human replication result |
| `mouse_gene` | mapped one-to-one mouse orthologue |
| `Mouse_log2FC` | mouse fetal vs adult effect |
| `Mouse_FDR` | BH-adjusted mouse P value |
| `M1_pass` | mouse in vivo result |
| `one_to_one_orthologue` | Ensembl relationship status |
| `final_label` | final evidence category |
| `Final_conserved_fetal_high` | strict H1 ∩ H2 ∩ M1 flag |

Also produce:

- the dataset manifest and source checksums;
- the frozen epithelial label map and sample-inclusion table;
- pseudobulk QC tables and plots;
- complete HGCA, Gao, and mouse effect tables;
- the frozen Ensembl orthologue map and release;
- scripts, package versions, commands, and session information needed for a
  clean rerun.

## Repository separation

All executable analysis, downloaded data, intermediate objects, and large
results must live under:

`/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/oncofetal_signature/step2_developmental_axis/`

This GitHub repository receives only lightweight, reviewable artifacts:
analysis plans, source code, environment lock files, provenance manifests,
final tables, and figure files. It must not contain raw data, serialized
single-cell objects, or large intermediate matrices.

## Execution sequence and checkpoints

1. **Checkpoint A — dataset admission:** approve manifests, the mouse primary
   dataset, biological replicate counts, and the confounding assessment.
2. **Checkpoint B — annotation freeze:** approve epithelial inclusion and the
   combined adult proliferative label map.
3. Construct donor/animal-level pseudobulks and publish QC.
4. **Checkpoint C — model freeze:** approve the simple within-dataset models,
   thresholds, and the Gao donor-consistency operation.
5. Run and freeze H1, H2, and M1 independently.
6. Apply the frozen one-to-one orthologue and deterministic intersection rules.
7. **Checkpoint D — result freeze:** approve the auditable candidate table.
8. Stop Step 2 and proceed to the CRC-high axis; do not tune the developmental
   definition using downstream results.

## Sources verified for this plan

- Elmentaite R. et al. *Cells of the human intestinal tract mapped across space
  and time*. Nature 597, 250–255 (2021).
  <https://doi.org/10.1038/s41586-021-03852-1>
- Human Gut Cell Atlas download portal.
  <https://www.gutcellatlas.org/>
- Gao S. et al. *Tracing the temporal-spatial transcriptome landscapes of the
  human fetal digestive tract using single-cell RNA-sequencing*. Nature Cell
  Biology 20, 721–734 (2018).
  <https://doi.org/10.1038/s41556-018-0105-4>
- Gao fetal GEO series GSE95630.
  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95630>
- Gao adult GEO series GSE103154.
  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103154>
