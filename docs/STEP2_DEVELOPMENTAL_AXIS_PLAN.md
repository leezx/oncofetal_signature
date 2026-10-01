# Step 2 analysis plan — a conserved intestinal developmental axis

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: plan
- Origin Date: 2026-10-01
- Verification Status: PLAN REVIEW REQUIRED; no analysis has been run
- Version Label: step2_plan_v1

## Decision to be made

Identify genes that are genuinely higher in **in vivo fetal intestinal
epithelium** than in adult intestinal epithelium in both human and mouse, while
excluding genes whose apparent fetal enrichment is explained by the adult
crypt/stem/transit-amplifying (TA) compartment.

The primary product is a table of **Conserved Fetal-High Candidates** with an
auditable pass/fail result for every required contrast. Organoid and spheroid
data are validation datasets only and cannot establish developmental status.

## Biological hypothesis

A conserved intestinal developmental program exists that is enriched in fetal
epithelium, is attenuated in both differentiated and proliferative adult
epithelial compartments, and is directionally reproducible across species and
independent human cohorts.

## Scope and non-goals

This step establishes only the developmental axis. It does not:

- use CRC expression to select genes;
- optimize a score against TNFRSF12A/TWEAKR or any favored marker;
- infer a fetal program from organoids, spheroids, injury, YAP activation, or
  published gene lists;
- treat cells as biological replicates;
- merge studies before estimating within-study effects.

CRC enrichment, tissue specificity, literature support, and gene-module
coherence are later evidence axes.

## Dataset roles

| Role | Dataset | Required material | Planned use |
|---|---|---|---|
| Human discovery | Elmentaite et al. 2021, Space-Time Gut Cell Atlas; ArrayExpress E-MTAB-9543, E-MTAB-9536, E-MTAB-9532, E-MTAB-9533, E-MTAB-10386, and E-MTAB-8901 where relevant | raw UMI/count matrix plus donor, age, region, chemistry/library, and cell-type metadata | Primary fetal-versus-adult contrasts and adult proliferative negative control |
| Human independent validation | Gao et al. 2018; fetal GSE95630 and adult large intestine GSE103154 (SuperSeries GSE103239) | processed counts/TPM and barcode metadata; raw data only if necessary | Large-intestine-only replication. Adult n=2 means emphasis is effect direction, magnitude, and donor consistency rather than a standalone significance claim |
| Mouse primary | To be selected from the existing data inventory | in vivo, epithelial-resolved fetal and adult intestine; raw counts and biological replicate metadata | Required cross-species developmental contrast |
| Culture validation | Fetal/adult organoid or spheroid datasets, including GSE228519 subseries where useful | expression matrix and replicate metadata | Sensitivity/biological interpretation only; never a primary gate |

The Human Gut Cell Atlas reports more than 428,000 cells across fetal,
paediatric, and adult gut regions and supplies downloadable raw and normalized
objects. Gao et al. profile fetal digestive tract from 6–25 gestational weeks;
the paired adult GEO series contains 1,463 QC-passing cells from two adult large
intestines.

### Mouse dataset admission gate

The mouse dataset must satisfy all of the following before analysis:

1. freshly isolated **in vivo** intestinal epithelium at both fetal and adult
   stages;
2. at least three independent biological replicates per stage, or a documented
   limitation with no cell-level pseudo-replication;
3. raw gene-level counts and replicate-level metadata available;
4. comparable anatomical region(s), with small and large intestine kept
   separate when possible;
5. no treatment, injury, tumour, transgene activation, or culture as the
   developmental comparison;
6. epithelial identity can be established without selecting on the candidate
   genes being tested.

If no existing asset passes this gate, execution stops for dataset selection or
acquisition. An organoid dataset must not be promoted to primary evidence.

## Preflight audit and frozen metadata

Before differential expression (DE), create a dataset manifest containing:

- accession, paper, assay, genome build, and matrix provenance;
- donor/animal ID, library ID, developmental stage, exact age, anatomical
  region, sex when available, chemistry, and batch;
- original and harmonized epithelial cell labels;
- cells per donor/animal × region × compartment;
- exclusions and their reasons;
- whether raw counts, normalized values, or both are present;
- checksums for downloaded source files.

The manifest and a frozen cell-inclusion table are required inputs. Cell labels
must be harmonized without reference to the eventual candidate-gene results.
Paediatric samples are descriptive and are excluded from the primary fetal vs
adult contrast.

### Epithelial inclusion

Retain high-confidence epithelial cells after dataset-appropriate QC and
doublet removal. Confirm epithelial identity using a panel (for example EPCAM,
KRT8, KRT18, KRT19) and exclude immune, endothelial, and mesenchymal
contamination using multi-gene panels. Do not require TACSTD2, CLU, ANXA1,
LY6A, TNFRSF12A, or other proposed oncofetal markers for inclusion.

Adult epithelium is partitioned before testing into:

1. **adult proliferative**: crypt stem, cycling stem/progenitor, and TA cells;
2. **adult differentiated**: absorptive and secretory epithelial lineages;
3. **adult all-epithelium**: the union of adult epithelial cells.

Original annotations are preserved, and every harmonized label must map back to
an original label. Ambiguous cells are retained only in all-epithelium
sensitivity analyses, not in compartment-specific gates.

## Statistical unit and pseudobulk construction

The biological replicate is the donor (human) or animal (mouse), not the cell
or sequencing library.

1. Sum raw counts within donor/animal × broad anatomical region × developmental
   stage × epithelial compartment.
2. Combine technical libraries from the same biological replicate before DE.
3. Require a pre-specified minimum of 30 retained cells per pseudobulk and at
   least 10 detected counts for a gene in enough biological replicates to make
   the contrast estimable. Report sensitivity at 20 and 50 cells.
4. Do not pseudobulk normalized or integrated expression values.
5. Never use batch-corrected embeddings or imputed expression for DE.

When a donor contributes multiple regions, repeated observations must not be
treated as independent. The primary analysis will either sum within a matched
broad region per donor or use a repeated-measures model; the exact choice is
frozen after the metadata audit and before examining candidate results.

## Primary contrasts

Estimate effects independently within each dataset. Small intestine and large
intestine/colon are analyzed separately first; a pooled regional estimate is
secondary and must model region.

For each eligible region, fit a count-aware pseudobulk model (edgeR
quasi-likelihood is the default) with developmental stage as the coefficient of
interest and known technical covariates included only when estimable. Report
log2 fold change, 95% confidence interval, raw P value, and Benjamini–Hochberg
FDR.

Required contrasts are:

- fetal epithelium vs adult all-epithelium;
- fetal epithelium vs adult differentiated epithelium;
- fetal epithelium vs adult proliferative epithelium.

Cell-cycle scores are reported as diagnostics, not regressed from the primary
model. The proliferative-compartment contrast is the biological negative
control. A cell-cycle-regressed result may be shown only as sensitivity analysis.

## Candidate gates

Thresholds are frozen before DE and applied gene-by-gene without exceptions.

### Gate H1 — human discovery (HGCA)

A gene passes H1 only when:

- log2FC ≥ 0.5 and FDR < 0.05 for all three required contrasts in at least one
  anatomically matched intestinal region;
- the log2FC is positive in every other estimable human region;
- the leave-one-donor-out log2FC remains positive; and
- no single donor contributes more than 50% of total pseudobulk counts for that
  gene within either stage.

### Gate H2 — independent human validation (Gao)

Because only two adult donors are available, H2 is a replication gate rather
than a second discovery test. A gene passes when:

- fetal large-intestinal epithelium has log2FC ≥ 0.5 versus adult large-
  intestinal epithelium;
- the effect is positive against each adult donor separately; and
- the direction is positive in the adult proliferative comparison when the
  published cell labels permit that contrast.

Nominal P values/FDR are reported but are not used to rescue a weak or
discordant effect. If raw-count-compatible pseudobulk cannot be reconstructed,
the result is explicitly labeled `effect-only validation`.

### Gate M1 — mouse in vivo

The one-to-one mouse orthologue must show log2FC ≥ 0.5 and FDR < 0.05 for fetal
versus adult in vivo epithelium, remain positive in leave-one-animal-out
analysis, and—when adult annotations allow—remain positive versus adult
proliferative and adult differentiated compartments separately.

### Final labels

| Label | Definition |
|---|---|
| Conserved Fetal-High Candidate | passes H1, H2, and M1 |
| Conserved discovery-only holdout | passes H1 and M1 but fails or is untestable in H2; excluded from the final strict set |
| Human fetal-high only | passes H1/H2 but not M1 or lacks a one-to-one orthologue |
| Proliferation-associated reject | fetal-high vs adult all/differentiated but fails fetal vs adult proliferative |
| Culture-supported | annotation added after primary classification; cannot change a fail to pass |

Human–mouse mapping uses Ensembl one-to-one orthologues frozen to a documented
release. Ambiguous one-to-many mappings are reported but excluded from the
strict conserved set.

## Robustness and falsification analyses

The following are mandatory and cannot be substituted by a larger cell count:

- donor/animal-level expression plots with all biological replicates visible;
- leave-one-donor/animal-out effect estimates;
- small-intestine and colon effects shown separately;
- early and late fetal strata where sample numbers permit;
- downsampling cells per donor to test cell-number dominance;
- label sensitivity using narrow versus broad adult proliferative definitions;
- comparison of results with and without ambiguous epithelial cells;
- enrichment of canonical S/G2M genes among passes and rejects;
- detection of chemistry/batch–stage confounding;
- comparison with culture datasets only after the primary list is frozen.

Hard stop conditions:

- stage is perfectly confounded with an unmodelled processing batch in the
  proposed primary contrast;
- fewer than three usable biological replicates exist in an HGCA stratum or the
  selected mouse dataset;
- raw counts or donor/animal identifiers are unavailable for a primary dataset;
- epithelial or adult proliferative labels cannot be defended independently of
  the target genes.

Under a hard stop, report the limitation and redesign the contrast; do not fall
back to cell-level Wilcoxon tests.

## Multiple testing and reporting rules

- BH correction is performed over all tested genes separately for each primary
  dataset × region × contrast family.
- The conjunction across contrasts/datasets is a deterministic evidence rule;
  P values are not pooled across studies with incompatible assays.
- Report effect sizes and confidence intervals even for genes failing FDR.
- Show all proposed markers (including TACSTD2, CLU, ANXA1, LY6A and
  TNFRSF12A) in the audit table regardless of outcome, without overriding gates.
- Freeze thresholds, label mappings, and exclusions before viewing final gene
  ranks.

## Planned implementation and separation of repositories

All executable analysis, downloaded data, intermediate objects, and large
results must live under:

`/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/oncofetal_signature/step2_developmental_axis/`

Proposed structure:

```text
step2_developmental_axis/
├── raw/                 # immutable source files and checksums
├── metadata/            # manifests, label maps, inclusion tables
├── processed/           # filtered objects and pseudobulk counts
├── scripts/             # numbered, reproducible scripts
├── results/
│   ├── tables/
│   ├── figures/
│   └── qc/
└── logs/                # commands, package versions, session information
```

This GitHub repository receives only lightweight, reviewable artifacts:
analysis plans, source code, environment lock files, final tables, figure files,
and provenance manifests. It must not contain raw data, serialized single-cell
objects, or large intermediate matrices.

## Required deliverables

| Deliverable | Minimum content | Acceptance criterion |
|---|---|---|
| Dataset manifest | accessions, source URLs, checksums, assay/build, replicate and region counts | every analyzed matrix traceable to a source |
| Cell-label map | original label, harmonized label, evidence/rule | no target-gene-based selection |
| Pseudobulk QC | cells, library size, detected genes per replicate/stratum | exclusions documented before DE |
| Within-dataset DE tables | gene ID/symbol, base expression, log2FC, CI, P, FDR, contrast, region | one row per tested gene/contrast; no cell-level P values |
| Orthologue map | human and mouse stable IDs, relationship type, Ensembl release | strict set uses one-to-one only |
| Candidate audit table | H1/H2/M1 subcriteria and final label | every gate is machine-readable and independently checkable |
| Robustness report | leave-one-out, region, cell-number, label, and batch sensitivities | conclusions stable or limitations explicit |
| Figures | replicate-level effects, contrast UpSets, heatmap/forest plot | donor/animal values visible; no cell-count-inflated error bars |
| Reproducibility record | commands, package versions, seeds, file hashes | clean rerun produces the same candidate table |

## Planned figure set

1. Dataset and contrast schematic with discovery/validation roles.
2. Pseudobulk sample map and QC by donor/animal, stage, region, and compartment.
3. HGCA effect-size concordance across adult all, differentiated, and
   proliferative contrasts.
4. Cross-species human–mouse log2FC plot with Gao validation status.
5. Donor-level forest/strip plots for all final candidates and prespecified
   markers.
6. Candidate flow diagram showing removals at H1, proliferation control, H2,
   orthology, and M1.

## Execution sequence and review checkpoints

1. **Checkpoint A — dataset admission:** approve manifests, mouse primary
   dataset, replicate counts, and confounding assessment.
2. **Checkpoint B — annotation freeze:** approve epithelial inclusion and adult
   proliferative/differentiated label map.
3. Construct pseudobulks and publish QC; no gene ranking yet.
4. **Checkpoint C — model freeze:** approve contrasts, covariates, filtering,
   and thresholds.
5. Run HGCA H1, then lock the human discovery table.
6. Run Gao H2 and mouse M1 independently.
7. Apply the deterministic conjunction and orthology rules.
8. Run robustness/falsification analyses and generate deliverables.
9. **Checkpoint D — interpretation:** approve the final strict and holdout
   tables before any CRC-axis analysis.

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
