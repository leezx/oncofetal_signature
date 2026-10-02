# Extended CIOC — plan (rules v0.1 PROPOSED; not frozen, not applied)

## Architecture

| Level | Name | Definition | Status |
|---|---|---|---|
| 1 | **CIOC** | Literature-31 candidates passing Gates H, M, C (method v4.0 + data-QC amendment) | frozen |
| 2 | **Genome-wide fetal–CRC candidates** | Same frozen gates applied to every gene tested in ≥ 1 gate contrast | computed; frozen as the discovery universe |
| 2a | Cross-species fetal–CRC candidates | H ∧ M ∧ C | discovery universe for Level 3 |
| 2b | Human conserved fetal–CRC candidates | H ∧ C (mouse ignored; no cross-species evidence) | reported only |
| 3 | **Extended CIOC** | Level 2a ∩ Gate E (epithelial compatibility) ∩ CIOC program coherence | rules to be frozen before computation |

- Level 2 is **not a signature** and is not used for scoring. The name
  "Extended CIOC" is reserved for Level 3.
- Level 2 is produced by `scripts/build_genomewide_candidates.py`. Its
  membership is restricted (Gate C uses Joanito).
- **Internal validation.** Level 1 and Level 2 start from different
  nomination universes (31 literature genes vs every tested gene). The
  genome-wide search recovers all eight CIOC genes, and no other
  Literature-31 candidate passes H, M and C.

## Why Level 3 needs more than DE gates

Level 2a contains genes with strong, concordant H/M/C support from collagen,
smooth-muscle, endothelial and lymphoid programs. Examples are COL1A1, SPARC,
COL6A1, CALD1, MYL9, NOTCH3, PDGFB and LCK. DE gates cannot tell whether
the signal comes from:
- non-epithelial cells retained in "epithelial" pseudobulks or annotations;
- or mesenchymal-like transcription in developing or malignant epithelium.

Level 2a also mixes proliferation, translation, YAP, wound, EMT,
metabolism, stress and inflammatory programs. Two further requirements
address this:
- epithelial compatibility (Gate E);
- co-variation with the CIOC state.

Genes are not removed by hand. H, M and C are not modified.

## Gate E — epithelial compatibility (design stage)

**Intent:** exclude genes whose apparent signal is predominantly
attributable to non-epithelial compartments. It does **not** require
epithelial specificity, because shared programs (for example SPP1) are
legitimate.

**Candidate metrics, per gene, from an atlas with all compartments:**
- detection in epithelial (malignant or fetal) cells;
- epithelial share of expression across epithelial, fibroblast/stromal,
  endothelial, myeloid and lymphoid compartments.

Each metric is computed on donor-level compartment means.

**Order of work** (thresholds are not chosen in advance):
1. Plot the epithelial-vs-non-epithelial distribution for the Level 2a
   genes.
   - Canonical lineage markers are shown as calibrators: epithelial EPCAM,
     KRT8, CDH1; fibroblast COL1A2, DCN; endothelial PECAM1, VWF; immune
     PTPRC.
   - Candidate genes are not labelled on the plot.
2. Set the rule from technical detectability and the calibrators.
3. Freeze the rule.
4. Apply the rule.

**Reference requirement (open decision).** The reference must contain
epithelial and non-epithelial cells. It must also represent the states in
which these genes are expressed: malignant and/or fetal epithelium. It
should be independent of the gate datasets.
- **Locally available:** Tabula Sapiens Large Intestine is independent and
  has all compartments. It is **adult normal**, however. Oncofetal genes are
  by definition low in adult epithelium, so in this atlas genuine
  oncofetal-epithelial genes would look non-epithelial. The atlas can
  therefore serve only as a secondary check, not as the reference.
- **Not local:**
  - An independent CRC atlas with all compartments and non-overlapping
    patients. Joanito includes the SMC and KUL3 cohorts, so Lee 2020 is not
    independent. Options are Khaliq 2022 (GSE200997), Che 2021 (GSE178318)
    and Becker 2022 (HTAN, GSE201348). The accessions still need to be verified before download.
  - A fetal atlas with all compartments.

## CIOC program coherence (design stage)

**Question:** when the eight-gene CIOC state is high, is the candidate also
systematically high?

- **Data:** independent malignant epithelial cells, not Joanito or Pelka.
  Ideally the same independent CRC atlas used for Gate E.
- **Candidate statistic:** within-tumour (patient-blocked) correlation
  between the candidate and the CIOC score.
  - The CIOC score is computed leaving the candidate out, which matters
    only for the CIOC genes themselves.
  - The correlation is computed on malignant epithelial cells or
    metacells, and summarised across patients.
- **Controls:**
  - Use metacells or pseudobulk to limit dropout.
  - Regress or stratify by cell-cycle score so that proliferation does not
    drive coherence.
  - Build a null from random gene sets matched on expression level.
- **Rule:** to be frozen before computation. The final set size is not
  predetermined.

## Reference data (decided)

- **Primary:** Khaliq 2022 (GSE200997). 16 CRC + 7 adjacent normal, 10x,
  every compartment.
- **Replication:** Che 2021 (GSE178318). 6 patients, primary CRC + liver
  metastasis; PBMC excluded.
- **Independence:** neither dataset is used by any gate, and neither
  overlaps Joanito (SMC, KUL3, SG cohorts) or Pelka.
- **Secondary check only:** Tabula Sapiens Large Intestine (adult normal).
  It is reported as an annotation, never as a selection criterion.
- **Processing:**
  - `scripts/ext_01_prepare_atlases.py` builds raw-count AnnData objects.
  - `scripts/ext_02_annotate_compartments.py` annotates compartments.
    Neither deposit has cell types, so each Leiden cluster is assigned by
    canonical marker modules, scored as mean marker detection. No module
    contains a Level 2 candidate gene.
  - Khaliq: 1,236 ambiguous cells are excluded.
  - Tumour-tissue epithelium is treated as malignant. No CNV inference was
    run.

## Gate E — blinded descriptive profile (done)

Produced by `scripts/ext_03_gate_e_profile.py`.
- **Units:** patient × compartment pseudobulks from tumour-tissue cells.
  - Compartments: epithelial, stromal (fibroblast + SMC/pericyte),
    endothelial, myeloid, T/NK, B/plasma, mast.
  - A unit needs ≥ 20 cells, and a compartment needs ≥ 3 patients. Khaliq
    mast cells (2 patients) are dropped.
- **Ratio:** log2((epithelial CPM + 1) / (top non-epithelial compartment
  CPM + 1)), using patient medians.
- **Detectability:** the fraction of tumour epithelial cells with ≥ 1 UMI,
  as the patient median.
- **Blinding:** candidates are plotted unlabelled, and only quantiles and
  threshold counts were inspected.

**Calibrators (log2 ratio, Khaliq / Che):**

| Class | Genes | log2 ratio | Epithelial detection |
|---|---|---|---|
| Epithelial | EPCAM, KRT8, CDH1, KRT20, CEACAM5, CDX2 | +4.9 to +6.4 | 0.22–0.95 |
| Ubiquitous | GAPDH, ACTB, B2M | +0.3 to −2.5 | 0.73–0.99 |
| Non-epithelial, ambient-prone | LYZ, CD68 | −4.2 to −5.3 | 0.07–0.35 |
| Non-epithelial | PTPRC, CD3E, VIM, PECAM1, ACTA2, VWF, COL1A2, DCN | −5.1 to −11.8 | 0.00–0.30 |

Ubiquitous genes sit at −2.5 to +0.3. A max over six non-epithelial
compartments biases the ratio downward, and immune cells have different
library composition, so a truly shared gene is not expected at 0.

**Level 2 candidates (no gene identities inspected):**
- **Median log2 ratio:** −1.9 (Khaliq), −1.4 (Che). The background of
  expressed genes is centred near −0.5.
- **The candidates are shifted towards non-epithelial attribution.** For
  the candidates below 0, the top non-epithelial compartment is most often
  stromal, then endothelial.
- **Epithelial detectability is low in Khaliq.** Its epithelial cells have
  a median of 789 detected genes, against about 1,000 in other compartments.
  Che is more sensitive.

## Gate E — proposed rule v0.1 (for approval; not applied)

**E1 — non-epithelial attribution.**
- In an atlas, E1 fails if log2 ratio < −3.
- The threshold is calibrated between the ubiquitous calibrators (≥ −2.5)
  and the ambient-prone non-epithelial markers (≤ −4.2).
- The gene **fails E1 only if it fails in both atlases** (replicated
  non-epithelial attribution). If it is measured in only one atlas, that
  atlas decides.

**E2 — technical detectability.**
- The epithelial detection fraction must be ≥ 0.05 in at least one atlas.
- **Rationale:**
  - Ambient-free non-epithelial markers sit at ≤ 0.03 in epithelium.
  - In each atlas, 5% of tumour epithelial cells is several hundred cells,
    the minimum needed for the coherence step to have signal.
  - "At least one atlas" allows for the shallower Khaliq epithelium.

**Gate E pass = E2 AND NOT E1-fail.**
- Tabula Sapiens adult-normal attribution is reported beside each gene and
  never selects.

## CIOC program coherence — proposed rule v0.1 (for approval; not applied)

- **Cells:** tumour-tissue epithelial cells in Khaliq and Che. A patient
  needs ≥ 200 epithelial cells.
- **Units:** metacells within each patient: k-means on the patient's
  epithelial PCA, k = n_cells / 20, aggregated counts, log-CPM. A patient
  needs ≥ 10 metacells.
- **CIOC score:** the mean within-patient z-score of the eight CIOC genes.
  When the candidate is itself a CIOC gene, it is left out of the score.
- **Statistic:**
  - Within-patient Spearman ρ between the candidate and the CIOC score,
    after regressing out S and G2M cell-cycle scores.
  - Patients are combined by their Fisher-z mean per atlas.
- **Null:**
  - 1,000 random 8-gene sets, matched to the CIOC genes by expression bin.
  - The candidate's mean ρ with each random-set score gives an empirical
    P per gene and atlas.
  - BH is applied across Gate-E-passing candidates within each atlas.
- **Pass:** mean ρ > 0 in both atlases, and BH FDR < 0.05 in at least one.
  This mirrors Gate H: concordant in ≥ 2, supported in ≥ 1.

## Extended CIOC

**Extended CIOC = Level 2a ∩ Gate E ∩ coherence.**
- The set size is not predetermined.
- **Order of work:**
  1. Approve the rules.
  2. Freeze them (commit).
  3. Apply them once.
  4. Report without revision.
