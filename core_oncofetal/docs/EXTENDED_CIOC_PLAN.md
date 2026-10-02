# Extended CIOC — plan (draft; rules NOT yet frozen)

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

## Next decisions (owner: user)

1. Which independent CRC atlas (all compartments, non-overlapping patients)
   to download for Gate E and coherence.
2. Whether a fetal all-compartment reference is required for Gate E, or
   whether malignant plus Tabula Sapiens is sufficient.
