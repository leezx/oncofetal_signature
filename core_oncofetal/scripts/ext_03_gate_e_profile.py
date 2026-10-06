#!/usr/bin/env python3
"""Extended CIOC, step 3: Gate E descriptive profile (BLINDED; no rule applied).

For each independent CRC atlas (Khaliq 2022, Che 2021), patient x compartment
pseudobulks are built from tumour-tissue cells (Khaliq tumour samples; Che
primary CRC + liver metastasis), excluding Ambiguous clusters:
  Epithelial (tumour-derived epithelial cells; malignancy not inferred),
  Stromal, Endothelial, Myeloid, T_NK, B_plasma, Mast.
A patient contributes to a compartment only with >= 20 cells; a compartment is
used only with >= 3 such patients (Khaliq mast: 2 patients, excluded).
Per gene and dataset:
  epi_CPM             median over patients of epithelial pseudobulk CPM
  top_nonepi_CPM      max over non-epithelial compartments of the patient-median CPM
  top_nonepi          that compartment
  log2_epi_vs_nonepi  log2((epi_CPM + 1) / (top_nonepi_CPM + 1))
  epi_detect          median over patients of the fraction of epithelial cells
                      with >= 1 UMI
The plot shows the distribution for the Level 2 cross-species candidates
WITHOUT gene labels, against all expressed genes, with canonical lineage
markers as calibrators. The per-gene table is written for later rule
application only; it is restricted (Level 2 membership) and is not used to
choose thresholds.
Outputs (DATA work dir 2026-10-02_extended_cioc_v0.1):
  GateE_gene_profile.csv.gz (restricted), GateE_calibrators.csv,
  GateE_distribution.png, GateE_compartment_counts.csv
Usage (repo root): python3 core_oncofetal/scripts/ext_03_gate_e_profile.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]
SC = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq")
WORK = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_extended_cioc_v0.1")
RAW = {"Khaliq2022": SC / "GSE200997_Khaliq2022/processed/v0.1/Khaliq2022_raw.h5ad",
       "Che2021": SC / "GSE178318_Che2021/processed/v0.1/Che2021_raw.h5ad"}
COMPS = ["Epithelial", "Stromal", "Endothelial", "Myeloid", "T_NK", "B_plasma", "Mast"]
MIN_CELLS = 20
MIN_PATIENTS = 3
CALIB = {"epithelial": ["EPCAM", "KRT8", "CDH1", "KRT20", "CEACAM5", "CDX2"],
         "non-epithelial": ["COL1A2", "DCN", "ACTA2", "PECAM1", "VWF", "PTPRC", "CD3E", "LYZ", "CD68"],
         "shared/ubiquitous": ["ACTB", "B2M", "GAPDH", "VIM"]}


def profile(name):
    a = sc.read_h5ad(RAW[name])
    o = pd.read_csv(WORK / f"{name}_cell_compartments.csv.gz", index_col=0)
    o = o[(o.compartment != "Ambiguous") & (o.tissue != "Normal")]
    a = a[o.index]
    a.obs["comp"] = o.compartment.values
    a.obs["unit"] = a.obs.patient.astype(str) + "|" + a.obs.comp.astype(str)
    units = a.obs.unit.value_counts()
    keep = units[units >= MIN_CELLS].index
    a = a[a.obs.unit.isin(keep)]
    u = pd.Categorical(a.obs.unit)
    ind = sp.csr_matrix((np.ones(a.n_obs), (u.codes, np.arange(a.n_obs))), shape=(len(u.categories), a.n_obs))
    X = a.X.tocsr()
    pb = np.asarray((ind @ X).todense())
    cpm = pb / pb.sum(axis=1, keepdims=True) * 1e6
    det = np.asarray((ind @ (X > 0).astype(np.float32)).todense()) / np.asarray(ind.sum(axis=1))
    meta = pd.DataFrame([c.split("|") for c in u.categories], columns=["patient", "comp"])
    counts = meta.assign(n_cells=units.reindex(u.categories).values)
    med = {c: np.median(cpm[(meta.comp == c).values], axis=0) for c in COMPS
           if (meta.comp == c).sum() >= MIN_PATIENTS}
    epi_det = np.median(det[(meta.comp == "Epithelial").values], axis=0)
    non = pd.DataFrame({c: v for c, v in med.items() if c != "Epithelial"}, index=a.var_names)
    df = pd.DataFrame({"epi_CPM": med["Epithelial"], "top_nonepi_CPM": non.max(axis=1).values,
                       "top_nonepi": non.idxmax(axis=1).values, "epi_detect": epi_det}, index=a.var_names)
    df["log2_epi_vs_nonepi"] = np.log2((df.epi_CPM + 1) / (df.top_nonepi_CPM + 1))
    for c, v in med.items():
        df[f"CPM_{c}"] = v
    return df, counts.assign(dataset=name)


def main():
    calls = pd.read_csv(ROOT / "core_oncofetal/results/Genomewide_fetal_CRC_candidates_calls.csv")
    l2 = set(calls.loc[calls["Cross-species fetal–CRC candidate"] == "YES", "Gene"])
    prof, cnt = {}, []
    for n in RAW:
        prof[n], c = profile(n)
        cnt.append(c)
    pd.concat(cnt).to_csv(WORK / "GateE_compartment_counts.csv", index=False)
    long = pd.concat({n: p for n, p in prof.items()}, names=["dataset", "gene"]).reset_index()
    long["level2_candidate"] = long.gene.isin(l2)
    long.to_csv(WORK / "GateE_gene_profile.csv.gz", index=False)
    cal = long[long.gene.isin(sum(CALIB.values(), []))].copy()
    cal["calibrator_class"] = cal.gene.map({g: k for k, v in CALIB.items() for g in v})
    cal[["dataset", "gene", "calibrator_class", "epi_CPM", "top_nonepi_CPM", "top_nonepi", "log2_epi_vs_nonepi",
         "epi_detect"]].round(3).to_csv(WORK / "GateE_calibrators.csv", index=False)

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    col = {"epithelial": "tab:blue", "non-epithelial": "tab:red", "shared/ubiquitous": "tab:green"}
    for i, n in enumerate(RAW):
        p = prof[n]
        bg = p[(p.epi_CPM >= 1) | (p.top_nonepi_CPM >= 1)]
        cand = p[p.index.isin(l2)]
        bins = np.linspace(-10, 10, 81)
        a0 = ax[i, 0]
        a0.hist(bg.log2_epi_vs_nonepi, bins=bins, density=True, color="lightgrey", label=f"expressed genes (n={len(bg)})")
        a0.hist(cand.log2_epi_vs_nonepi, bins=bins, density=True, histtype="step", lw=2, color="purple",
                label=f"Level 2 candidates (n={len(cand)}, unlabelled)")
        for k, gs in CALIB.items():
            for g in gs:
                if g in p.index:
                    x = p.loc[g, "log2_epi_vs_nonepi"]
                    a0.axvline(x, color=col[k], lw=0.8, alpha=0.7)
                    a0.text(x, a0.get_ylim()[1] * 0.95, g, rotation=90, fontsize=6, color=col[k], va="top")
        a0.set_xlabel("log2((epithelial CPM+1)/(top non-epithelial compartment CPM+1))")
        a0.set_title(f"{n}: epithelial vs top non-epithelial compartment")
        a0.legend(fontsize=7, loc="upper left")
        a1 = ax[i, 1]
        b2 = np.linspace(0, 1, 51)
        a1.hist(bg.epi_detect, bins=b2, density=True, color="lightgrey", label="expressed genes")
        a1.hist(cand.epi_detect, bins=b2, density=True, histtype="step", lw=2, color="purple", label="Level 2 candidates")
        for k, gs in CALIB.items():
            for g in gs:
                if g in p.index:
                    a1.axvline(p.loc[g, "epi_detect"], color=col[k], lw=0.8, alpha=0.7)
        a1.set_xlabel("fraction of tumour epithelial cells detecting the gene (patient median)")
        a1.set_title(f"{n}: epithelial detectability")
        a1.legend(fontsize=7)
        a2 = ax[i, 2]
        a2.scatter(np.log2(bg.epi_CPM + 1), np.log2(bg.top_nonepi_CPM + 1), s=2, c="lightgrey", rasterized=True)
        a2.scatter(np.log2(cand.epi_CPM + 1), np.log2(cand.top_nonepi_CPM + 1), s=6, c="purple", alpha=0.6)
        for k, gs in CALIB.items():
            for g in gs:
                if g in p.index:
                    x, y = np.log2(p.loc[g, "epi_CPM"] + 1), np.log2(p.loc[g, "top_nonepi_CPM"] + 1)
                    a2.scatter(x, y, s=18, c=col[k])
                    a2.text(x, y, g, fontsize=6, color=col[k])
        lim = [0, max(a2.get_xlim()[1], a2.get_ylim()[1])]
        a2.plot(lim, lim, "k--", lw=0.6)
        a2.set_xlabel("log2(epithelial CPM + 1)")
        a2.set_ylabel("log2(top non-epithelial CPM + 1)")
        a2.set_title(f"{n}: candidates unlabelled; calibrators labelled")
    fig.tight_layout()
    fig.savefig(WORK / "GateE_distribution.png", dpi=150)

    # blinded summary: distribution quantiles only
    rows = []
    for n, p in prof.items():
        cand = p[p.index.isin(l2)]
        q = cand.log2_epi_vs_nonepi.quantile([0.05, 0.1, 0.25, 0.5, 0.75]).round(2).to_dict()
        rows.append(dict(dataset=n, n_candidates=len(cand),
                         **{f"log2ratio_q{int(k*100)}": v for k, v in q.items()},
                         **{f"n_log2ratio_below_{t}": int((cand.log2_epi_vs_nonepi < t).sum()) for t in (-1, -2, -3, -4)},
                         **{f"n_epi_detect_below_{t}": int((cand.epi_detect < t).sum()) for t in (0.02, 0.05, 0.1)},
                         top_nonepi_counts=cand[cand.log2_epi_vs_nonepi < 0].top_nonepi.value_counts().to_dict()))
    s = pd.DataFrame(rows)
    s.to_csv(WORK / "GateE_blinded_distribution_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(s.T.to_string())
    print(cal.pivot_table(index=["calibrator_class", "gene"], columns="dataset",
                          values=["log2_epi_vs_nonepi", "epi_detect"]).round(2).to_string())


if __name__ == "__main__":
    main()
