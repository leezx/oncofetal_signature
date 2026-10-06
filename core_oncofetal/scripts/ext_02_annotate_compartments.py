#!/usr/bin/env python3
"""Extended CIOC, step 2: compartment annotation of the independent CRC
atlases (Khaliq 2022, Che 2021). Neither GEO deposit provides cell types.

Per dataset: tissue cells only (Che PBMC excluded); QC (>= 200 genes,
>= 500 UMI, mito < 25%); log-normalise; 2,000 HVGs; PCA 30; kNN; Leiden 1.0.
Each cluster is assigned the compartment whose canonical marker module has
the highest module detection score = mean, over the module's marker genes, of
the fraction of cluster cells with the marker detected (modules contain no
Level 2 candidate genes; checked). Rules:
  - Hepatocyte is allowed only for clusters with >= 50% liver-metastasis cells
    (colon tumour cells can express APOA1/HP/SERPINA1);
  - Fibroblast and SMC_pericyte as the top two -> 'Stromal' (both non-epithelial
    stromal; kept separate otherwise);
  - 'Ambiguous' (excluded downstream) if the best score < 0.3 or the runner-up
    is >= 70% of it.
(score_genes was tried first and under-scored abundant T cells; replaced.)
Outputs (DATA work dir 2026-10-02_extended_cioc_v0.1):
  <ds>_cell_compartments.csv.gz, <ds>_cluster_summary.csv, <ds>_umap_compartments.png
Usage (repo root): python3 core_oncofetal/scripts/ext_02_annotate_compartments.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc

SC = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq")
WORK = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_extended_cioc_v0.1")
MARKERS = {
    "Epithelial": "EPCAM KRT8 KRT18 KRT19 CDH1 KRT20 CDH17 ELF3",
    "Fibroblast": "COL1A2 COL3A1 DCN LUM PDGFRA FBLN1",
    "SMC_pericyte": "ACTA2 DES MYH11 RGS5 PDGFRB CNN1",
    "Endothelial": "PECAM1 VWF CDH5 CLDN5 KDR PLVAP",
    "Myeloid": "CD68 LYZ CD14 C1QA CSF1R FCGR3A",
    "T_NK": "CD3E CD3D CD2 NKG7 GZMA CD7",
    "B_plasma": "MS4A1 CD79A CD19 JCHAIN MZB1 IGHG1",
    "Mast": "TPSAB1 CPA3 KIT MS4A2",
    "Hepatocyte": "ALB APOA1 APOC3 TTR HP SERPINA1",
}
DATASETS = {"Khaliq2022": SC / "GSE200997_Khaliq2022/processed/v0.1/Khaliq2022_raw.h5ad",
            "Che2021": SC / "GSE178318_Che2021/processed/v0.1/Che2021_raw.h5ad"}


def annotate(name, path):
    a = sc.read_h5ad(path)
    a = a[a.obs.tissue != "PBMC"].copy()
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    n0 = a.n_obs
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.total_counts >= 500) & (a.obs.pct_counts_mt < 25)].copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="samples")
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30)
    sc.pp.neighbors(b, n_neighbors=15)
    sc.tl.leiden(b, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(b)
    a.obs["leiden"] = b.obs.leiden.values
    a.obsm["X_umap"] = b.obsm["X_umap"]
    det = {}
    for k, v in MARKERS.items():
        genes = [g for g in v.split() if g in a.var_names]
        x = (a[:, genes].layers["counts"] > 0).astype(np.float32)
        a.obs[f"det_{k}"] = np.asarray(x.mean(axis=1)).ravel()
    cm = a.obs.groupby("leiden", observed=True)[[f"det_{k}" for k in MARKERS]].mean()
    cm = cm.rename(columns=lambda c: c.replace("det_", ""))
    lm = a.obs.groupby("leiden", observed=True).tissue.apply(lambda s: (s == "LM").mean())
    cm.loc[lm.reindex(cm.index) < 0.5, "Hepatocyte"] = -np.inf
    lab = []
    for cl, r in cm.iterrows():
        o = r.sort_values(ascending=False)
        (b1, s1), (b2, s2) = list(o.items())[:2]
        if {b1, b2} == {"Fibroblast", "SMC_pericyte"} and s1 >= 0.3:
            lab.append("Stromal")
        elif s1 < 0.3 or s2 >= 0.7 * s1:
            lab.append("Ambiguous")
        else:
            lab.append(b1)
    lab = np.array(lab)
    summ = cm.replace(-np.inf, np.nan).round(2).rename(columns=lambda c: "det_" + c)
    summ.insert(0, "compartment", lab)
    summ.insert(1, "n_cells", a.obs.leiden.value_counts().reindex(summ.index).values)
    summ.insert(2, "frac_tumour_tissue", a.obs.groupby("leiden", observed=True).tissue.apply(
        lambda s: (s != "Normal").mean()).round(2).values)
    summ.insert(3, "frac_LM", lm.reindex(summ.index).round(2).values)
    summ.insert(4, "n_patients", a.obs.groupby("leiden", observed=True).patient.nunique().values)
    a.obs["compartment"] = a.obs.leiden.map(dict(zip(summ.index, lab))).astype(str)
    WORK.mkdir(parents=True, exist_ok=True)
    summ.to_csv(WORK / f"{name}_cluster_summary.csv")
    a.obs[["samples", "patient", "tissue", "leiden", "compartment", "n_genes_by_counts", "total_counts",
           "pct_counts_mt"]].to_csv(WORK / f"{name}_cell_compartments.csv.gz")
    sc.settings.figdir = WORK
    sc.pl.umap(a, color=["compartment", "tissue", "patient"], save=f"_{name}_compartments.png", show=False, ncols=3)
    print(f"{name}: {n0} tissue cells -> {a.n_obs} after QC; {a.obs.leiden.nunique()} clusters")
    print(a.obs.groupby(["compartment", "tissue"], observed=True).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    for n, p in DATASETS.items():
        annotate(n, p)
