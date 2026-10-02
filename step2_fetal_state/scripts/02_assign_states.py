#!/usr/bin/env python3
"""Epithelial state labels (plan v1), identical procedure for Fawkner and Gao.

Reference = the Fawkner authors' sub-cluster marker lists (Mendeley sheets 3
and 10), top 30 per author cluster (p_val_adj < 0.05), pooled into 9 states
(config/state_map.tsv). All analysed genes and aliases (config/genes.tsv) are
removed from the marker sets. Each Leiden cluster gets the state with the
highest cluster-mean score (mean per-gene z of the state's markers). Labels are
written before any analysed gene is tabulated; the label check uses only the
control genes. Addendum v1.2: a cluster keeps its state only if the winning
score is > 0.2 and leads the runner-up by >= 0.1; otherwise "Unresolved".
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

CONTROLS = ["EPCAM", "OLFM4", "LGR5", "ASCL2", "MKI67"]
RES = {"Fawkner2021": 1.0, "Gao2018": 1.5}
MIN_SCORE, MIN_MARGIN = 0.2, 0.1


def marker_sets(xlsx, state_map, excluded, top_n=30):
    sets = {}
    for sheet, rows in state_map.groupby("source_sheet"):
        m = pd.read_excel(xlsx, sheet_name=sheet)
        m = m[m.p_val_adj < 0.05]
        for _, r in rows.iterrows():
            top = m[m.cluster == r.author_cluster].sort_values("avg_logFC", ascending=False).gene.head(top_n)
            if top.empty:
                raise SystemExit(f"No markers for {r.author_cluster}")
            sets.setdefault(r.state, set()).update(top)
    return {s: sorted(g - excluded) for s, g in sets.items()}


def label(a, name, sets, out):
    a = a.copy()
    if name == "Fawkner2021":
        sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers["lognorm"] = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, random_state=0)
    sc.pp.neighbors(b, n_neighbors=15, random_state=0)
    sc.tl.leiden(b, resolution=RES[name], random_state=0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(b, random_state=0)
    a.obs["leiden"] = b.obs.leiden.astype(str).to_numpy()
    a.obsm["X_umap"] = b.obsm["X_umap"]
    # State scores on per-gene z (all cells of the dataset).
    used = {s: [g for g in gs if g in a.var_names] for s, gs in sets.items()}
    genes = sorted(set().union(*used.values()))
    z = a[:, genes].X
    z = np.asarray(z.todense()) if hasattr(z, "todense") else np.asarray(z)
    z = (z - z.mean(axis=0)) / np.where(z.std(axis=0) > 0, z.std(axis=0), 1)
    zdf = pd.DataFrame(z, columns=genes, index=a.obs_names)
    scores = pd.DataFrame({s: zdf[g].mean(axis=1) for s, g in used.items() if g})
    cl = scores.groupby(a.obs.leiden).mean()
    cl["state"] = cl.idxmax(axis=1)
    srt = np.sort(cl.drop(columns="state").to_numpy(), axis=1)
    cl["margin"] = srt[:, -1] - srt[:, -2]
    cl["best_score"] = srt[:, -1]
    cl["n_cells"] = a.obs.leiden.value_counts().reindex(cl.index)
    # Addendum v1.2: low-confidence clusters are Unresolved.
    cl["forced_state"] = cl.state
    cl.loc[(cl.best_score <= MIN_SCORE) | (cl.margin < MIN_MARGIN), "state"] = "Unresolved"
    a.obs["state"] = a.obs.leiden.map(cl.state)
    cl.round(3).to_csv(out / f"{name}_cluster_state_scores.csv")
    pd.DataFrame([dict(state=s, n_markers_used=len(g), markers=";".join(g)) for s, g in used.items()]).to_csv(
        out / f"{name}_state_marker_sets_used.csv", index=False)
    # Label check: control genes by state (mean log-normalised, % detected).
    rows = []
    for st, idx in a.obs.groupby("state").indices.items():
        for g in CONTROLS:
            if g in a.var_names:
                v = np.asarray(a.layers["lognorm"][idx][:, a.var_names.get_loc(g)].todense()).ravel() \
                    if hasattr(a.layers["lognorm"], "todense") else a.layers["lognorm"][idx, a.var_names.get_loc(g)]
                rows.append(dict(state=st, gene=g, n_cells=len(idx), mean_lognorm=v.mean(), pct_detected=100 * (v > 0).mean()))
    pd.DataFrame(rows).round(3).to_csv(out / f"{name}_state_label_check_controls.csv", index=False)
    return a.obs[["leiden", "state"]], a.obsm["X_umap"], cl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--markers-xlsx", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    a = ap.parse_args()
    repo, work, tables = pathlib.Path(a.repo), pathlib.Path(a.work_dir), pathlib.Path(a.tables_dir)
    cfg = repo / "step2_fetal_state/config"
    genes = pd.read_csv(cfg / "genes.tsv", sep="\t", keep_default_na=False)
    excluded = set(genes.symbol) | set(genes.alias[genes.alias != ""])
    sets = marker_sets(a.markers_xlsx, pd.read_csv(cfg / "state_map.tsv", sep="\t"), excluded)
    pd.set_option("display.width", 220)
    for name, f in (("Fawkner2021", "fawkner_epithelial_cells.h5ad"), ("Gao2018", "gao_epithelial_cells.h5ad")):
        cells = ad.read_h5ad(work / f)
        labels, umap, cl = label(cells, name, sets, tables)
        labels["umap1"], labels["umap2"] = umap[:, 0], umap[:, 1]
        labels.to_csv(work / f"{name}_state_labels.csv.gz")
        print(name, cells.shape)
        print(cl[["forced_state", "state", "best_score", "margin", "n_cells"]].sort_values("state").to_string())
        print(labels.state.value_counts().to_string())


if __name__ == "__main__":
    main()
