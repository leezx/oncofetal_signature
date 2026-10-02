#!/usr/bin/env python3
"""Joanito 2022 (Synapse syn26844071) epithelial raw-count pseudobulks.

Cells are assigned to analysis groups only through the frozen
config/joanito_label_map.tsv (sample.origin x iCMS -> Malignant / Normal /
EXCLUDED). One pseudobulk per patient x group; Tumor and Tumor-2 samples of a
patient are pooled. `cohort` (the authors' dataset) is kept as a covariate.
"""
import argparse
import pathlib

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_h5(path):
    with h5py.File(path, "r") as f:
        keys = []
        f.visit(keys.append)
        grp = next(k for k in ("matrix", "X") if k in f)
        g = f[grp]
        if "features" in g:  # 10x v3 layout
            names = g["features"]["name"][:].astype(str)
            # The released file has an empty `id` group (symbols only,
            # GRCh38_ensembl93); fall back to symbols as identifiers.
            fid = g["features"]["id"]
            ids = fid[:].astype(str) if isinstance(fid, h5py.Dataset) else names.copy()
        else:  # 10x v2 layout
            names = g["gene_names"][:].astype(str)
            ids = g["genes"][:].astype(str)
        mat = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        barcodes = g["barcodes"][:].astype(str)
    return mat, barcodes, ids, names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--label-map", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=50)
    a = ap.parse_args()
    raw = pathlib.Path(a.raw_dir)
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(raw / "Epithelial_metadata.csv", index_col=0)
    lmap = pd.read_csv(a.label_map, sep="\t")
    meta = meta.reset_index().merge(
        lmap[["sample_origin", "iCMS", "analysis_group"]],
        left_on=["sample.origin", "iCMS"], right_on=["sample_origin", "iCMS"], how="left",
        validate="many_to_one")
    if meta.analysis_group.isna().any():
        raise SystemExit("Unmapped sample.origin x iCMS combinations; update the frozen label map first")
    cell_col = meta.columns[0]
    assign = meta[meta.analysis_group != "EXCLUDED"].copy()
    assign["pseudobulk"] = assign["patient.ID"] + "__" + assign.analysis_group

    mat, barcodes, ids, names = read_h5(raw / "Epithelial_Count_matrix.h5")
    if mat.shape[1] != len(barcodes):  # stored genes x cells vs cells x genes
        mat = mat.T.tocsc()
    if mat.shape[0] != len(names):
        mat = mat.T.tocsc()
    col_index = pd.Series(np.arange(len(barcodes)), index=barcodes)
    missing = set(assign[cell_col]) - set(col_index.index)
    if missing:
        raise SystemExit(f"{len(missing)} metadata cells absent from the count matrix")

    pbs = sorted(assign.pseudobulk.unique())
    codes = pd.Categorical(assign.pseudobulk, categories=pbs).codes
    cols = col_index.loc[assign[cell_col]].to_numpy()
    ind = sp.csr_matrix((np.ones(len(cols)), (cols, codes)), shape=(mat.shape[1], len(pbs)))
    counts = (mat.tocsr() @ ind).toarray()
    if not np.allclose(counts, np.round(counts)):
        raise SystemExit("Matrix is not integer raw counts")

    df = pd.DataFrame(counts.astype(np.int64), index=ids, columns=pbs)
    dup = df.index.duplicated(keep=False)
    if dup.any():  # duplicated symbols cannot be matched unambiguously
        print(f"Dropping {dup.sum()} rows with duplicated identifiers")
        df = df[~dup]
        names = names[~dup]
    df.insert(0, "symbol", names)
    df.index.name = "gene_id"
    df.to_csv(out / "Joanito_epithelial_pseudobulk_counts.csv.gz")

    samples = (assign.groupby("pseudobulk")
               .agg(patient=("patient.ID", "first"), group=("analysis_group", "first"),
                    cohort=("dataset", "first"), n_cells=(cell_col, "size"),
                    n_samples=("sample.ID", "nunique"), msi=("msi", "first"))
               .reset_index())
    # Subtype composition per malignant pseudobulk: descriptive only, not in
    # any model; lets a later marker failure be traced to iCMS composition.
    comp = (assign[assign.analysis_group == "Malignant"]
            .groupby("pseudobulk").iCMS.value_counts(normalize=True).unstack(fill_value=0) * 100)
    samples["pct_iCMS2"] = samples.pseudobulk.map(comp.get("iCMS2", pd.Series(dtype=float))).round(1)
    samples["pct_iCMS3"] = samples.pseudobulk.map(comp.get("iCMS3", pd.Series(dtype=float))).round(1)
    samples["library_size"] = df[samples.pseudobulk].sum().to_numpy()
    samples["eligible"] = samples.n_cells >= a.min_cells
    samples.to_csv(out / "Joanito_pseudobulk_samples.csv", index=False)
    print(samples.groupby(["group", "eligible"]).size())


if __name__ == "__main__":
    main()
