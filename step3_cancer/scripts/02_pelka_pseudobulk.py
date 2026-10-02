#!/usr/bin/env python3
"""Pelka 2021 (GSE178341) epithelial raw-count pseudobulks.

Groups (frozen in ANALYSIS_PLAN.md):
  Tumor          epithelial cells (clTopLevel == "Epi") from SPECIMEN_TYPE T
  Normal         epithelial cells from SPECIMEN_TYPE N
  Normal_prolif  normal-specimen cells in cE01/cE02/cE03 (Stem/TA-like);
                 only used if the Joanito proliferation control is unavailable
One pseudobulk per patient (PID) x group; pseudobulks with < min_cells are
written to the inclusion table but flagged ineligible.
"""
import argparse
import pathlib

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

PROLIF = {"cE01", "cE02", "cE03"}


def read_10x_h5(path):
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        mat = sp.csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:])
        )
        barcodes = g["barcodes"][:].astype(str)
        feats = g["features"]
        ids = feats["id"][:].astype(str)
        names = feats["name"][:].astype(str)
    return mat, barcodes, ids, names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=50)
    a = ap.parse_args()
    raw = pathlib.Path(a.raw_dir)
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(raw / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
    clus = pd.read_csv(raw / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz").rename(
        columns={"sampleID": "cellID"}
    )
    cells = meta.merge(clus, on="cellID", validate="one_to_one")
    epi = cells[cells.clTopLevel == "Epi"].copy()
    epi["group"] = epi.SPECIMEN_TYPE.map({"T": "Tumor", "N": "Normal"})
    prolif = epi[(epi.SPECIMEN_TYPE == "N") & epi.cl295v11SubShort.isin(PROLIF)].copy()
    prolif["group"] = "Normal_prolif"
    assign = pd.concat([epi, prolif])
    assign["pseudobulk"] = assign.PID + "__" + assign.group

    mat, barcodes, ids, names = read_10x_h5(raw / "GSE178341_crc10x_full_c295v4_submit.h5")
    keep_feat = ~(pd.Series(names).str.match(r"^(HASH|ADT|TCR)-").to_numpy())
    col_index = pd.Series(np.arange(len(barcodes)), index=barcodes)
    missing = set(assign.cellID) - set(col_index.index)
    if missing:
        raise SystemExit(f"{len(missing)} annotated cells absent from the count matrix")

    pbs = sorted(assign.pseudobulk.unique())
    pb_codes = pd.Categorical(assign.pseudobulk, categories=pbs).codes
    cols = col_index.loc[assign.cellID].to_numpy()
    indicator = sp.csr_matrix(
        (np.ones(len(cols)), (cols, pb_codes)), shape=(mat.shape[1], len(pbs))
    )
    counts = (mat.tocsr()[keep_feat] @ indicator).toarray()
    if not np.allclose(counts, np.round(counts)):
        raise SystemExit("Matrix is not integer raw counts")

    df = pd.DataFrame(counts.astype(np.int64), index=ids[keep_feat], columns=pbs)
    df.insert(0, "symbol", names[keep_feat])
    df.index.name = "gene_id"
    df.to_csv(out / "Pelka_epithelial_pseudobulk_counts.csv.gz")

    samples = (
        assign.groupby("pseudobulk")
        .agg(patient=("PID", "first"), group=("group", "first"), n_cells=("cellID", "size"),
             MMRStatus=("MMRStatus", lambda s: s.dropna().iloc[0] if s.notna().any() else ""),
             chemistry=("SINGLECELL_TYPE", lambda s: ";".join(sorted(s.unique()))))
        .reset_index()
    )
    samples["library_size"] = df[samples.pseudobulk].sum().to_numpy()
    samples["eligible"] = samples.n_cells >= a.min_cells
    samples.to_csv(out / "Pelka_pseudobulk_samples.csv", index=False)
    print(samples.groupby(["group", "eligible"]).size())


if __name__ == "__main__":
    main()
