#!/usr/bin/env python3
"""Extended CIOC, step 1: build raw-count AnnData objects for the two
independent all-compartment CRC atlases (Gate E and CIOC coherence).

  Khaliq 2022 (GSE200997): 16 CRC + 7 adjacent normal, 10x; GEO gives a dense
    gene x cell UMI CSV and sample/condition/location/MSI/CMS annotation (no
    cell types).
  Che 2021 (GSE178318): CRC primary, liver metastasis and PBMC, 10x MTX;
    barcodes carry <barcode>_<patient>_<tissue>; no cell annotation.
Output: DATA/scRNAseq/<dataset>/processed/v0.1/<name>_raw.h5ad
Usage (repo root): python3 core_oncofetal/scripts/ext_01_prepare_atlases.py
"""
import gzip
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

SC = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/scRNAseq")


def khaliq():
    d = SC / "GSE200997_Khaliq2022"
    ann = pd.read_csv(d / "raw/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz", index_col=0)
    genes, blocks = [], []
    reader = pd.read_csv(d / "raw/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz", index_col=0,
                         chunksize=2000)
    cells = None
    for ch in reader:
        cells = ch.columns if cells is None else cells
        genes.extend(ch.index.astype(str))
        blocks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32)))
    X = sp.vstack(blocks).T.tocsr()
    obs = ann.reindex(cells)
    obs["in_geo_annotation"] = obs.index.isin(ann.index)
    obs["samples"] = ["_".join(c.split("_")[:2]) for c in cells]  # B_/T_<patient>_<barcode>
    obs["patient"] = obs.samples.str.split("_").str[1]
    obs["tissue"] = obs.samples.str[0].map({"T": "CRC", "B": "Normal"})
    obs = obs.astype({c: str for c in ("Condition", "Location", "MSI_Status", "bulk_prediction", "prediction")})
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="symbol")))
    a.var_names_make_unique()
    a.obs["dataset"] = "Khaliq2022"
    out = d / "processed/v0.1"
    out.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out / "Khaliq2022_raw.h5ad", compression="gzip")
    print("Khaliq", a.shape, a.obs.tissue.value_counts().to_dict(), a.obs.patient.nunique(), "patients")


def che():
    d = SC / "GSE178318_Che2021"
    X = scipy.io.mmread(gzip.open(d / "raw/GSE178318_matrix.mtx.gz")).T.tocsr().astype(np.float32)
    bc = pd.read_csv(d / "raw/GSE178318_barcodes.tsv.gz", header=None)[0]
    g = pd.read_csv(d / "raw/GSE178318_genes.tsv.gz", sep="\t", header=None)
    parts = bc.str.split("_")
    obs = pd.DataFrame({"patient": parts.str[-2].values, "tissue": parts.str[-1].values}, index=bc.values)
    obs["samples"] = obs.patient + "_" + obs.tissue
    var = pd.DataFrame({"gene_id": g[0].values}, index=pd.Index(g[1].astype(str).values, name="symbol"))
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.obs["dataset"] = "Che2021"
    out = d / "processed/v0.1"
    out.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out / "Che2021_raw.h5ad", compression="gzip")
    print("Che", a.shape, a.obs.tissue.value_counts().to_dict(), a.obs.patient.nunique(), "patients")


if __name__ == "__main__":
    khaliq()
    che()
