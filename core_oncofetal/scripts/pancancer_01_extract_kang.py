#!/usr/bin/env python3
"""Pan-cancer epithelial reactivation annotation, step 1: data extraction (no statistics).

Kang et al., Nat Commun 2024 tumour-normal pan-cancer atlas (Zenodo 10.5281/zenodo.10651059,
atlas_dataset.tar: one log1p(CP10k) .h5ad.xz per dataset; obs Dataset, Organ_origin, Sample,
Patient, Tissue, Cancer type, cnv_status, Celltype).
For each dataset: record the metadata vocabulary, and build pseudobulks per
Sample x Tissue x Cancer type x Organ_origin x Celltype x cnv_status = mean CP10k
(expm1 of log1p values) over cells, plus cell counts, for the annotated gene universe
(338 cross-species candidates + Literature-31 + calibrators), symbols resolved through HGNC.
Outputs (DATA work dir 2026-10-05_pancancer_epithelial_v0.1):
  Kang_metadata_vocabulary.csv, Kang_pseudobulk_meanCP10k.csv.gz, Kang_gene_resolution.csv
Usage (repo root): python3 core_oncofetal/scripts/pancancer_01_extract_kang.py
"""
import lzma
import pathlib
import shutil
import tarfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA")
TAR = DATA / "scRNAseq/Kang2024_pancancer_tumor_normal_atlas/raw/atlas_dataset.tar"
WORK = DATA / "2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-05_pancancer_epithelial_v0.1"
TMP = WORK / "tmp"
HGNC = DATA / "1.Databases/HGNC_gene_id_mapping/raw/hgnc_custom_download.tsv"
OUT = ROOT / "core_oncofetal/results"
CALIB = ["EPCAM", "KRT8", "CDH1", "MKI67", "TOP2A", "CEACAM5", "CDX2", "PTPRC", "COL1A2", "LYZ"]
KEYS = ["Dataset", "Sample", "Patient", "Tissue", "Cancer type", "Organ_origin", "Celltype", "cnv_status"]


def universe():
    l2 = pd.read_csv(OUT / "Genomewide_fetal_CRC_candidates_calls.csv")
    c338 = l2.loc[l2["Cross-species fetal–CRC candidate"] == "YES", "Gene"].tolist()
    lit = pd.read_csv(ROOT / "step2_fetal/config/literature_candidates_31.tsv", sep="\t").gene.tolist()
    return list(dict.fromkeys(c338 + lit + CALIB))


def resolver(var_names):
    h = pd.read_csv(HGNC, sep="\t", dtype=str, usecols=["Approved symbol", "Previous symbols"])
    prev = {a: [p.strip() for p in str(p).split(",") if p.strip() and p != "nan"]
            for a, p in zip(h["Approved symbol"], h["Previous symbols"])}
    cur = {p: a for a, ps in prev.items() for p in ps}  # previous symbol -> approved symbol
    vs = set(var_names)
    return lambda g: (g if g in vs else next((p for p in prev.get(g, []) if p in vs), None)
                      or (cur.get(g) if cur.get(g) in vs else None))


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(exist_ok=True)
    genes = universe()
    vocab, pbs, res_rows = [], [], []
    with tarfile.open(TAR) as tf:
        for m in tf:
            if not m.name.endswith(".h5ad.xz"):
                continue
            ds = pathlib.Path(m.name).name.replace(".h5ad.xz", "")
            h5 = TMP / f"{ds}.h5ad"
            with tf.extractfile(m) as src, lzma.open(src) as xz, open(h5, "wb") as dst:
                shutil.copyfileobj(xz, dst, length=16 << 20)
            a = ad.read_h5ad(h5, backed="r")
            res = resolver(a.var_names)
            sym = {g: res(g) for g in genes}
            res_rows += [dict(Dataset=ds, gene=g, dataset_symbol=s or "") for g, s in sym.items()]
            idx = [list(a.var_names).index(s) for s in sym.values() if s]
            have = [g for g, s in sym.items() if s]
            obs = a.obs[[k for k in KEYS if k in a.obs.columns]].copy()
            for k in KEYS:
                if k not in obs.columns:
                    obs[k] = "NA"
                obs[k] = obs[k].astype(str)
            for k in KEYS[3:]:
                for v, n in obs[k].value_counts().items():
                    vocab.append(dict(Dataset=ds, field=k, value=v, n_cells=n))
            X = a.X[:, idx] if not hasattr(a.X, "to_memory") else a.X[:, idx]
            X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
            X = X.expm1() if hasattr(X, "expm1") else sp.csr_matrix(np.expm1(X.toarray()))
            grp = obs.groupby(KEYS, observed=True).indices
            for key, rows in grp.items():
                mean = np.asarray(X[rows].mean(axis=0)).ravel()
                pbs.append(dict(zip(KEYS, key), n_cells=len(rows), **dict(zip(have, mean))))
            a.file.close()
            h5.unlink()
            print(ds, a.n_obs, "cells;", len(grp), "groups;", len(have), "genes", flush=True)
    pd.DataFrame(vocab).to_csv(WORK / "Kang_metadata_vocabulary.csv", index=False)
    pd.DataFrame(pbs).to_csv(WORK / "Kang_pseudobulk_meanCP10k.csv.gz", index=False)
    pd.DataFrame(res_rows).to_csv(WORK / "Kang_gene_resolution.csv", index=False)
    TMP.rmdir()


if __name__ == "__main__":
    main()
