#!/usr/bin/env python3
"""H-new2 fetal arm: Fawkner-Corbett 2021 EPCAM+ pools (GSE158702), hashtag-
demultiplexed into fetal samples; one raw-count pseudobulk per pool x hashtag.

Cell QC: >= 500 UMIs and >= 200 detected genes. Hashtag call: the highest HTO
count is >= 10 and >= 2x the second highest; otherwise the cell is dropped
(negative/doublet). Rules frozen in docs/QUALIFICATION_PLAN.md.
"""
import argparse
import gzip
import pathlib
import tarfile

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

# GEX archive -> HTO archive from the same 10x reaction. The GEO sample titles
# do not pair the libraries correctly; pairing was established empirically as
# the HTO library containing >= 99% of the pool's QC-passing cell barcodes
# (EPI2 <-> HTO3: 99.3%, EPI3 <-> HTO5: 99.7%, EPI pool 4 <-> "HTO_stromal_4":
# 99.1%; every other HTO library < 13%). Technical pairing only; no genes used.
POOLS = {
    "GSM4808339_EPI1_RUN3": "GSM4808349_HTO1",
    "GSM4808340_EPI2_RUN3": "GSM4808351_HTO3",
    "GSM4808341_EPI3_RUN3": "GSM4808353_HTO5",
    "GSM4808345_EPI_run2": "GSM4808356_HTO_stromal_4",
}


def read_mtx(archive, folder):
    t = tarfile.open(archive, "r:gz")
    get = lambda name: t.extractfile(next(m for m in t.getmembers() if m.name.endswith(f"{folder}/{name}")))
    mat = sp.csc_matrix(scipy.io.mmread(gzip.open(get("matrix.mtx.gz"))))
    bcs = pd.read_csv(gzip.open(get("barcodes.tsv.gz")), header=None)[0].astype(str)
    feats = pd.read_csv(gzip.open(get("features.tsv.gz")), header=None, sep="\t")
    return mat, bcs.str.replace(r"-1$", "", regex=True).to_numpy(), feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=50)
    a = ap.parse_args()
    raw, work = pathlib.Path(a.raw_dir), pathlib.Path(a.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    pbs, rows, genes = {}, [], None
    for gex_name, hto_name in POOLS.items():
        pool = gex_name.split("_", 1)[1]
        gex, gbc, gfeat = read_mtx(raw / f"{gex_name}.tar.gz", "raw_feature_bc_matrix")
        if genes is None:
            genes = gfeat
        elif not (genes[0].to_numpy() == gfeat[0].to_numpy()).all():
            raise SystemExit("Feature order differs between pools")
        umis = np.asarray(gex.sum(axis=0)).ravel()
        ngenes = np.asarray((gex > 0).sum(axis=0)).ravel()
        qc = (umis >= 500) & (ngenes >= 200)
        hto, hbc, hfeat = read_mtx(raw / f"{hto_name}.tar.gz", "umi_count")
        tags = hfeat[0].astype(str).to_numpy()
        # CITE-seq-Count appends an "unmapped" row that is not in features.tsv.
        hto = hto[: len(tags)].toarray()
        top_idx = hto.argmax(axis=0)
        srt = np.sort(hto, axis=0)
        top, second = srt[-1], srt[-2] if hto.shape[0] > 1 else np.zeros(hto.shape[1])
        called = (top >= 10) & (top >= 2 * second)
        call = pd.Series(np.where(called, tags[top_idx], ""), index=hbc)
        cell_tag = call.reindex(gbc).fillna("").to_numpy()
        use = qc & (cell_tag != "")
        n_qc, n_tagged = int(qc.sum()), int(use.sum())
        for tag in sorted(set(cell_tag[use])):
            sel = use & (cell_tag == tag)
            unit = f"{pool}__{tag.split('-')[0]}"
            pbs[unit] = np.asarray(gex[:, sel].sum(axis=1)).ravel().astype(np.int64)
            rows.append(dict(pseudobulk=unit, patient=unit, group="Fetal", pool=pool, hashtag=tag,
                             n_cells=int(sel.sum()), library_size=int(pbs[unit].sum()),
                             pool_cells_passing_qc=n_qc, pool_cells_hashtag_called=n_tagged))
    counts = pd.DataFrame(pbs)
    counts.insert(0, "symbol", genes[1].to_numpy())
    counts.index = genes[0].to_numpy()
    counts.index.name = "gene_id"
    counts.to_csv(work / "Fawkner_scRNA_fetal_pseudobulk_counts.csv.gz")
    units = pd.DataFrame(rows)
    units["eligible"] = units.n_cells >= a.min_cells
    units.to_csv(work / "Fawkner_scRNA_fetal_units.csv", index=False)
    print(units[["pseudobulk", "n_cells", "library_size", "eligible"]].to_string(index=False))
    print(units.groupby("pool")[["pool_cells_passing_qc", "pool_cells_hashtag_called"]].first())


if __name__ == "__main__":
    main()
