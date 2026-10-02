#!/usr/bin/env python3
"""H-new2 adult arm: Burclaff 2022 (GSE185224) healthy adult epithelium.

All cells of the authors' annotated epithelial object, one pseudobulk per
donor (all six regions pooled). Counts are taken from the per-donor full-gene
Cell Ranger matrices (Gene Expression features) for exactly those annotated
cells: the published h5ad keeps only 23,170 genes (e.g. SPP1, ANKRD1, SOX17
were removed by the authors' gene filter), which would bias marker coverage. The adult pseudobulks are joined
to the Fawkner fetal pseudobulks by Ensembl gene ID into one count matrix for
edgeR. Rules frozen in docs/QUALIFICATION_PLAN.md.
"""
import argparse
import pathlib

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--fetal-counts", required=True)
    ap.add_argument("--fetal-units", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=50)
    a = ap.parse_args()
    work = pathlib.Path(a.work_dir)
    obs = ad.read_h5ad(a.h5ad, backed="r").obs
    raw_dir = pathlib.Path(a.h5ad).parents[2] / "raw"
    pbs, rows, gene_ids = {}, [], None
    for d in sorted(obs.donor.astype(str).unique()):
        n = d.split()[-1]
        with h5py.File(raw_dir / f"GSE185224_Donor{n}_filtered_feature_bc_matrix.h5", "r") as f:
            g = f["matrix"]
            mat = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
            bcs = pd.Series(np.arange(mat.shape[1]), index=[b.decode() for b in g["barcodes"][:]])
            gex = np.array([t.decode() == "Gene Expression" for t in g["features/feature_type"][:]])
            ids = np.array([i.decode() for i in g["features/id"][:]])[gex]
        cells = pd.Series([c.rsplit("-", 1)[0] for c in obs.index[obs.donor.astype(str) == d]])
        found = cells[cells.isin(bcs.index)]
        if gene_ids is None:
            gene_ids = ids
        elif not (gene_ids == ids).all():
            raise SystemExit("Feature order differs between donors")
        unit = f"Burclaff_{d.replace(' ', '')}"
        pbs[unit] = np.asarray(mat[gex][:, bcs.loc[found].to_numpy()].sum(axis=1)).ravel().astype(np.int64)
        sub = obs[obs.donor.astype(str) == d]
        rows.append(dict(pseudobulk=unit, patient=unit, group="Adult", n_cells=int(len(found)),
                         annotated_cells=int(len(cells)),
                         regions=";".join(sorted(set(sub["region"].astype(str)))),
                         library_size=int(pbs[unit].sum())))
    adult = pd.DataFrame(pbs, index=gene_ids)
    adult = adult[~adult.index.duplicated(keep=False)]

    fetal = pd.read_csv(a.fetal_counts, index_col=0)
    fetal = fetal[~fetal.index.duplicated(keep=False)]
    shared = fetal.index.intersection(adult.index)
    merged = pd.concat([fetal.loc[shared], adult.loc[shared]], axis=1)
    merged.index.name = "gene_id"
    merged.to_csv(work / "Hnew2_Fawkner_fetal_Burclaff_adult_counts.csv.gz")

    units = pd.concat([pd.read_csv(a.fetal_units), pd.DataFrame(rows).assign(eligible=lambda x: x.n_cells >= a.min_cells)],
                      ignore_index=True)
    units.to_csv(work / "Hnew2_units.csv", index=False)
    print(f"{len(shared):,} shared Ensembl genes (fetal {len(fetal):,}, adult {len(adult):,})")
    print(pd.DataFrame(rows).to_string(index=False))
    print(units.groupby(["group", "eligible"]).size())


if __name__ == "__main__":
    main()
