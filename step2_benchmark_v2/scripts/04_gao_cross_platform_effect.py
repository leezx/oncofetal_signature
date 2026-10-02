#!/usr/bin/env python3
"""H-new3, H-new3-debug and Gao-original (recomputed): cross-platform,
effect-only fetal vs adult contrasts with Gao 2018 fetal epithelium.

Unit value = mean per-cell expression on a per-million scale: TPM for Gao
fetal (GSE95630), UMI CPM for Wang 2020 (GSE125970), the released UMI-TPM for
GSE103154. Effect = log2((median fetal unit + 1) / (median adult unit + 1)).
Support = Welch t-test on log2(unit + 1); platform is confounded with stage,
so P is descriptive. BH FDR over all genes present in both arms and non-zero
in at least one unit. Rules frozen in docs/QUALIFICATION_PLAN.md.
"""
import argparse
import gzip
import pathlib

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests


def gao_fetal_units(gao_dir, tissues):
    ann = pd.read_excel(gao_dir / "41556_2018_105_MOESM4_ESM.xlsx")
    meta = ann[ann.Tissue.isin(tissues) & (ann.CellType == "Epithelial")].copy()
    parsed = meta.Sample.str.extract(r"_(\d+W)_embryo(\d+)(?:\.\d+)?_")
    if parsed.isna().any().any():
        raise SystemExit("Could not parse embryo IDs")
    meta["unit"] = parsed[0] + "_embryo" + parsed[1]
    meta["matrix_sample"] = meta.Sample.str.replace(r"_(\d+)W_", r"_\1w_", regex=True)
    path = gao_dir / "GSE95630_Digestion_TPM_new.txt.gz"
    expr = pd.read_csv(path, sep="\t", usecols=["Gene"] + meta.matrix_sample.tolist()).set_index("Gene")
    units = pd.DataFrame({u: expr[g.matrix_sample].mean(axis=1) for u, g in meta.groupby("unit")})
    info = meta.groupby("unit").agg(n_cells=("Sample", "size"),
                                    tissues=("Tissue", lambda s: ";".join(sorted(set(s))))).reset_index()
    return units, info


def wang_adult_units(wang_dir, samples):
    info = pd.read_csv(wang_dir / "GSE125970_cell_info.txt.gz", sep="\t")
    info = info[info.Sample_ID.isin(samples)]
    # Chunked sparse read: a column-filtered read of the 14,537-cell matrix is
    # prohibitively slow. Every matrix column is an annotated epithelial cell.
    blocks, genes, cells = [], [], None
    for ch in pd.read_csv(wang_dir / "GSE125970_raw_UMIcounts.txt.gz", sep="\t", index_col=0, chunksize=2000):
        cells = ch.columns
        genes.extend(ch.index)
        blocks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32)))
    mat = sp.vstack(blocks).tocsc()
    cpm = mat.multiply(1e6 / np.asarray(mat.sum(axis=0)).ravel()).tocsc()
    col = pd.Series(np.arange(len(cells)), index=cells)
    units = pd.DataFrame({s: np.asarray(cpm[:, col.loc[g.UniqueCell_ID].to_numpy()].mean(axis=1)).ravel()
                          for s, g in info.groupby("Sample_ID")}, index=genes)
    units = units[~units.index.duplicated(keep=False)]
    meta = info.groupby("Sample_ID").size().rename("n_cells").reset_index().rename(columns={"Sample_ID": "unit"})
    return units, meta


def gse103154_units(gao_dir):
    expr = pd.read_csv(gao_dir / "GSE103154_All_Merge_umi_tpm_gene.txt.gz", sep="\t").set_index("Gene")
    units = pd.DataFrame({d: expr[[c for c in expr.columns if c.startswith(d + "_")]].mean(axis=1)
                          for d in ("P1", "P2")})
    meta = pd.DataFrame({"unit": ["P1", "P2"],
                         "n_cells": [sum(c.startswith(d + "_") for c in expr.columns) for d in ("P1", "P2")]})
    return units, meta


def effect(fetal, adult):
    genes = fetal.index.intersection(adult.index)
    f, a = fetal.loc[genes], adult.loc[genes]
    expressed = (f.sum(axis=1) + a.sum(axis=1)) > 0
    f, a = f[expressed], a[expressed]
    lf, la = np.log2(f + 1), np.log2(a + 1)
    t = stats.ttest_ind(lf, la, axis=1, equal_var=False)
    p = np.nan_to_num(t.pvalue, nan=1.0)
    res = pd.DataFrame({
        "symbol": f.index,
        "log2FC": np.log2((f.median(axis=1) + 1) / (a.median(axis=1) + 1)).to_numpy(),
        "PValue": p,
        "fetal_median": f.median(axis=1).to_numpy(),
        "adult_median": a.median(axis=1).to_numpy(),
        "fetal_units_above_adult_max": (f.gt(a.max(axis=1), axis=0)).sum(axis=1).to_numpy(),
        "n_fetal_units": f.shape[1], "n_adult_units": a.shape[1],
    })
    res["FDR"] = multipletests(res.PValue, method="fdr_bh")[1]
    return res.sort_values("PValue"), f, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gao-dir", required=True)
    ap.add_argument("--wang-dir", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--aliases", required=True)
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    a = ap.parse_args()
    gao, wang = pathlib.Path(a.gao_dir), pathlib.Path(a.wang_dir)
    de, tab = pathlib.Path(a.de_dir), pathlib.Path(a.tables_dir)
    de.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(a.candidates, sep="\t").gene
    alias = pd.read_csv(a.aliases, sep="\t")
    keys = {g: [g] + alias.loc[alias.current_hgnc_symbol == g, "source_symbol"].tolist() for g in cand}

    fetal_all, fetal_all_info = gao_fetal_units(gao, ["SI", "LI"])
    fetal_li, fetal_li_info = gao_fetal_units(gao, ["LI"])
    wang_all, wang_all_info = wang_adult_units(wang, ["Ileum-1", "Ileum-2", "Colon-1", "Colon-2", "Rectum-1", "Rectum-2"])
    wang_lc = wang_all[["Colon-1", "Colon-2", "Rectum-1", "Rectum-2"]]
    adult_gao, adult_gao_info = gse103154_units(gao)

    contrasts = {
        "Hnew3_GaoSILI_vs_Wang": (fetal_all, wang_all, fetal_all_info, wang_all_info),
        "Hnew3debug_GaoLI_vs_WangColonRectum": (fetal_li, wang_lc, fetal_li_info,
                                                 wang_all_info[wang_all_info.unit.isin(wang_lc.columns)]),
        "GaoOriginal_GaoLI_vs_GSE103154": (fetal_li, adult_gao, fetal_li_info, adult_gao_info),
    }
    for name, (f, ad, fi, ai) in contrasts.items():
        res, fu, au = effect(f, ad)
        res.to_csv(de / f"{name}.csv", index=False)
        units = pd.concat([fi.assign(group="Fetal"), ai.assign(group="Adult")], ignore_index=True)
        units.to_csv(tab / f"{name}_units.csv", index=False)
        rows = []
        for g, ks in keys.items():
            hit = next((k for k in ks if k in fu.index), None)
            if hit is None:
                continue
            rows.append(pd.concat([fu.loc[hit], au.loc[hit]]).rename(g))
        pd.DataFrame(rows).round(3).to_csv(tab / f"{name}_marker_unit_values.csv")
        tn = res.set_index("symbol").loc["TNFRSF12A"]
        print(f"{name}: {len(res):,} genes; fetal units {fu.shape[1]}, adult units {au.shape[1]}; "
              f"TNFRSF12A log2FC {tn.log2FC:.2f} P {tn.PValue:.3g}")


if __name__ == "__main__":
    main()
