#!/usr/bin/env python3
"""Cell-level objects for the fetal epithelial state analysis.

Fawkner-Corbett 2021: QC + hashtag calls exactly as benchmark v2
(02_fawkner_scrna_pseudobulk.py), hashtags mapped to samples with
config/sample_key_fawkner.tsv, then checked by sex (XIST vs Y genes) and
region (SATB2) and per-cell sex discordance (plan addendum v1.1); ineligible
units and sex-discordant cells are dropped. Gao 2018: author-labelled SI/LI epithelial cells (TPM).
Only the check genes named here are read at this stage (plan v1).
"""
import argparse
import importlib.util
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D"]


def load_v2_module(repo):
    path = repo / "step2_benchmark_v2/scripts/02_fawkner_scrna_pseudobulk.py"
    spec = importlib.util.spec_from_file_location("fawkner_v2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fawkner_cells(raw, repo, key):
    v2 = load_v2_module(repo)
    mats, obs, genes = [], [], None
    for gex_name, hto_name in v2.POOLS.items():
        pool = gex_name.split("_", 1)[1]
        gex, gbc, gfeat = v2.read_mtx(raw / f"{gex_name}.tar.gz", "raw_feature_bc_matrix")
        if genes is None:
            genes = gfeat
        umis = np.asarray(gex.sum(axis=0)).ravel()
        ngenes = np.asarray((gex > 0).sum(axis=0)).ravel()
        qc = (umis >= 500) & (ngenes >= 200)
        hto, hbc, hfeat = v2.read_mtx(raw / f"{hto_name}.tar.gz", "umi_count")
        tags = hfeat[0].astype(str).to_numpy()
        hto = hto[: len(tags)].toarray()
        srt = np.sort(hto, axis=0)
        called = (srt[-1] >= 10) & (srt[-1] >= 2 * srt[-2])
        call = pd.Series(np.where(called, tags[hto.argmax(axis=0)], ""), index=hbc)
        cell_tag = call.reindex(gbc).fillna("").to_numpy()
        use = qc & (cell_tag != "")
        # Hashtag label as in the key: "HTOk" for pools 1-3, leading letter for pool 4.
        lab = pd.Series(cell_tag[use]).str.split("-").str[0]
        if pool == "EPI_run2":
            lab = lab.str[0]
        mats.append(gex[:, use].T.tocsr())
        obs.append(pd.DataFrame({"pool": pool, "hashtag_totalseq": lab.to_numpy(),
                                 "n_umi": umis[use], "n_genes": ngenes[use]},
                                index=[f"{pool}_{b}" for b in gbc[use]]))
    o = pd.concat(obs)
    o = o.reset_index().merge(key, on=["pool", "hashtag_totalseq"], how="left").set_index("index")
    if o["sample"].isna().any():
        raise SystemExit("Hashtags missing from sample key: " + str(o[o["sample"].isna()][["pool", "hashtag_totalseq"]].drop_duplicates()))
    a = ad.AnnData(sp.vstack(mats).tocsr().astype(np.float32), obs=o)
    a.var_names = genes[1].astype(str).to_numpy()
    a.var["gene_id"] = genes[0].astype(str).to_numpy()
    a.var_names_make_unique()
    a.obs["dataset"] = "Fawkner2021"
    a.obs["unit"] = a.obs["sample"]
    return a


def cell_sex(a):
    X = a.X.tocsc()
    x = np.asarray(X[:, a.var_names.get_loc("XIST")].todense()).ravel()
    y = np.asarray(X[:, [a.var_names.get_loc(g) for g in Y_GENES]].sum(axis=1)).ravel()
    return np.where((x > 0) & (y == 0), "XX", np.where((y > 0) & (x == 0), "XY", "untyped"))


def fawkner_checks(a, min_cells, max_discordance):
    cpm = a.X.multiply(1e6 / np.asarray(a.X.sum(axis=1))).tocsr()
    rows = []
    for unit, idx in a.obs.groupby("unit", observed=True).indices.items():
        o = a.obs.iloc[idx[0]]
        val = lambda g: float(np.asarray(cpm[idx][:, a.var_names.get_loc(g)].mean()))
        xist, ychr, satb2 = val("XIST"), sum(val(g) for g in Y_GENES), val("SATB2")
        rows.append(dict(unit=unit, pool=o.pool, hashtag=o.hashtag_totalseq, donor=o.donor, pcw=o.pcw,
                         region=o.region, key_sex=o.sex, n_cells=len(idx), XIST_cpm=xist,
                         Ygenes_cpm=ychr, SATB2_cpm=satb2,
                         observed_sex="XX" if xist > ychr else "XY"))
    c = pd.DataFrame(rows)
    typed = a.obs[a.obs.cell_sex != "untyped"]
    disc = (typed.cell_sex != typed.sex).groupby(typed.unit, observed=True).mean()
    c["cell_sex_discordant_frac"] = c.unit.map(disc)
    c["sex_check"] = c.observed_sex == c.key_sex
    # Region (addendum v1.1): within donor, among units with >= min_cells,
    # colon/hindgut SATB2 must exceed every TI sample.
    def region_ok(r):
        if r.n_cells < min_cells:
            return np.nan
        d = c[(c.donor == r.donor) & (c.n_cells >= min_cells)]
        ti, col = d[d.region == "TI"].SATB2_cpm, d[d.region != "TI"].SATB2_cpm
        if len(ti) == 0 or len(col) == 0:
            return np.nan
        return bool(r.SATB2_cpm > ti.max()) if r.region != "TI" else bool(r.SATB2_cpm < col.min())
    c["region_check"] = c.apply(region_ok, axis=1)
    c["discordance_check"] = c.cell_sex_discordant_frac <= max_discordance
    c["eligible"] = (c.n_cells >= min_cells) & c.sex_check & (c.region_check != False) & c.discordance_check
    c["flag"] = np.where(c.unit == "ABZ2", "elevated discordance; sensitivity excludes", "")
    return c


def gao_cells(gao_dir):
    ann = pd.read_excel(gao_dir / "41556_2018_105_MOESM4_ESM.xlsx")
    m = ann[ann.Tissue.isin(["SI", "LI"]) & (ann.CellType == "Epithelial")].copy()
    m["week"] = m.Sample.str.extract(r"_(\d+)W_")[0].astype(int)
    m["donor"] = m.Sample.str.extract(r"_(\d+W_embryo\d+)")[0]
    m["region"] = m.Tissue
    m["unit"] = m.donor + "_" + m.region
    m["matrix_sample"] = m.Sample.str.replace(r"_(\d+)W_", r"_\1w_", regex=True)
    expr = pd.read_csv(gao_dir / "GSE95630_Digestion_TPM_new.txt.gz", sep="\t",
                       usecols=["Gene"] + m.matrix_sample.tolist()).set_index("Gene")
    expr = expr[~expr.index.duplicated(keep="first")]
    a = ad.AnnData(sp.csr_matrix(expr[m.matrix_sample].to_numpy(dtype=np.float32).T),
                   obs=m.set_index("Sample")[["Tissue", "Group", "week", "donor", "region", "unit"]])
    a.var_names = expr.index.astype(str)
    a.obs["dataset"] = "Gao2018"
    a.obs["pcw"] = a.obs.week
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--fawkner-raw", required=True)
    ap.add_argument("--gao-raw", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--min-cells", type=int, default=50)
    ap.add_argument("--max-discordance", type=float, default=0.10)
    a = ap.parse_args()
    repo, work, tables = map(pathlib.Path, (a.repo, a.work_dir, a.tables_dir))
    work.mkdir(parents=True, exist_ok=True)
    key = pd.read_csv(repo / "step2_fetal_state/config/sample_key_fawkner.tsv", sep="\t")
    f = fawkner_cells(pathlib.Path(a.fawkner_raw), repo, key)
    f.obs["cell_sex"] = cell_sex(f)
    checks = fawkner_checks(f, a.min_cells, a.max_discordance)
    checks.to_csv(tables / "Fawkner_hashtag_sample_checks.csv", index=False)
    pd.set_option("display.width", 220)
    print(checks.round(1).to_string(index=False))
    f.obs["unit_eligible"] = f.obs.unit.map(checks.set_index("unit").eligible).astype(bool)
    # Addendum v1.1: drop sex-discordant cells and ineligible units.
    keep = f.obs.unit_eligible & ((f.obs.cell_sex == "untyped") | (f.obs.cell_sex == f.obs.sex))
    f = f[keep.to_numpy()].copy()
    f.write_h5ad(work / "fawkner_epithelial_cells.h5ad")
    g = gao_cells(pathlib.Path(a.gao_raw))
    g.write_h5ad(work / "gao_epithelial_cells.h5ad")
    print("Fawkner", f.shape, "Gao", g.shape)


if __name__ == "__main__":
    main()
