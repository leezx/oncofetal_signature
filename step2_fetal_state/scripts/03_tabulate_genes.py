#!/usr/bin/env python3
"""Tabulate the analysed genes on frozen state labels (plan v1, addenda v1.1-1.2).

Unit = donor x region (Fawkner sample; Gao embryo x SI/LI).
- unit / unit x state value: Fawkner pseudobulk log2(CPM + 1); Gao log2(mean TPM + 1)
- % positive: cells with count/TPM > 0
- decomposition: unit mean per-cell expression (Fawkner CP10k, Gao TPM) equals
  sum_s fraction_s x state mean_s exactly; counterfactuals hold either the
  composition or the state means at the dataset average.
Descriptive statistics only (Spearman with age; within-donor region differences).
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

MIN_UNIT = {"Fawkner2021": 50, "Gao2018": 20}
MIN_UNIT_STATE = {"Fawkner2021": 10, "Gao2018": 5}
SENSITIVITY_DROP = {"Fawkner2021": ["ABZ2"]}


def resolve(names, genes):
    out = {}
    for _, r in genes.iterrows():
        sym = r.symbol if r.symbol in names else (r.alias if r.alias and r.alias in names else None)
        if sym is None:
            raise SystemExit(f"{r.symbol} not found")
        out[r.symbol] = sym
    return out


def dense(x):
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def tabulate(name, a, labels, gmap):
    a = a[labels.index].copy()
    a.obs["state"] = labels.state.to_numpy()
    cols = [gmap[g] for g in gmap]
    X = dense(a[:, cols].X)
    if name == "Fawkner2021":
        lib = np.asarray(a.X.sum(axis=1)).ravel()
        per_cell = X / lib[:, None] * 1e4          # CP10k for the decomposition
    else:
        lib = None
        per_cell = X                                # TPM
    obs = a.obs.copy()
    meta_cols = ["dataset", "donor", "pcw", "region"]

    def summarise(idx):
        r = {}
        if name == "Fawkner2021":
            cpm = X[idx].sum(axis=0) / lib[idx].sum() * 1e6
            r.update({f"{g}": np.log2(v + 1) for g, v in zip(gmap, cpm)})
        else:
            r.update({f"{g}": np.log2(v + 1) for g, v in zip(gmap, X[idx].mean(axis=0))})
        r.update({f"{g}_pct_pos": 100 * v for g, v in zip(gmap, (X[idx] > 0).mean(axis=0))})
        r.update({f"{g}_mean_percell": v for g, v in zip(gmap, per_cell[idx].mean(axis=0))})
        return r

    units, unit_state = [], []
    for unit, idx in obs.groupby("unit", observed=True).indices.items():
        o = obs.iloc[idx[0]]
        base = dict(dataset=name, unit=unit, **{c: o[c] for c in meta_cols[1:]}, n_cells=len(idx))
        units.append({**base, **summarise(idx)})
        for st, sidx in obs.iloc[idx].groupby("state", observed=True).indices.items():
            unit_state.append({**base, "state": st, "n_cells_state": len(sidx),
                               "state_fraction": len(sidx) / len(idx), **summarise(idx[sidx])})
    u = pd.DataFrame(units)
    us = pd.DataFrame(unit_state)
    u["eligible"] = u.n_cells >= MIN_UNIT[name]
    us["eligible"] = us.unit.map(u.set_index("unit").eligible) & (us.n_cells_state >= MIN_UNIT_STATE[name])
    u["in_sensitivity"] = u.eligible & ~u.unit.isin(SENSITIVITY_DROP.get(name, []))
    return u, us


def decomposition(name, u, us, gene, keep):
    """Observed vs counterfactual unit means (all states incl. Unresolved)."""
    us = us[us.unit.isin(keep)]
    f = us.pivot(index="unit", columns="state", values="state_fraction").fillna(0)
    m = us.pivot(index="unit", columns="state", values=f"{gene}_mean_percell")
    mbar = m.mean(axis=0)
    m_filled = m.apply(lambda c: c.fillna(mbar[c.name]))
    fbar = f.mean(axis=0)
    obs = (f * m_filled).sum(axis=1)
    fixed_comp = (m_filled * fbar).sum(axis=1)
    fixed_expr = (f * mbar).sum(axis=1)
    d = pd.DataFrame({"observed": obs, "fixed_composition": fixed_comp, "fixed_expression": fixed_expr})
    d = d.join(u.set_index("unit")[["donor", "pcw", "region"]])
    d.insert(0, "gene", gene)
    d.insert(0, "dataset", name)
    return d.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    a = ap.parse_args()
    repo, work, out = pathlib.Path(a.repo), pathlib.Path(a.work_dir), pathlib.Path(a.tables_dir)
    genes = pd.read_csv(repo / "step2_fetal_state/config/genes.tsv", sep="\t", keep_default_na=False)
    U, US, DEC, COR, REG = [], [], [], [], []
    for name, f in (("Fawkner2021", "fawkner_epithelial_cells.h5ad"), ("Gao2018", "gao_epithelial_cells.h5ad")):
        cells = ad.read_h5ad(work / f)
        labels = pd.read_csv(work / f"{name}_state_labels.csv.gz", index_col=0)
        gmap = resolve(set(cells.var_names), genes)
        u, us = tabulate(name, cells, labels, gmap)
        U.append(u); US.append(us)
        for scope, keep in (("all", u[u.eligible].unit), ("sensitivity", u[u.in_sensitivity].unit)):
            if scope == "sensitivity" and name not in SENSITIVITY_DROP:
                continue
            ue = u[u.unit.isin(keep)]
            for g in gmap:
                rho, p = stats.spearmanr(ue.pcw, ue[g])
                COR.append(dict(dataset=name, scope=scope, level="unit", state="all", gene=g,
                                n_units=len(ue), n_donors=ue.donor.nunique(), spearman_rho_age=rho, P=p))
                use = us[us.eligible & us.unit.isin(keep)]
                for st, s in use.groupby("state"):
                    if len(s) >= 4 and s.pcw.nunique() > 2:
                        rho, p = stats.spearmanr(s.pcw, s[g])
                        COR.append(dict(dataset=name, scope=scope, level="unit_x_state", state=st, gene=g,
                                        n_units=len(s), n_donors=s.donor.nunique(), spearman_rho_age=rho, P=p))
            for g in ("TNFRSF12A", "TACSTD2", "CLU", "ANXA1"):
                d = decomposition(name, u, us, g, keep)
                d.insert(2, "scope", scope)
                DEC.append(d)
        # Within-donor region contrasts (eligible units): colon/LI minus TI/SI.
        ue = u[u.eligible]
        small = {"Fawkner2021": "TI", "Gao2018": "SI"}[name]
        for donor, d in ue.groupby("donor"):
            if (d.region == small).sum() == 1 and (d.region != small).any():
                ti = d[d.region == small].iloc[0]
                for _, r in d[d.region != small].iterrows():
                    REG.append(dict(dataset=name, donor=donor, pcw=r.pcw, region=r.region, vs=small,
                                    **{g: r[g] - ti[g] for g in gmap}))
    u, us = pd.concat(U), pd.concat(US)
    u.round(4).to_csv(out / "Unit_gene_values.csv", index=False)
    us.round(4).to_csv(out / "Unit_state_gene_values.csv", index=False)
    pd.concat(DEC).round(4).to_csv(out / "Composition_decomposition.csv", index=False)
    pd.DataFrame(COR).round(4).to_csv(out / "Age_correlations.csv", index=False)
    pd.DataFrame(REG).round(4).to_csv(out / "Within_donor_region_differences.csv", index=False)
    # Summaries: decomposition variance (log2 with pseudocount = 1% of median
    # observed unit value), early (<= 8 weeks) vs later units, co-variation.
    rows = []
    for (d, scope, g), x in pd.concat(DEC).groupby(["dataset", "scope", "gene"]):
        pc = 0.01 * x.observed.median()
        r = dict(dataset=d, scope=scope, gene=g, n_units=len(x))
        for k in ("observed", "fixed_composition", "fixed_expression"):
            r[f"var_log2_{k}"] = np.var(np.log2(x[k] + pc))
            r[f"rho_age_{k}"] = stats.spearmanr(x.pcw, x[k])[0]
        rows.append(r)
    pd.DataFrame(rows).round(4).to_csv(out / "Composition_decomposition_summary.csv", index=False)
    rows = []
    for d, x in u[u.eligible].groupby("dataset"):
        for g in [c for c in genes.symbol]:
            e, l = x[x.pcw <= 8][g], x[x.pcw >= 9][g]
            rows.append(dict(dataset=d, gene=g, n_units_le8=len(e), median_le8=e.median(), n_units_ge9=len(l),
                             median_ge9=l.median(), fold_range_all=2 ** (x[g].max() - x[g].min()),
                             fold_range_ge9=2 ** (l.max() - l.min()),
                             rho_age_ge9=stats.spearmanr(x[x.pcw >= 9].pcw, l)[0],
                             rho_with_TNFRSF12A=stats.spearmanr(x.TNFRSF12A, x[g])[0]))
    pd.DataFrame(rows).round(4).to_csv(out / "Early_vs_later_units.csv", index=False)
    comp = us.pivot_table(index=["dataset", "unit", "donor", "pcw", "region"], columns="state",
                          values="state_fraction", fill_value=0).round(4)
    comp.to_csv(out / "Unit_state_composition.csv")
    print("units", u.groupby("dataset").eligible.sum().to_dict())


if __name__ == "__main__":
    main()
