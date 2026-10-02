#!/usr/bin/env python3
"""Joanito cohort x group confounding check (review decision v3).

`prepare`: write the patient-level cohort x group contingency table and a
sensitivity sample sheet restricted to cohorts that contain both eligible
malignant and eligible normal patients.
`compare`: compare full-model (~ cohort + group) and sensitivity S1 results.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def prepare(a):
    smp = pd.read_csv(a.samples)
    el = smp[smp.eligible]
    tab = el.groupby(["cohort", "group"]).size().unstack(fill_value=0)
    tab = tab.rename(columns={"Malignant": "malignant_n", "Normal": "normal_n"})
    both = el.pivot_table(index=["cohort", "patient"], columns="group", values="n_cells").dropna()
    tab["paired_patients"] = both.reset_index().groupby("cohort").size().reindex(tab.index, fill_value=0)
    tab["in_sensitivity"] = (tab.malignant_n > 0) & (tab.normal_n > 0)
    tab.reset_index().to_csv(a.contingency, index=False)
    keep = tab.index[tab.in_sensitivity]
    smp.assign(eligible=smp.eligible & smp.cohort.isin(keep)).to_csv(a.sensitivity_samples, index=False)
    print(tab.to_string())


def compare(a):
    full = pd.read_csv(a.full).set_index("gene_id")
    sens = pd.read_csv(a.sensitivity).set_index("gene_id")
    j = full.join(sens, lsuffix="_full", rsuffix="_sens", how="inner")
    pass_full = (j.log2FC_full >= a.min_log2fc) & (j.FDR_full < a.max_fdr)
    pass_sens = (j.log2FC_sens >= a.min_log2fc) & (j.FDR_sens < a.max_fdr)
    rho = spearmanr(j.log2FC_full, j.log2FC_sens).statistic
    out = pd.DataFrame([{
        "n_jointly_tested": len(j),
        "spearman_log2FC": rho,
        "S1_pass_full": int(pass_full.sum()),
        "S1_pass_sensitivity": int(pass_sens.sum()),
        "S1_pass_both": int((pass_full & pass_sens).sum()),
        "full_pass_retained_fraction": (pass_full & pass_sens).sum() / max(pass_full.sum(), 1),
        "full_pass_same_sign_in_sensitivity": float((np.sign(j.log2FC_sens[pass_full]) > 0).mean()),
    }])
    out.to_csv(a.out, index=False)
    print(out.T.to_string(header=False))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--samples", required=True)
    p.add_argument("--contingency", required=True)
    p.add_argument("--sensitivity-samples", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--full", required=True)
    c.add_argument("--sensitivity", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--min-log2fc", type=float, default=0.5)
    c.add_argument("--max-fdr", type=float, default=0.05)
    a = ap.parse_args()
    prepare(a) if a.cmd == "prepare" else compare(a)


if __name__ == "__main__":
    main()
