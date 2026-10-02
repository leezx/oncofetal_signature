#!/usr/bin/env python3
"""Addendum v1.4: Core-gate dataset selection benchmark.

Design eligibility (E1-E5) first, 31-marker literature recovery second.
Values come from the public marker summary (no recomputation). TNFRSF12A is
not an admission criterion. No gate is selected here.
Usage (repo root): python3 step2_benchmark_v2/scripts/11_core_gate_recovery.py
"""
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "marker_summary/results/Literature_31_marker_summary_statistics_public.csv"
OUTS = [ROOT / "step2_benchmark_v2/results/tables", ROOT / "marker_summary/results"]

# key, species, role, description, E1 same study, E2 in vivo, E3 epithelium, E4 >=3/arm, E5 replicate P/FDR, n fetal, n adult, note
C = [
    ("HGCA_late", "Human", "candidate A", "HGCA fetal >= 9 PCW vs adult epithelium (donor pseudobulk, edgeR)",
     1, 1, 1, 1, 1, 10, 7, "Gut Cell Atlas; fetal (HDBR) and adult (transplant donors) cohorts collected within one atlas framework"),
    ("HGCA", "Human", "reference (all fetal)", "HGCA all fetal vs adult epithelium",
     1, 1, 1, 1, 1, 16, 7, "includes 6 donors < 9 PCW"),
    ("GaoOriginal_late", "Human", "candidate B", "Gao fetal LI >= 9 W vs GSE103154 adult LI (effect-only)",
     1, 1, 1, 0, 0, 7, 2, "same research programme but fetal TPM vs adult UMI platform; adult n = 2; Welch P on 7 vs 2"),
    ("Hnew2_late", "Human", "candidate C", "Fawkner fetal >= 9 PCW vs Burclaff adult (donor pseudobulk, edgeR)",
     0, 1, 1, 1, 1, 20, 3, "cross-study: stage collinear with study/platform"),
    ("Hnew2", "Human", "reference (all fetal)", "Fawkner all fetal vs Burclaff adult",
     0, 1, 1, 1, 1, 22, 3, "cross-study"),
    ("Hnew3_late", "Human", "reference", "Gao fetal SI+LI >= 9 W vs Wang adult (effect-only)",
     0, 1, 1, 1, 0, 9, 6, "cross-study, cross-platform"),
    ("Hbulk1", "Human", "reference", "Roadmap fetal SI vs HPA adult SI (bulk)",
     0, 1, 0, 1, 1, 6, 6, "cross-study whole tissue"),
    ("Senger", "Human", "reference", "Senger fetal vs adult enterospheres",
     1, 0, 1, 1, 1, 6, 3, "culture"),
    ("Mouse_GSE230581", "Mouse", "candidate (Pikkupeura in vivo)", "Pikkupeura E16.5 proximal SI epithelium vs adult crypt epithelium (in vivo, edgeR)",
     1, 1, 1, 1, 1, 3, 3, "= Pikkupeura 2023 Fig 1D-G; adult arm is crypt-only epithelium"),
    ("GSE44433", "Mouse", "candidate (independent in vivo)", "Hemmerling WT E17.5 vs 8-week LCM ileal epithelium (microarray)",
     1, 1, 1, 1, 1, 5, 5, "older microarray; incomplete probe coverage"),
    ("Pikkupeura_LN", "Mouse", "reference (culture)", "Pikkupeura FEnS vs adult organoids, laminin",
     1, 0, 1, 1, 1, np.nan, np.nan, "culture; FDR only"),
    ("Pikkupeura_collagen", "Mouse", "reference (culture)", "Pikkupeura FEnS vs adult organoids, collagen",
     1, 0, 1, 1, 1, np.nan, np.nan, "culture; FDR only"),
]
NOT_APPLICABLE = {"Human": {"LY6A", "REG3B"}, "Mouse": {"SPRR1A"}}


def main():
    s = pd.read_csv(SUMMARY).set_index("marker")
    rows, per = [], []
    for key, sp, role, desc, e1, e2, e3, e4, e5, nf, na, note in C:
        lfc, fdr = s[f"{key}__log2FC"], s[f"{key}__FDR"]
        applicable = [g for g in s.index if g not in NOT_APPLICABLE[sp]]
        lfc, fdr = lfc[applicable], fdr[applicable]
        meas = lfc.notna()
        pos = meas & (lfc > 0)
        pos_sig = pos & (fdr < 0.05)
        pos_sig_eff = pos_sig & (lfc >= 0.5)
        n_app, n_meas = len(applicable), int(meas.sum())
        elig = all([e1, e2, e3, e4, e5])
        r = dict(contrast=key, species=sp, role=role, description=desc,
                 E1_same_study=bool(e1), E2_in_vivo=bool(e2), E3_epithelium_resolved=bool(e3),
                 E4_ge3_replicates_per_arm=bool(e4), E5_replicate_level_FDR=bool(e5),
                 design_eligible=elig, n_fetal=nf, n_adult=na, design_note=note,
                 panel_applicable=n_app, markers_measured=n_meas)
        for lab, m in (("fetal_positive", pos), ("fetal_positive_FDR05", pos_sig),
                       ("fetal_positive_FDR05_lfc0.5", pos_sig_eff)):
            r[f"n_{lab}"] = int(m.sum())
            r[f"pct_{lab}_of_measured"] = round(100 * m.sum() / n_meas, 1) if n_meas else np.nan
            r[f"pct_{lab}_of_applicable"] = round(100 * m.sum() / n_app, 1)
        r["n_adult_high_FDR05"] = int((meas & (lfc < 0) & (fdr < 0.05)).sum())
        r["TNFRSF12A_log2FC_info_only"] = s.loc["TNFRSF12A", f"{key}__log2FC"]
        r["TNFRSF12A_FDR_info_only"] = s.loc["TNFRSF12A", f"{key}__FDR"]
        rows.append(r)
        for g in s.index:
            v, q = s.loc[g, f"{key}__log2FC"], s.loc[g, f"{key}__FDR"]
            if g in NOT_APPLICABLE[sp]:
                cls = "not applicable (no one2one orthologue)"
            elif pd.isna(v):
                cls = "not measured"
            elif v > 0 and q < 0.05 and v >= 0.5:
                cls = "fetal-high (FDR<0.05, log2FC>=0.5)"
            elif v > 0 and q < 0.05:
                cls = "fetal-high (FDR<0.05, log2FC<0.5)"
            elif v > 0:
                cls = "fetal-positive, n.s." if pd.notna(q) else "fetal-positive (no FDR)"
            elif q < 0.05:
                cls = "adult-high (FDR<0.05)"
            else:
                cls = "adult-direction, n.s."
            per.append(dict(marker=g, contrast=key, species=sp, design_eligible=elig,
                            log2FC=v, PValue=s.loc[g, f"{key}__PValue"], FDR=q, class_=cls))
    rec = pd.DataFrame(rows)
    p = pd.DataFrame(per).rename(columns={"class_": "class"})
    wide = p.pivot(index="marker", columns="contrast", values=["log2FC", "FDR", "class"])
    wide.columns = [f"{c}__{s_}" for s_, c in wide.columns]
    order = [k for k, *_ in C]
    wide = wide[[f"{k}__{s_}" for k in order for s_ in ("log2FC", "FDR", "class")]].reindex(s.index)
    for out in OUTS:
        rec.to_csv(out / "Core_gate_dataset_recovery.csv", index=False)
        wide.round(4).to_csv(out / "Core_gate_marker_by_dataset.csv")
    pd.set_option("display.width", 250)
    print(rec[["contrast", "design_eligible", "markers_measured", "n_fetal_positive", "pct_fetal_positive_of_measured",
               "n_fetal_positive_FDR05", "pct_fetal_positive_FDR05_of_measured", "n_fetal_positive_FDR05_lfc0.5",
               "pct_fetal_positive_FDR05_lfc0.5_of_measured", "n_adult_high_FDR05"]].to_string(index=False))


if __name__ == "__main__":
    main()
