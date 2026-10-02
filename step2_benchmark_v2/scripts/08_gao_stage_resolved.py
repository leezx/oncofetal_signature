#!/usr/bin/env python3
"""Addendum v1.2: stage-resolved Gao contrasts (breakpoint 9 weeks, frozen).

Reuses module 04 unchanged (units, effect statistic). The fetal arm is split
by the reported embryo week: early < 9 W, mid/late >= 9 W.
  GaoOriginal_{early,late}: Gao fetal LI vs GSE103154 adult
  Hnew3_{early,late}:       Gao fetal SI+LI vs Wang adult
  GaoLI_late_vs_early, GaoSILI_late_vs_early: within-fetal (same platform)
"""
import argparse
import importlib.util
import pathlib

import pandas as pd

BREAK = 9


def load_04(script_dir):
    spec = importlib.util.spec_from_file_location("gao04", script_dir / "04_gao_cross_platform_effect.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def split(units, info):
    week = info.set_index("unit").index.str.extract(r"^(\d+)W", expand=False).astype(int)
    wk = pd.Series(week.to_numpy(), index=info.unit)
    early, late = wk[wk < BREAK].index, wk[wk >= BREAK].index
    return (units[early], info[info.unit.isin(early)].assign(week=wk[early].to_numpy())), \
           (units[late], info[info.unit.isin(late)].assign(week=wk[late].to_numpy()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gao-dir", required=True)
    ap.add_argument("--wang-dir", required=True)
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    a = ap.parse_args()
    m = load_04(pathlib.Path(__file__).resolve().parent)
    gao, wang = pathlib.Path(a.gao_dir), pathlib.Path(a.wang_dir)
    de, tab = pathlib.Path(a.de_dir), pathlib.Path(a.tables_dir)
    li, li_info = m.gao_fetal_units(gao, ["LI"])
    sili, sili_info = m.gao_fetal_units(gao, ["SI", "LI"])
    adult_gao, adult_gao_info = m.gse103154_units(gao)
    wang_all, wang_info = m.wang_adult_units(wang, ["Ileum-1", "Ileum-2", "Colon-1", "Colon-2", "Rectum-1", "Rectum-2"])
    (li_e, li_e_i), (li_l, li_l_i) = split(li, li_info)
    (si_e, si_e_i), (si_l, si_l_i) = split(sili, sili_info)
    contrasts = {
        "GaoOriginal_early_GaoLI_vs_GSE103154": (li_e, adult_gao, li_e_i.assign(group="Fetal_early"), adult_gao_info.assign(group="Adult")),
        "GaoOriginal_late_GaoLI_vs_GSE103154": (li_l, adult_gao, li_l_i.assign(group="Fetal_midlate"), adult_gao_info.assign(group="Adult")),
        "Hnew3_early_GaoSILI_vs_Wang": (si_e, wang_all, si_e_i.assign(group="Fetal_early"), wang_info.assign(group="Adult")),
        "Hnew3_late_GaoSILI_vs_Wang": (si_l, wang_all, si_l_i.assign(group="Fetal_midlate"), wang_info.assign(group="Adult")),
        "GaoLI_late_vs_early": (li_l, li_e, li_l_i.assign(group="Fetal_midlate"), li_e_i.assign(group="Fetal_early")),
        "GaoSILI_late_vs_early": (si_l, si_e, si_l_i.assign(group="Fetal_midlate"), si_e_i.assign(group="Fetal_early")),
    }
    for name, (case, ref, ci, ri) in contrasts.items():
        res, fu, au = m.effect(case, ref)
        res.to_csv(de / f"{name}.csv", index=False)
        pd.concat([ci, ri], ignore_index=True).to_csv(tab / f"{name}_units.csv", index=False)
        tn = res.set_index("symbol").loc["TNFRSF12A"]
        print(f"{name}: case {fu.shape[1]} vs ref {au.shape[1]} units; TNFRSF12A {tn.log2FC:+.2f} P {tn.PValue:.3g}")


if __name__ == "__main__":
    main()
