#!/usr/bin/env python3
"""Core oncofetal workbook (method v2.0, docs/Method_core_oncofetal_validated.md).

Axes: L literature (config/literature_provenance_31.tsv); H human developmental
replicated (HGCA >=9 PCW, H-new2 >=9 PCW, Gao >=9 W: >=2/3 fetal-positive and
>=1 with log2FC >= 0.5 & FDR < 0.05); M mouse in vivo GSE230581 (log2FC >= 0.5
& FDR < 0.05); C CRC replicated (Joanito and Pelka both log2FC >= 0.5 &
FDR < 0.05). Core = L & H & M & C pass. Calls: pass / fail / not evaluable (NE).
Stringent intersection set (former v1.0) is reported as a sensitivity subset.
Values come from the marker summary (no recomputation). The full workbook
(Joanito, C call, Core labels) is restricted and git-ignored; the public one
omits all Joanito-derived content.
Usage (repo root): python3 core_oncofetal/scripts/build_core_gates.py
"""
import pathlib
import shutil

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

ROOT = pathlib.Path(__file__).resolve().parents[2]
MS = ROOT / "marker_summary/results/Literature_31_marker_summary_statistics.csv"
PROV = ROOT / "core_oncofetal/config/literature_provenance_31.tsv"
OUT = ROOT / "core_oncofetal/results"
RESTRICTED = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/results/"
                          "2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/core_oncofetal")
LFC, FDR = 0.5, 0.05
NA_HUMAN, NA_MOUSE = {"LY6A", "REG3B"}, {"SPRR1A"}

# key, column label, species, axis/role, restricted
EVIDENCE = [
    ("HGCA_late", "HGCA >=9 PCW vs adult (human, primary)", "Human", "H", False),
    ("Hnew2_late", "H-new2 Fawkner >=9 PCW vs Burclaff adult (human)", "Human", "H", False),
    ("GaoOriginal_late", "Gao LI >=9 W vs GSE103154 adult (human)", "Human", "H", False),
    ("Mouse_GSE230581", "GSE230581 E16.5 vs adult crypt (mouse in vivo)", "Mouse", "M", False),
    ("GSE44433", "GSE44433 E17.5 vs 8-week LCM (mouse replication)", "Mouse", "support", False),
    ("Joanito", "Joanito malignant vs normal epithelium (CRC)", "Human", "C", True),
    ("Pelka", "Pelka tumour vs normal epithelium (CRC replication)", "Human", "C", False),
    ("TCGA", "TCGA tumour vs normal, bulk (support)", "Human", "support", False),
]
SUPPORT = [
    ("HGCA", "HGCA all fetal vs adult", "Human", False),
    ("Pikkupeura_LN", "Pikkupeura culture laminin", "Mouse", False),
    ("Pikkupeura_collagen", "Pikkupeura culture collagen", "Mouse", False),
    ("Joanito_sens", "Joanito cohorts with both groups", "Human", True),
    ("Pelka_P", "Pelka tumour vs normal stem/TA (P)", "Human", False),
    ("TCGA_paired", "TCGA paired", "Human", False),
    ("TCGA_GTEx", "TCGA vs GTEx colon", "Human", False),
]


def applicable(g, sp):
    return not ((sp == "Human" and g in NA_HUMAN) or (sp == "Mouse" and g in NA_MOUSE))


def supported(v, q):
    return pd.notna(v) and v >= LFC and pd.notna(q) and q < FDR


def cell(g, sp, v, q):
    if not applicable(g, sp):
        return "NA (no one2one orthologue)"
    if pd.isna(v):
        return "NA (not measured)"
    star = "**" if supported(v, q) else ("*" if pd.notna(q) and q < FDR else "")
    qs = f"FDR {q:.2g}" if pd.notna(q) else "no FDR"
    return f"{v:+.2f} ({qs}){star}"


def calls(g, s, prov):
    get = lambda k: (s.loc[g, f"{k}__log2FC"], s.loc[g, f"{k}__FDR"])
    out = {"L": prov.loc[g, "L_call"]}
    # H
    if g in NA_HUMAN:
        out["H"] = "NE"
    else:
        hv = [get(k) for k in ("HGCA_late", "Hnew2_late", "GaoOriginal_late")]
        meas = [x for x in hv if pd.notna(x[0])]
        if len(meas) < 2:
            out["H"] = "NE"
        else:
            npos = sum(v > 0 for v, _ in meas)
            out["H"] = "pass" if npos >= 2 and any(supported(v, q) for v, q in meas) else "fail"
        out["H_n_measured"], out["H_n_positive"] = len(meas), sum(v > 0 for v, _ in meas)
        out["H_n_supported"] = sum(supported(v, q) for v, q in meas)
    # M
    v, q = get("Mouse_GSE230581")
    out["M"] = "NE" if (g in NA_MOUSE or pd.isna(v)) else ("pass" if supported(v, q) else "fail")
    # C
    if g in NA_HUMAN:
        out["C"] = "NE"
    else:
        cv = [get(k) for k in ("Joanito", "Pelka")]
        failed = any(pd.notna(v) and not supported(v, q) for v, q in cv)
        out["C"] = "fail" if failed else ("NE" if any(pd.isna(v) for v, _ in cv) else "pass")
    # stringent intersection (v1.0): G1 HGCA_late, G2 mouse, G3 Joanito, G4 Pelka
    st = []
    for k in ("HGCA_late", "Mouse_GSE230581", "Joanito", "Pelka"):
        v, q = get(k)
        sp = "Mouse" if k == "Mouse_GSE230581" else "Human"
        st.append("NE" if (not applicable(g, sp) or pd.isna(v)) else ("pass" if supported(v, q) else "fail"))
    out["stringent_dev"] = st[0] == "pass" and st[1] == "pass"
    out["stringent_all"] = all(x == "pass" for x in st)
    out["stringent_calls"] = st
    return out


def label(c):
    if all(c[a] == "pass" for a in "LHMC"):
        return "CORE"
    failed = [a for a in "LHMC" if c[a] == "fail"]
    ne = [a for a in "LHMC" if c[a] == "NE"]
    parts = []
    if failed:
        parts.append("fails " + "+".join(failed))
    if ne:
        parts.append("not evaluable " + "+".join(ne))
    return "not Core: " + "; ".join(parts)


def build(s, prov, restricted):
    ev = [e for e in EVIDENCE if restricted or not e[4]]
    sup = [e for e in SUPPORT if restricted or not e[3]]
    mat, longr, supr = [], [], []
    for g in s.index:
        c = calls(g, s, prov)
        p = prov.loc[g]
        row = {"marker": g, "literature_label": p.literature_label, "curator code": p.curator_code,
               "named studies (n)": int(p.n_named_studies), "L literature": c["L"]}
        for key, lab, sp, axis, _ in ev:
            row[lab] = cell(g, sp, s.loc[g, f"{key}__log2FC"], s.loc[g, f"{key}__FDR"])
            if key == "GaoOriginal_late":
                row["H human developmental"] = c["H"]
            if key == "Mouse_GSE230581":
                row["M mouse in vivo"] = c["M"]
            if key == "Pelka" and restricted:
                row["C CRC replicated"] = c["C"]
            longr.append(dict(marker=g, evidence=lab, axis=axis, log2FC=s.loc[g, f"{key}__log2FC"],
                              PValue=s.loc[g, f"{key}__PValue"], FDR=s.loc[g, f"{key}__FDR"],
                              applicable=applicable(g, sp), supported=supported(s.loc[g, f"{key}__log2FC"],
                                                                                s.loc[g, f"{key}__FDR"])))
        if restricted:
            row["Core label (v2.0)"] = label(c)
            row["Stringent intersection set (v1.0 sensitivity)"] = "yes" if c["stringent_all"] else "no"
        row["Developmental intersection (HGCA & GSE230581)"] = "yes" if c["stringent_dev"] else "no"
        mat.append(row)
        srow = {"marker": g}
        for key, lab, sp, _ in sup:
            srow[lab] = cell(g, sp, s.loc[g, f"{key}__log2FC"], s.loc[g, f"{key}__FDR"])
            longr.append(dict(marker=g, evidence=lab, axis="support", log2FC=s.loc[g, f"{key}__log2FC"],
                              PValue=s.loc[g, f"{key}__PValue"], FDR=s.loc[g, f"{key}__FDR"],
                              applicable=applicable(g, sp),
                              supported=supported(s.loc[g, f"{key}__log2FC"], s.loc[g, f"{key}__FDR"])))
        supr.append(srow)
    return pd.DataFrame(mat), pd.DataFrame(supr), pd.DataFrame(longr)


def rules(restricted):
    r = [
        ("L literature", "mandatory", "curator code contains A (L-A), or contains B with >= 2 named studies (L-B)"),
        ("H human developmental", "mandatory, replicated",
         "HGCA >=9 PCW, H-new2 >=9 PCW, Gao >=9 W: >= 2 measured & fetal-positive AND >= 1 with log2FC >= 0.5 & FDR < 0.05; NE if < 2 measured"),
        ("M mouse in vivo", "mandatory", "GSE230581 log2FC >= 0.5 & FDR < 0.05; NE if no orthologue / not measured"),
        ("C CRC replicated", "mandatory" + ("" if restricted else " (Joanito: restricted)"),
         "Joanito and Pelka both log2FC >= 0.5 & FDR < 0.05; fail if a measured contrast lacks support; NE if otherwise unmeasured"),
        ("Core (v2.0)", "", "L & H & M & C all pass"),
        ("Stringent intersection set (v1.0)", "sensitivity",
         "HGCA >=9 PCW & GSE230581 & Joanito & Pelka each log2FC >= 0.5 & FDR < 0.05 (not measured = NE)"),
        ("Cell notation", "", "log2FC (FDR); ** log2FC >= 0.5 & FDR < 0.05; * FDR < 0.05 below effect threshold; NA = not measured / no orthologue"),
        ("Supportive evidence", "never selecting", "GSE44433, TCGA, Pikkupeura cultures, HGCA all fetal, Joanito sensitivity, Pelka P"),
    ]
    return pd.DataFrame(r, columns=["item", "status", "rule"])


def write_xlsx(path, sheets):
    wb = Workbook()
    wb.remove(wb.active)
    thin = Side(style="thin", color="BFBFBF")
    colours = {"pass": "C6EFCE", "fail": "F8CBAD", "NE": "D9D9D9", "CORE": "00B050", "yes": "C6EFCE",
               "not Core": "F2F2F2"}
    for name, df in sheets.items():
        ws = wb.create_sheet(name)
        for j, c in enumerate(df.columns, 1):
            x = ws.cell(1, j, c)
            axis = c.split(" ")[0] in ("L", "H", "M", "C") or c.startswith(("Core", "Stringent", "Developmental"))
            x.font = Font(name="Arial", bold=True, color="FFFFFF", size=9)
            x.fill = PatternFill("solid", start_color="7030A0" if axis else "1F4E78")
            x.alignment = Alignment(wrap_text=True, vertical="center")
        for i, r in enumerate(df.itertuples(index=False), 2):
            for j, v in enumerate(r, 1):
                v = None if (isinstance(v, float) and np.isnan(v)) else (v.item() if hasattr(v, "item") else v)
                x = ws.cell(i, j, v)
                x.font = Font(name="Arial", size=9, bold=(j == 1))
                x.border = Border(top=thin, bottom=thin, left=thin, right=thin)
                x.alignment = Alignment(wrap_text=True, vertical="top")
                if isinstance(v, str):
                    if v.endswith("**"):
                        x.fill = PatternFill("solid", start_color="E2F0D9")
                    elif v.startswith("NA"):
                        x.fill = PatternFill("solid", start_color="EDEDED")
                    for k, col in colours.items():
                        if v == k or (k in ("CORE", "not Core") and v.startswith(k)):
                            x.fill = PatternFill("solid", start_color=col)
        for j in range(1, len(df.columns) + 1):
            ws.column_dimensions[get_column_letter(j)].width = 12 if j == 1 else (60 if name == "Rules" else 20)
        ws.row_dimensions[1].height = 60
        ws.freeze_panes = "B2"
    wb.save(path)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(MS).set_index("marker")
    prov = pd.read_csv(PROV, sep="\t").set_index("gene")
    for restricted in (True, False):
        mat, sup, longr = build(s, prov, restricted)
        sheets = {"Evidence_matrix_31": mat, "Supportive_31": sup}
        if restricted:
            core = mat[mat["Core label (v2.0)"] == "CORE"]
            sheets["Core_genes_v2"] = core
            sheets["Stringent_set_v1"] = mat[mat["Stringent intersection set (v1.0 sensitivity)"] == "yes"]
        sheets["Developmental_intersection"] = mat[mat["Developmental intersection (HGCA & GSE230581)"] == "yes"]
        sheets["Literature_provenance"] = prov.reset_index()
        sheets["Long_numeric"] = longr
        sheets["Rules"] = rules(restricted)
        sfx = "" if restricted else "_public"
        write_xlsx(OUT / f"Core_oncofetal_gate_statistics{sfx}.xlsx", sheets)
        mat.to_csv(OUT / f"Core_oncofetal_gate_statistics{sfx}.csv", index=False)
        longr.to_csv(OUT / ("Core_oncofetal_gate_long.csv" if restricted else "Core_oncofetal_gate_long_public.csv"),
                     index=False)
        if restricted:
            RESTRICTED.mkdir(parents=True, exist_ok=True)
            for f in ("Core_oncofetal_gate_statistics.xlsx", "Core_oncofetal_gate_statistics.csv",
                      "Core_oncofetal_gate_long.csv"):
                shutil.copy2(OUT / f, RESTRICTED / f)
            pd.set_option("display.width", 250)
            print(mat[["marker", "L literature", "H human developmental", "M mouse in vivo", "C CRC replicated",
                       "Core label (v2.0)"]].to_string(index=False))


if __name__ == "__main__":
    main()
