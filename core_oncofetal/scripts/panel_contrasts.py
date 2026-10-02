#!/usr/bin/env python3
"""Re-run the edgeR gate contrasts with the panel-gene filter exemption
(CIOC data-QC amendment) and extract the 31-gene panel values.

Contrasts (original configuration reproduced):
  HGCA_late        HGCA fetal >=9 PCW vs adult        design_robust, ~ group, min 3
  Hnew2_late       Fawkner >=9 PCW vs Burclaff adult  generic, ~ group, min 2
  Mouse_GSE230581  E16.5 epithelium vs adult crypt     design_robust, ~ group, min 3
  Pelka            tumour vs normal epithelium         generic, ~ group, min 3
  Joanito          malignant vs normal epithelium      generic, ~ cohort + group, min 3 (restricted)
Gao (effect-only) has no expression filter; its panel values are taken as is,
with all-zero genes labelled zero_expression.
Outputs: DATA work dir (full DE), core_oncofetal/results/Panel_gate_values_{restricted,public}.csv
and a reproduction check against the original DE tables.
Usage (repo root): python3 core_oncofetal/scripts/panel_contrasts.py
"""
import pathlib
import subprocess

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA")
BM2 = DATA / "2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step2_human_benchmark_v0.2"
WORK = DATA / "2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_core_oncofetal_panel_v0.1"
OUT = ROOT / "core_oncofetal/results"
R = ROOT / "core_oncofetal/scripts/panel_edger.R"
MOUSE_HIST = {"Ccn1": ["Cyr61"], "Ccn2": ["Ctgf"]}
IMPLAUSIBLE_ZERO = {("GaoOriginal_late", "MIF"): "MIF zero in every Gao and GSE103154 unit (ubiquitous gene): likely source quantification artefact"}


def run(name, counts, samples, case, ref, covar, min_units, panel, mode):
    WORK.mkdir(parents=True, exist_ok=True)
    cf, sf, pf, of = (WORK / f"{name}_counts.csv.gz", WORK / f"{name}_samples.csv",
                      WORK / f"{name}_panel.txt", WORK / f"{name}_panel_DE.csv")
    counts.to_csv(cf, index=False)
    samples.to_csv(sf, index=False)
    pf.write_text("\n".join(panel))
    subprocess.run(["Rscript", str(R), str(cf), str(sf), case, ref, covar, str(of), str(min_units), str(pf), mode],
                   check=True)
    return pd.read_csv(of)


def main():
    cand = pd.read_csv(ROOT / "step2_fetal/config/literature_candidates_31.tsv", sep="\t").gene.tolist()
    alias = pd.read_csv(ROOT / "step3_cancer/config/symbol_aliases.tsv", sep="\t")
    hkeys = {g: [g] + alias.loc[alias.current_hgnc_symbol == g, "source_symbol"].tolist() for g in cand}
    hpanel = sorted({k for ks in hkeys.values() for k in ks})
    ms = pd.read_csv(ROOT / "marker_summary/results/Literature_31_marker_summary_statistics.csv").set_index("marker")
    mkeys = {g: [m] + MOUSE_HIST.get(m, []) for g, m in ms.mouse_gene_used.items() if isinstance(m, str)}
    mpanel = sorted({k for ks in mkeys.values() for k in ks})
    res = {}

    # HGCA >= 9 PCW vs adult
    man = pd.read_csv(DATA / "scRNAseq/HGCA_Elmentaite2021/processed/v0.1/HGCA_gene_manifest.csv")
    hc = pd.read_csv(DATA / "scRNAseq/HGCA_Elmentaite2021/processed/v0.1/HGCA_primary_pseudobulk_counts.tsv.gz",
                     sep="\t", index_col=0)
    hc.insert(0, "symbol", hc.index)
    idmap = dict(zip(man.gene, man.gene_ids))
    hc.insert(0, "gene_id", [idmap.get(x, x) for x in hc.index])
    hu = pd.read_csv(ROOT / "step2_benchmark_v2/results/tables/HGCA_stage_resolved_units.csv")
    hs = pd.DataFrame({"pseudobulk": hu.donor_id, "group": hu.stage_group, "eligible": True})
    res["HGCA_late"] = run("HGCA_late", hc.reset_index(drop=True), hs, "Fetal_midlate", "Adult", "none", 3, hpanel,
                           "design_robust")

    # H-new2 >= 9 PCW
    nc = pd.read_csv(BM2 / "pseudobulk/Hnew2_Fawkner_fetal_Burclaff_adult_counts.csv.gz")
    nu = pd.read_csv(ROOT / "step2_benchmark_v2/results/tables/Hnew2_stage_resolved_units.csv")
    ns = pd.DataFrame({"pseudobulk": nu.pseudobulk, "group": nu.stage_group, "eligible": nu.eligible})
    res["Hnew2_late"] = run("Hnew2_late", nc, ns, "Fetal_midlate", "Adult", "none", 2, hpanel, "generic")

    # Mouse GSE230581 in vivo
    mc = pd.read_csv(DATA / "bulkRNAseq/GSE230581/raw/GSE230581_in_vivo_counts.tsv.gz", sep="\t")
    mc = mc.rename(columns={"ENTREZID": "gene_id", "SYMBOL": "symbol"})
    mc["symbol"] = mc.symbol.fillna(mc.gene_id.astype(str))
    msamp = [c for c in mc.columns if c not in ("gene_id", "symbol")]
    msd = pd.DataFrame({"pseudobulk": msamp, "group": ["Fetal" if c.startswith("E16") else "Adult" for c in msamp],
                        "eligible": True})
    res["Mouse_GSE230581"] = run("Mouse_GSE230581", mc, msd, "Fetal", "Adult", "none", 3, mpanel, "design_robust")

    # Pelka and Joanito
    pb = DATA / "scRNAseq/GSE178341/processed/v0.1/step3_epithelial_pseudobulk"
    res["Pelka"] = run("Pelka", pd.read_csv(pb / "Pelka_epithelial_pseudobulk_counts.csv.gz"),
                       pd.read_csv(pb / "Pelka_pseudobulk_samples.csv"), "Tumor", "Normal", "none", 3, hpanel, "generic")
    jb = DATA / "scRNAseq/Joanito2022_syn26844071/processed/v0.1/step3_epithelial_pseudobulk"
    res["Joanito"] = run("Joanito", pd.read_csv(jb / "Joanito_epithelial_pseudobulk_counts.csv.gz"),
                         pd.read_csv(jb / "Joanito_pseudobulk_samples.csv"), "Malignant", "Normal", "cohort", 3,
                         hpanel, "generic")

    # Gao effect-only: values as published in the DE table; all-zero genes -> zero_expression
    gde = pd.read_csv(BM2 / "de/GaoOriginal_late_GaoLI_vs_GSE103154.csv")

    rows, repro = [], []
    original = {"HGCA_late": (BM2 / "de/HGCA_late_fetal_vs_adult.csv", "symbol", "log2FC"),
                "Hnew2_late": (BM2 / "de/Hnew2_late_Fawkner_vs_Burclaff.csv", "symbol", "log2FC"),
                "Mouse_GSE230581": (ROOT / "step2_fetal/results/tables/mouse_in_vivo_fetal_vs_adult_DEG.csv",
                                    "mouse_gene", "Mouse_log2FC"),
                "Pelka": (ROOT / "step3_cancer/results/tables/Pelka_tumor_vs_normal_epithelium_DEG.csv", "symbol", "log2FC"),
                "Joanito": (ROOT / "step3_cancer/results/tables/Joanito_Malignant_vs_Normal.csv", "symbol", "log2FC")}
    for key, de in res.items():
        mouse = key == "Mouse_GSE230581"
        od = pd.read_csv(original[key][0]).drop_duplicates(original[key][1]).set_index(original[key][1])
        for g in cand:
            keys = mkeys.get(g, []) if mouse else hkeys[g]
            if (mouse and g == "SPRR1A") or (not mouse and g in ("LY6A", "REG3B")):
                rows.append(dict(contrast=key, gene=g, NA_reason="no_1to1_orthologue"))
                continue
            hit = de[de.symbol.isin(keys)]
            if hit.empty:
                rows.append(dict(contrast=key, gene=g, NA_reason="not_in_annotation"))
                continue
            r = hit.iloc[0]
            reason = "zero_expression" if pd.isna(r.log2FC) else ""
            low = (r.mean_CPM_case < 1 and r.mean_CPM_ref < 1) if pd.notna(r.log2FC) else np.nan
            rows.append(dict(contrast=key, gene=g, source_symbol=r.symbol, log2FC=r.log2FC, PValue=r.PValue,
                             FDR=r.FDR, mean_CPM_case=r.mean_CPM_case, mean_CPM_ref=r.mean_CPM_ref,
                             forced_by_panel_exemption=bool(r.forced_by_panel_exemption),
                             low_count_flag=low, NA_reason=reason))
            ok = next((k for k in keys if k in od.index), None)
            if ok is not None and pd.notna(r.log2FC):
                repro.append(dict(contrast=key, gene=g, original_log2FC=od.loc[ok, original[key][2]],
                                  panel_log2FC=r.log2FC))
    for g in cand:
        if g in ("LY6A", "REG3B"):
            rows.append(dict(contrast="GaoOriginal_late", gene=g, NA_reason="no_1to1_orthologue"))
            continue
        hit = gde[gde.symbol.isin(hkeys[g])]
        if hit.empty:
            note = IMPLAUSIBLE_ZERO.get(("GaoOriginal_late", g))
            rows.append(dict(contrast="GaoOriginal_late", gene=g,
                             NA_reason="not_available_in_source" if note else "zero_expression", note=note or ""))
        else:
            r = hit.iloc[0]
            rows.append(dict(contrast="GaoOriginal_late", gene=g, source_symbol=r.symbol, log2FC=r.log2FC,
                             PValue=r.PValue, FDR=r.FDR, forced_by_panel_exemption=False, NA_reason=""))
    v = pd.DataFrame(rows)
    v.to_csv(OUT / "Panel_gate_values_restricted.csv", index=False)
    v[v.contrast != "Joanito"].to_csv(OUT / "Panel_gate_values_public.csv", index=False)
    rp = pd.DataFrame(repro)
    rp["abs_diff"] = (rp.original_log2FC - rp.panel_log2FC).abs()
    rp.to_csv(WORK / "reproduction_check.csv", index=False)
    print("reproduction (genes in original DE): max |diff| by contrast")
    print(rp.groupby("contrast").abs_diff.max().round(4).to_string())
    pd.set_option("display.width", 200)
    print(v[(v.forced_by_panel_exemption == True) | (v.NA_reason != "")].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
