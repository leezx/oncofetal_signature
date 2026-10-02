#!/usr/bin/env python3
"""31-marker x contrast benchmark matrix and frozen qualification verdicts.

Rule 1: TNFRSF12A log2FC > 0 and P < 0.05 (Pikkupeura retains only FDR, which
        is used in its place and is the stricter quantity).
Rule 2: one-sided exact binomial test that > 50% of measured markers have
        log2FC > 0, P < 0.05. LY6A/REG3B are not counted in human contrasts.
A contrast qualifies only if both hold. No signature is built here.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import binomtest

PRIORITY = ["TACSTD2", "CLU", "ANXA1", "TNFRSF12A"]

# key, label, species, source description
CONTRASTS = [
    ("HGCA", "HGCA fetal vs adult epithelium (original H1)", "Human", "existing"),
    ("GaoOriginal", "Gao fetal LI vs GSE103154 adult (recomputed)", "Human", "new module 04"),
    ("Senger", "Senger fetal vs adult enterospheres (GSE101531)", "Human", "existing benchmark"),
    ("Hnew1", "H-new1 Fawkner Visium fetal vs adult epithelial spots", "Human", "new"),
    ("Hnew1_colon", "H-new1 colon-only sensitivity", "Human", "new"),
    ("Hnew2", "H-new2 Fawkner fetal scRNA vs Burclaff adult", "Human", "new"),
    ("Hbulk1", "H-bulk1 Roadmap fetal SI vs HPA adult SI+duodenum (Finkbeiner design)", "Human", "addendum v1.1"),
    ("Hbulk2", "H-bulk2 Roadmap fetal SI vs HPA adult duodenum (Senger primary-tissue subset)", "Human", "addendum v1.1"),
    ("Hnew3", "H-new3 Gao fetal SI+LI vs Wang adult", "Human", "new"),
    ("Hnew3_debug", "H-new3-debug Gao fetal LI vs Wang colon+rectum", "Human", "new"),
    ("Pikkupeura_LN", "Pikkupeura fetal vs adult culture, laminin (GSE160449)", "Mouse", "existing benchmark"),
    ("Pikkupeura_collagen", "Pikkupeura fetal vs adult culture, collagen (GSE160449)", "Mouse", "existing benchmark"),
]


def by_symbol(df, sym, lfc, p, fdr, keys):
    idx = df.drop_duplicates(sym).set_index(sym)
    out = {}
    for g, ks in keys.items():
        hit = next((k for k in ks if k in idx.index), None)
        out[g] = (idx.loc[hit, lfc], idx.loc[hit, p] if p else np.nan, idx.loc[hit, fdr] if fdr else np.nan) \
            if hit is not None else (np.nan, np.nan, np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    repo, de, out = pathlib.Path(a.repo), pathlib.Path(a.de_dir), pathlib.Path(a.out_dir)
    cand = pd.read_csv(repo / "step2_fetal/config/literature_candidates_31.tsv", sep="\t")
    alias = pd.read_csv(repo / "step3_cancer/config/symbol_aliases.tsv", sep="\t")
    keys = {g: [g] + alias.loc[alias.current_hgnc_symbol == g, "source_symbol"].tolist() for g in cand.gene}
    hgca = pd.read_csv(repo / "step2_fetal/results/tables/human_HGCA_fetal_vs_adult_DEG.csv")
    bm = pd.read_csv(repo / "step2_dataset_benchmark/results/tables/Literature_31_four_dataset_benchmark.csv")
    bm_keys = {g: [g] for g in cand.gene}

    vals = {
        "HGCA": by_symbol(hgca, "gene", "HGCA_log2FC", "HGCA_PValue", "HGCA_FDR", keys),
        "GaoOriginal": by_symbol(pd.read_csv(de / "GaoOriginal_GaoLI_vs_GSE103154.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Senger": by_symbol(bm, "gene", "Human_Senger_log2FC", "Human_Senger_PValue", "Human_Senger_FDR", bm_keys),
        "Hnew1": by_symbol(pd.read_csv(de / "Hnew1_Fawkner_spatial_units.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hnew1_colon": by_symbol(pd.read_csv(de / "Hnew1_Fawkner_spatial_units_colon_only.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hnew2": by_symbol(pd.read_csv(de / "Hnew2_Fawkner_Burclaff.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hbulk1": by_symbol(pd.read_csv(de / "Hbulk1_Roadmap_vs_HPA_SI.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hbulk2": by_symbol(pd.read_csv(de / "Hbulk2_Roadmap_vs_HPA_duodenum.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hnew3": by_symbol(pd.read_csv(de / "Hnew3_GaoSILI_vs_Wang.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Hnew3_debug": by_symbol(pd.read_csv(de / "Hnew3debug_GaoLI_vs_WangColonRectum.csv"), "symbol", "log2FC", "PValue", "FDR", keys),
        "Pikkupeura_LN": by_symbol(bm, "gene", "Mouse_Pikkupeura_LN_log2FC", None, "Mouse_Pikkupeura_LN_FDR", bm_keys),
        "Pikkupeura_collagen": by_symbol(bm, "gene", "Mouse_Pikkupeura_collagen_log2FC", None, "Mouse_Pikkupeura_collagen_FDR", bm_keys),
    }
    order = PRIORITY + [g for g in cand.gene if g not in PRIORITY]
    mat = pd.DataFrame({"marker": order, "priority_control": [g in PRIORITY for g in order]})
    for k, *_ in CONTRASTS:
        for i, s in enumerate(["log2FC", "PValue", "FDR"]):
            mat[f"{k}__{s}"] = [vals[k][g][i] for g in order]
    out.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out / "Benchmark_v2_31_marker_matrix.csv", index=False)

    rows = []
    for k, label, species, src in CONTRASTS:
        lfc = mat.set_index("marker")[f"{k}__log2FC"]
        if species == "Human":
            lfc = lfc.drop(["LY6A", "REG3B"])
        meas = lfc.dropna()
        n_pos = int((meas > 0).sum())
        binom_p = binomtest(n_pos, len(meas), 0.5, alternative="greater").pvalue if len(meas) else np.nan
        tn = vals[k]["TNFRSF12A"]
        support = tn[1] if not np.isnan(tn[1]) else tn[2]
        r1 = bool(tn[0] > 0 and support < 0.05) if not np.isnan(tn[0]) else False
        r2 = bool(binom_p < 0.05)
        rows.append(dict(contrast=k, description=label, species=species, source=src,
                         TNFRSF12A_log2FC=tn[0], TNFRSF12A_P=tn[1], TNFRSF12A_FDR=tn[2],
                         rule1_TNFRSF12A_fetal_high_supported=r1,
                         markers_measured=len(meas), markers_fetal_positive=n_pos,
                         fraction_fetal_positive=round(n_pos / len(meas), 3) if len(meas) else np.nan,
                         binomial_P_one_sided=binom_p, rule2_panel_bias=r2,
                         **{f"{g}_log2FC": vals[k][g][0] for g in ["TACSTD2", "CLU", "ANXA1"]},
                         qualifies=r1 and r2))
    q = pd.DataFrame(rows)
    olfm4 = {}
    srcs = {"HGCA": (hgca, "gene", "HGCA_log2FC"), "Senger": None, "Pikkupeura_LN": None, "Pikkupeura_collagen": None}
    for k, *_ in CONTRASTS:
        f = {"GaoOriginal": "GaoOriginal_GaoLI_vs_GSE103154.csv", "Hnew1": "Hnew1_Fawkner_spatial_units.csv",
             "Hnew1_colon": "Hnew1_Fawkner_spatial_units_colon_only.csv", "Hnew2": "Hnew2_Fawkner_Burclaff.csv",
             "Hbulk1": "Hbulk1_Roadmap_vs_HPA_SI.csv", "Hbulk2": "Hbulk2_Roadmap_vs_HPA_duodenum.csv",
             "Hnew3": "Hnew3_GaoSILI_vs_Wang.csv", "Hnew3_debug": "Hnew3debug_GaoLI_vs_WangColonRectum.csv"}.get(k)
        if f:
            d = pd.read_csv(de / f).drop_duplicates("symbol").set_index("symbol")
            olfm4[k] = d.loc["OLFM4", "log2FC"] if "OLFM4" in d.index else np.nan
        elif k == "HGCA":
            d = hgca.drop_duplicates("gene").set_index("gene")
            olfm4[k] = d.loc["OLFM4", "HGCA_log2FC"] if "OLFM4" in d.index else np.nan
        else:
            olfm4[k] = np.nan
    q["OLFM4_log2FC_sanity"] = q.contrast.map(olfm4)
    q["OLFM4_fetal_lt_adult"] = q.OLFM4_log2FC_sanity < 0
    # Review decisions recorded after the frozen verdicts (never change them).
    dec = repo / "step2_benchmark_v2/config/review_decisions.tsv"
    if dec.exists():
        q = q.merge(pd.read_csv(dec, sep="\t"), on="contrast", how="left")
    q.to_csv(out / "Benchmark_v2_qualification.csv", index=False)
    pd.set_option("display.width", 200)
    print(q[["contrast", "TNFRSF12A_log2FC", "TNFRSF12A_P", "rule1_TNFRSF12A_fetal_high_supported",
             "markers_fetal_positive", "markers_measured", "binomial_P_one_sided", "rule2_panel_bias",
             "qualifies"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
