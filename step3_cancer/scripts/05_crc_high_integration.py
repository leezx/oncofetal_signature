#!/usr/bin/env python3
"""Dataset admission QA, gene-level CRC-high evidence table, and the
31-candidate descriptive audit (rules frozen in docs/ANALYSIS_PLAN.md).

Joanito S1 is optional at run time: if its DE tables are absent, S1 columns
are left empty, every gene that would need S1 is labelled `S1_pending`, and
no CRC_high call is made.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd

UP = ["CDH3", "CLDN1", "FOXQ1", "KRT23", "LGR5", "ASCL2", "MYC", "TESC"]
DOWN = ["CA1", "CA2", "CA4", "GUCA2A", "GUCA2B", "CLCA4", "AQP8", "SLC26A3", "MS4A12", "CEACAM7"]


def admission(name, df, lfc_col, gated=True):
    rows = []
    for panel, genes, want in (("tumour_up", UP, 1), ("tumour_down", DOWN, -1)):
        sub = df[df.symbol.isin(genes)].drop_duplicates("symbol")
        measured = sub[lfc_col].notna()
        correct = int((np.sign(sub.loc[measured, lfc_col]) == want).sum())
        n = int(measured.sum())
        rows.append(dict(contrast=name, panel=panel, n_panel=len(genes), n_measured=n,
                         n_expected_direction=correct,
                         fraction=correct / n if n else np.nan,
                         wrong_direction=";".join(sub.loc[measured & (np.sign(sub[lfc_col]) != want), "symbol"]),
                         not_measured=";".join(sorted(set(genes) - set(sub.loc[measured, "symbol"])))))
    out = pd.DataFrame(rows)
    out["panel_pass"] = out.fraction >= 0.8
    # Admission applies only to gated contrasts (T1, S1, S2). The
    # proliferation control compares against normal stem/TA cells, for which
    # the differentiation-marker panel has no expected direction.
    out["contrast_admitted"] = bool(out.panel_pass.all()) if gated else "not_applicable"
    return out


def prefixed(path, prefix):
    df = pd.read_csv(path)
    df = df.rename(columns={c: f"{prefix}_{c}" for c in ["log2FC", "logCPM", "PValue", "FDR"]})
    return df.drop(columns=[f"{prefix}_logCPM", f"{prefix}_PValue"]).drop_duplicates("gene_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--aliases", required=True, help="current HGNC -> source annotation symbol map")
    ap.add_argument("--orthology", required=True, help="non-human candidate orthology notes")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-log2fc", type=float, default=0.5)
    ap.add_argument("--max-fdr", type=float, default=0.05)
    a = ap.parse_args()
    de = pathlib.Path(a.de_dir)
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    tcga = pd.read_csv(de / "TCGA_tumor_vs_normal_DEG.csv")
    tcga = tcga[tcga.T1_tested].drop_duplicates("gene_id")
    s2 = prefixed(de / "Pelka_Tumor_vs_Normal.csv", "Pelka")
    pelka_ctrl = prefixed(de / "Pelka_Tumor_vs_Normal_prolif.csv", "PelkaProlifCtrl")

    joanito_f = de / "Joanito_Malignant_vs_Normal.csv"
    sens_f = de / "Joanito_Malignant_vs_Normal_sensitivity_paired_cohorts.csv"
    have_s1 = joanito_f.exists()

    qa = [admission("T1_TCGA", tcga.rename(columns={"T1_log2FC": "lfc"}), "lfc"),
          admission("S2_Pelka", s2, "Pelka_log2FC"),
          admission("Ctrl_Pelka_vs_normal_prolif", pelka_ctrl, "PelkaProlifCtrl_log2FC", gated=False)]
    if have_s1:
        s1 = prefixed(joanito_f, "Joanito")
        qa.append(admission("S1_Joanito", s1, "Joanito_log2FC"))
    qa = pd.concat(qa, ignore_index=True)
    qa.to_csv(out / "Dataset_admission_QA.csv", index=False)
    # Public copy without Joanito rows (Synapse terms forbid disclosing
    # Joanito-derived material; this repository is public).
    qa[~qa.contrast.str.contains("Joanito")].to_csv(out / "Dataset_admission_QA_TCGA_Pelka.csv", index=False)

    ev = tcga[["gene_id", "symbol", "gene_type", "T1_log2FC", "T1_FDR", "T1_pass",
               "Paired_log2FC", "Paired_FDR", "Paired_concordant",
               "GTEx_log2FC", "GTEx_FDR", "GTEx_concordant"]]
    ev = ev.merge(s2.drop(columns="symbol"), on="gene_id", how="outer")
    ev = ev.merge(pelka_ctrl.drop(columns="symbol"), on="gene_id", how="outer")
    sym = pd.concat([tcga[["gene_id", "symbol"]], s2[["gene_id", "symbol"]]]).drop_duplicates("gene_id")
    ev = ev.drop(columns="symbol").merge(sym, on="gene_id", how="left")
    ev["match_key"] = "ensembl_gene_id"

    ev["S2_pass"] = (ev.Pelka_log2FC > 0) & (ev.Pelka_FDR < a.max_fdr)
    # Independent progenitor-specificity check (plan v3): Pelka tumour
    # epithelium vs Pelka normal stem/TA (cE01-03), direction only. It is a
    # separate biological check, not part of S1 and not a substitute for
    # Joanito labels; it keeps cycling/crypt-progenitor genes out of CRC_high.
    ev["Progenitor_check"] = ev.PelkaProlifCtrl_log2FC > 0
    if have_s1:
        # Joanito releases symbols only (GRCh38 Ensembl 93): match by symbol
        # against the source symbol of the Ensembl-keyed table.
        s1 = s1.drop(columns="gene_id").drop_duplicates("symbol")
        ev = ev.merge(s1, on="symbol", how="left")
        ev.loc[ev.Joanito_log2FC.notna(), "match_key"] = "ensembl_gene_id; Joanito by symbol"
        if sens_f.exists():
            sens = prefixed(sens_f, "JoanitoSens").drop(columns="gene_id").drop_duplicates("symbol")
            ev = ev.merge(sens, on="symbol", how="left")
            ev["S1_sensitivity_pass"] = (ev.JoanitoSens_log2FC >= a.min_log2fc) & (ev.JoanitoSens_FDR < a.max_fdr)
        ev["S1_pass"] = (ev.Joanito_log2FC >= a.min_log2fc) & (ev.Joanito_FDR < a.max_fdr)
    else:
        for c in ["Joanito_log2FC", "Joanito_FDR", "S1_pass"]:
            ev[c] = np.nan

    # Plan v2: malignant-epithelial evidence defines CRC-high; TCGA bulk is
    # orthogonal population-level support with no veto, because
    # stromal/immune composition can dilute or invert epithelial
    # re-expression in bulk tumour RNA.
    ev["bulk_support"] = ev.T1_pass.fillna(False).astype(bool)
    t1 = ev.bulk_support
    s2p = ev.S2_pass.fillna(False).astype(bool)
    chk = ev.Progenitor_check.fillna(False).astype(bool)
    if have_s1:
        s1p = ev.S1_pass.fillna(False).astype(bool)
        lab = np.select(
            [s1p & s2p & chk, s1p & s2p & ~chk, s1p & ~s2p, t1 & ~s1p],
            ["CRC_high", "Progenitor_associated_reject", "CRC_high_unreplicated",
             "Tumour_level_only"],
            default="Not_CRC_high")
    else:
        lab = np.where(s2p, "S2_pass_S1_pending", "S2_fail_S1_pending")
    ev["final_label"] = lab
    ev["CRC_high"] = ev.final_label == "CRC_high"

    # Report current HGNC symbols; keep the annotation's own symbol.
    alias = pd.read_csv(a.aliases, sep="\t")
    to_current = dict(zip(alias.source_symbol, alias.current_hgnc_symbol))
    ev = ev.rename(columns={"symbol": "source_symbol"})
    ev["hgnc_symbol"] = ev.source_symbol.map(to_current).fillna(ev.source_symbol)
    ev = ev.sort_values(["CRC_high", "S1_pass" if have_s1 else "S2_pass", "T1_log2FC"],
                        ascending=[False, False, False])
    lead = ["gene_id", "hgnc_symbol", "source_symbol", "match_key", "final_label", "CRC_high",
            "bulk_support"]
    ev = ev[lead + [c for c in ev.columns if c not in lead]]
    ev.to_csv(out / "CRC_high_evidence.csv", index=False)

    cand = pd.read_csv(a.candidates, sep="\t")
    orth = pd.read_csv(a.orthology, sep="\t").set_index("candidate")
    audit = cand.merge(ev, left_on="gene", right_on="hgnc_symbol", how="left")
    absent = audit.hgnc_symbol.isna()
    audit["annotation_status"] = np.where(
        audit.gene.isin(orth.index), audit.gene.map(orth.audit_status),
        np.where(absent, "symbol_absent_from_annotation", "present"))
    audit["T1_status"] = np.where(absent, audit.annotation_status,
                                  np.where(audit.T1_log2FC.isna(), "not_tested", "tested"))
    audit["Pelka_status"] = np.where(absent, audit.annotation_status,
                                     np.where(audit.Pelka_log2FC.isna(), "filtered_low_expression", "tested"))
    if have_s1:
        audit["Joanito_status"] = np.where(absent, audit.annotation_status,
                                           np.where(audit.Joanito_log2FC.isna(), "filtered_low_expression", "tested"))
    audit.to_csv(out / "Literature_31_CRC_axis_audit.csv", index=False)

    print(qa[["contrast", "panel", "n_measured", "n_expected_direction", "fraction", "contrast_admitted"]]
          .to_string(index=False))
    print(ev.final_label.value_counts().to_string())


if __name__ == "__main__":
    main()
