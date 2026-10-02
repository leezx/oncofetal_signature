#!/usr/bin/env python3
"""31-gene feature audit: why is a candidate NA in a gate contrast?

For every gate dataset x candidate gene, the gene is looked up in the raw
feature/count matrix by current HGNC symbol, historical symbol, Ensembl ID
(human) or mouse symbol / Entrez ID (mouse), using only the samples that the
contrast actually uses. Output fields:
  gene, dataset, canonical_symbol, source_symbol, ensembl_id,
  present_in_raw_features, present_in_raw_counts, total_count,
  n_samples_nonzero, n_samples, passed_expression_filter, present_in_DE_table,
  present_in_marker_summary, NA_reason
NA_reason (only these): not_in_annotation, no_1to1_orthologue, zero_expression,
filtered_low_expression, mapping_failure, not_available_in_source; empty when
the gene has a value. Joanito rows are restricted (output git-ignored).
Usage (repo root): python3 core_oncofetal/scripts/feature_audit_31.py
"""
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/Volumes/Stelligen_SSD/Stelligen/DATA")
BM2 = DATA / "2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step2_human_benchmark_v0.2"
OUT = ROOT / "core_oncofetal/results"
NO_HUMAN = {"LY6A", "REG3B"}
NO_MOUSE = {"SPRR1A"}
MOUSE_HIST = {"Ccn1": ["Cyr61"], "Ccn2": ["Ctgf"]}


def human_keys(gene, alias):
    return [gene] + alias.loc[alias.current_hgnc_symbol == gene, "source_symbol"].tolist()


def counts_lookup(mat, sym_col, id_col, keys, ensg, samples):
    """Return (source_symbol, row) from a count matrix using symbols then Ensembl ID."""
    for k in keys:
        hit = mat[mat[sym_col] == k]
        if len(hit):
            return k, hit[samples].sum(axis=0)
    if ensg and id_col is not None:
        hit = mat[mat[id_col].astype(str).str.split(".").str[0] == ensg]
        if len(hit):
            return f"{ensg} (by Ensembl ID; symbol {hit[sym_col].iloc[0]})", hit[samples].sum(axis=0)
    return None, None


def de_lookup(de, col, keys):
    s = set(de[col].astype(str))
    return next((k for k in keys if k in s), None)


def classify(present, vals, in_de, in_summary, applicable=True, source_has_gene=True):
    if not applicable:
        return "no_1to1_orthologue"
    if not present:
        return "not_in_annotation" if source_has_gene else "not_available_in_source"
    if in_summary:
        return ""
    if vals is not None and float(np.nansum(vals)) == 0:
        return "zero_expression"
    if in_de:
        return "mapping_failure"
    return "filtered_low_expression"


def main():
    cand = pd.read_csv(ROOT / "step2_fetal/config/literature_candidates_31.tsv", sep="\t").gene.tolist()
    alias = pd.read_csv(ROOT / "step3_cancer/config/symbol_aliases.tsv", sep="\t")
    ms = pd.read_csv(ROOT / "marker_summary/results/Literature_31_marker_summary_statistics.csv").set_index("marker")
    hgca_man = pd.read_csv(DATA / "scRNAseq/HGCA_Elmentaite2021/processed/v0.1/HGCA_gene_manifest.csv")
    ensg = {g: hgca_man.loc[hgca_man.gene.isin(human_keys(g, alias)), "gene_ids"].iloc[0]
            if hgca_man.gene.isin(human_keys(g, alias)).any() else "" for g in cand}
    mouse_sym = ms.mouse_gene_used.to_dict()
    rows = []

    def add(dataset, g, src, ens, vals, in_de, summary_key, applicable=True, filt_note=None):
        present = src is not None
        in_sum = pd.notna(ms.loc[g, f"{summary_key}__log2FC"]) if summary_key else False
        rows.append(dict(
            gene=g, dataset=dataset, canonical_symbol=g, source_symbol=src or "", ensembl_id=ens,
            present_in_raw_features=present, present_in_raw_counts=present,
            total_count=float(np.nansum(vals)) if vals is not None else np.nan,
            n_samples_nonzero=int((vals > 0).sum()) if vals is not None else np.nan,
            n_samples=len(vals) if vals is not None else np.nan,
            passed_expression_filter=(bool(in_de) if present else np.nan),
            present_in_DE_table=bool(in_de), present_in_marker_summary=bool(in_sum),
            NA_reason=classify(present, vals, in_de, in_sum, applicable),
            note=filt_note or ""))

    # HGCA >= 9 PCW vs adult (rows = symbols)
    hc = pd.read_csv(DATA / "scRNAseq/HGCA_Elmentaite2021/processed/v0.1/HGCA_primary_pseudobulk_counts.tsv.gz",
                     sep="\t", index_col=0)
    hu = pd.read_csv(ROOT / "step2_benchmark_v2/results/tables/HGCA_stage_resolved_units.csv")
    hs = hu.loc[hu.stage_group.isin(["Fetal_midlate", "Adult"]), "donor_id"].tolist()
    hc = hc.reset_index().rename(columns={"index": "symbol"})
    hc["gene_id"] = hc.symbol.map(dict(zip(hgca_man.gene, hgca_man.gene_ids)))
    hde = pd.read_csv(BM2 / "de/HGCA_late_fetal_vs_adult.csv")
    for g in cand:
        if g in NO_HUMAN:
            add("HGCA >=9 PCW vs adult", g, None, ensg[g], None, False, "HGCA_late", applicable=False)
            continue
        src, v = counts_lookup(hc, "symbol", "gene_id", human_keys(g, alias), ensg[g], hs)
        add("HGCA >=9 PCW vs adult", g, src, ensg[g], v, de_lookup(hde, "symbol", human_keys(g, alias)), "HGCA_late",
            filt_note="edgeR filterByExpr")

    # H-new2 >= 9 PCW (Fawkner) vs Burclaff adult
    nc = pd.read_csv(BM2 / "pseudobulk/Hnew2_Fawkner_fetal_Burclaff_adult_counts.csv.gz")
    nu = pd.read_csv(ROOT / "step2_benchmark_v2/results/tables/Hnew2_stage_resolved_units.csv")
    ns = nu.loc[nu.eligible.astype(str).isin(["True", "TRUE"]) & nu.stage_group.isin(["Fetal_midlate", "Adult"]),
                "pseudobulk"].tolist()
    nde = pd.read_csv(BM2 / "de/Hnew2_late_Fawkner_vs_Burclaff.csv")
    for g in cand:
        if g in NO_HUMAN:
            add("H-new2 >=9 PCW vs adult", g, None, ensg[g], None, False, "Hnew2_late", applicable=False)
            continue
        src, v = counts_lookup(nc, "symbol", "gene_id", human_keys(g, alias), ensg[g], ns)
        add("H-new2 >=9 PCW vs adult", g, src, ensg[g], v, de_lookup(nde, "symbol", human_keys(g, alias)),
            "Hnew2_late", filt_note="edgeR filterByExpr")

    # Gao fetal LI >= 9 W vs GSE103154 adult (effect-only; TPM-scale unit means)
    gde = pd.read_csv(BM2 / "de/GaoOriginal_late_GaoLI_vs_GSE103154.csv")
    gao_genes = pd.read_csv(DATA / "scRNAseq/GSE103239_Gao2018/raw/GSE95630_Digestion_TPM_new.txt.gz", sep="\t",
                            usecols=["Gene"]).Gene
    adult_genes = pd.read_csv(DATA / "scRNAseq/GSE103239_Gao2018/raw/GSE103154_All_Merge_umi_tpm_gene.txt.gz",
                              sep="\t", usecols=["Gene"]).Gene
    gunits = pd.read_csv(ROOT / "step2_benchmark_v2/results/tables/GaoOriginal_late_GaoLI_vs_GSE103154_units.csv")
    for g in cand:
        if g in NO_HUMAN:
            add("Gao LI >=9 W vs GSE103154", g, None, ensg[g], None, False, "GaoOriginal_late", applicable=False)
            continue
        keys = human_keys(g, alias)
        kf = next((k for k in keys if k in set(gao_genes)), None)
        ka = next((k for k in keys if k in set(adult_genes)), None)
        src = None if (kf is None or ka is None) else (kf if kf == ka else f"{kf} (fetal) / {ka} (adult)")
        add("Gao LI >=9 W vs GSE103154", g, src, ensg[g], None, de_lookup(gde, "symbol", keys), "GaoOriginal_late",
            filt_note=("fetal matrix missing symbol" if kf is None else "") + ("; adult matrix missing symbol" if ka is None else ""))

    # Mouse GSE230581 in vivo (ENTREZID, SYMBOL)
    mc = pd.read_csv(DATA / "bulkRNAseq/GSE230581/raw/GSE230581_in_vivo_counts.tsv.gz", sep="\t")
    msamp = [c for c in mc.columns if c not in ("ENTREZID", "SYMBOL")]
    mde = pd.read_csv(ROOT / "step2_fetal/results/tables/mouse_in_vivo_fetal_vs_adult_DEG.csv")
    for g in cand:
        mg = mouse_sym.get(g)
        if g in NO_MOUSE or not isinstance(mg, str):
            add("GSE230581 mouse in vivo", g, None, "", None, False, "Mouse_GSE230581", applicable=False)
            continue
        keys = [mg] + MOUSE_HIST.get(mg, [])
        src, v = counts_lookup(mc, "SYMBOL", None, keys, None, msamp)
        add("GSE230581 mouse in vivo", g, src, "", v, de_lookup(mde, "mouse_gene", keys), "Mouse_GSE230581",
            filt_note="edgeR filterByExpr")

    # Pelka tumour vs normal epithelium; Joanito malignant vs normal (restricted)
    for name, cf, sf, groups, dep, key in [
        ("Pelka tumour vs normal", DATA / "scRNAseq/GSE178341/processed/v0.1/step3_epithelial_pseudobulk/Pelka_epithelial_pseudobulk_counts.csv.gz",
         DATA / "scRNAseq/GSE178341/processed/v0.1/step3_epithelial_pseudobulk/Pelka_pseudobulk_samples.csv",
         ["Tumor", "Normal"], ROOT / "step3_cancer/results/tables/Pelka_tumor_vs_normal_epithelium_DEG.csv", "Pelka"),
        ("Joanito malignant vs normal", DATA / "scRNAseq/Joanito2022_syn26844071/processed/v0.1/step3_epithelial_pseudobulk/Joanito_epithelial_pseudobulk_counts.csv.gz",
         DATA / "scRNAseq/Joanito2022_syn26844071/processed/v0.1/step3_epithelial_pseudobulk/Joanito_pseudobulk_samples.csv",
         ["Malignant", "Normal"], ROOT / "step3_cancer/results/tables/Joanito_Malignant_vs_Normal.csv", "Joanito")]:
        cc = pd.read_csv(cf)
        sm = pd.read_csv(sf)
        ss = sm.loc[sm.eligible.astype(str).isin(["True", "TRUE"]) & sm.group.isin(groups), "pseudobulk"].tolist()
        de = pd.read_csv(dep)
        for g in cand:
            if g in NO_HUMAN:
                add(name, g, None, ensg[g], None, False, key, applicable=False)
                continue
            src, v = counts_lookup(cc, "symbol", "gene_id", human_keys(g, alias), ensg[g], ss)
            add(name, g, src, ensg[g], v, de_lookup(de, "symbol", human_keys(g, alias)), key,
                filt_note="edgeR filterByExpr")
    a = pd.DataFrame(rows)
    a.to_csv(OUT / "Feature_audit_31_restricted.csv", index=False)
    a[a.dataset != "Joanito malignant vs normal"].to_csv(OUT / "Feature_audit_31_public.csv", index=False)
    pd.set_option("display.width", 220)
    na = a[a.NA_reason != ""]
    print(na[["gene", "dataset", "source_symbol", "total_count", "n_samples_nonzero", "n_samples",
              "present_in_DE_table", "NA_reason", "note"]].to_string(index=False))
    print(a.groupby(["dataset", "NA_reason"]).size())


if __name__ == "__main__":
    main()
