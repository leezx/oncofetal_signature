#!/usr/bin/env python3
"""H-new1: Fawkner-Corbett 2021 Visium (GSE158328) epithelial-spot pseudobulks.

Per section: in-tissue spots; epithelial score = fraction of spot UMIs from a
fixed epithelial gene set (no candidate oncofetal marker). Epithelial-rich
spots = top quartile of the score within the section with EPCAM > 0.
One raw-count pseudobulk per section. Rules frozen in docs/QUALIFICATION_PLAN.md.
"""
import argparse
import gzip
import pathlib
import re
import tarfile

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

EPI = ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT20", "CLDN7", "VIL1"]
STROMA = ["VIM", "COL1A1", "COL3A1", "DCN", "ACTA2", "PTPRC"]


def section_meta(brief):
    meta, cur = {}, None
    for line in open(brief):
        if line.startswith("^SAMPLE"):
            cur = line.split("=")[1].strip()
            meta[cur] = {}
        elif line.startswith("!Sample_title"):
            meta[cur]["section"] = line.split("=")[1].strip().split()[0]
        elif line.startswith("!Sample_characteristics_ch1"):
            k, v = line.split("=", 1)[1].strip().split(":", 1)
            meta[cur][k.strip()] = v.strip()
    m = pd.DataFrame.from_dict(meta, orient="index").rename_axis("gsm").reset_index()
    m["group"] = np.where(m.tissue.str.startswith("Adult"), "Adult", "Fetal")
    m["region"] = np.where(m.tissue.str.contains("Colon", case=False), "colon", "small_intestine")
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-spots", type=int, default=20)
    a = ap.parse_args()
    raw, work, out = map(pathlib.Path, (a.raw_dir, a.work_dir, a.out_dir))
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    meta = section_meta(raw / "GSE158328_samples_brief.txt")

    outer = tarfile.open(raw / "GSE158328_RAW.tar")
    pbs, rows, genes = {}, [], None
    for member in outer.getmembers():
        gsm = re.match(r"(GSM\d+)_", member.name).group(1)
        sec = meta.set_index("gsm").loc[gsm]
        inner = tarfile.open(fileobj=outer.extractfile(member), mode="r:gz")
        get = lambda suffix: inner.extractfile(next(m for m in inner.getmembers() if m.name.endswith(suffix)))
        feats = pd.read_csv(gzip.open(get("raw_feature_bc_matrix/features.tsv.gz")), sep="\t", header=None)
        bcs = pd.read_csv(gzip.open(get("raw_feature_bc_matrix/barcodes.tsv.gz")), header=None)[0]
        mat = sp.csc_matrix(scipy.io.mmread(gzip.open(get("raw_feature_bc_matrix/matrix.mtx.gz"))))
        pos = pd.read_csv(get("spatial/tissue_positions_list.csv"), header=None)
        in_tissue = set(pos.loc[pos[1] == 1, 0])
        keep = bcs.isin(in_tissue).to_numpy()
        mat = mat[:, keep]
        sym = feats[1].to_numpy()
        if genes is None:
            genes = feats[[0, 1]].rename(columns={0: "gene_id", 1: "symbol"})
        elif not (genes.gene_id.to_numpy() == feats[0].to_numpy()).all():
            raise SystemExit("Feature order differs between sections")
        total = np.asarray(mat.sum(axis=0)).ravel()
        epi_idx = np.where(np.isin(sym, EPI))[0]
        str_idx = np.where(np.isin(sym, STROMA))[0]
        epcam = np.asarray(mat[np.where(sym == "EPCAM")[0]].sum(axis=0)).ravel()
        score = np.asarray(mat[epi_idx].sum(axis=0)).ravel() / np.maximum(total, 1)
        q75 = np.quantile(score, 0.75)
        sel = (score >= q75) & (epcam > 0) & (total > 0)
        pb = np.asarray(mat[:, sel].sum(axis=1)).ravel()
        stroma_frac = mat[str_idx][:, sel].sum() / max(pb.sum(), 1)
        unit = f"{sec.section}__{sec.group}"
        pbs[unit] = pb.astype(np.int64)
        rows.append(dict(pseudobulk=unit, patient=sec.section, group=sec.group, region=sec.region,
                         tissue=sec.tissue, pcw=sec.pcw, gsm=gsm, n_tissue_spots=int(keep.sum()),
                         n_cells=int(sel.sum()), epi_score_q75=round(float(q75), 5),
                         epithelial_umi_fraction=round(float(mat[epi_idx][:, sel].sum() / max(pb.sum(), 1)), 5),
                         stromal_immune_umi_fraction=round(float(stroma_frac), 5),
                         library_size=int(pb.sum())))
    counts = pd.DataFrame(pbs)
    counts.insert(0, "symbol", genes.symbol.to_numpy())
    counts.index = genes.gene_id.to_numpy()
    counts.index.name = "gene_id"
    counts.to_csv(work / "Fawkner_spatial_epithelial_pseudobulk_counts.csv.gz")
    units = pd.DataFrame(rows)
    units["eligible"] = units.n_cells >= a.min_spots  # n_cells = selected spots
    units.to_csv(work / "Fawkner_spatial_units.csv", index=False)
    colon = units.assign(eligible=units.eligible & (units.region == "colon"))
    colon.to_csv(work / "Fawkner_spatial_units_colon_only.csv", index=False)
    units.to_csv(out / "Fawkner_spatial_units.csv", index=False)
    print(units[["pseudobulk", "region", "pcw", "n_tissue_spots", "n_cells",
                 "epithelial_umi_fraction", "stromal_immune_umi_fraction", "eligible"]].to_string(index=False))


if __name__ == "__main__":
    main()
