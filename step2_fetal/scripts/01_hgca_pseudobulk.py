#!/usr/bin/env python3
"""Inspect HGCA metadata and create donor-level raw-count pseudobulks."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROLIFERATIVE_LABELS = {
    "Stem cells",
    "Proximal progenitor",
    "Distal progenitor",
    "TA",
    "Progenitor (NEUROG3+)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-cells", type=int, default=50)
    return parser.parse_args()


def sum_rows(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).ravel()
    return np.asarray(matrix).sum(axis=0)


def build_pseudobulk(adata, metadata: pd.DataFrame, mask: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = metadata.loc[mask].copy()
    rows = []
    manifest = []
    for donor, donor_meta in selected.groupby("Sample name", observed=True, sort=True):
        positions = metadata.index.get_indexer(donor_meta.index)
        values = sum_rows(adata.X[positions, :])
        stage = "fetal" if donor_meta["Diagnosis"].iloc[0] == "fetal" else "adult"
        rows.append(pd.Series(values, index=adata.var_names, name=donor))
        manifest.append(
            {
                "donor_id": donor,
                "stage": stage,
                "age": donor_meta["Age"].iloc[0],
                "n_cells": len(donor_meta),
                "regions": ";".join(sorted(donor_meta["Region code"].astype(str).unique())),
                "annotations": ";".join(sorted(donor_meta["annotation"].astype(str).unique())),
            }
        )
    matrix = pd.concat(rows, axis=1)
    if not np.allclose(matrix.to_numpy(), np.rint(matrix.to_numpy())):
        raise ValueError("HGCA raw-count matrix contains non-integer pseudobulk sums")
    return matrix.astype(np.int64), pd.DataFrame(manifest)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()

    base = (
        obs["Diagnosis"].isin(["fetal", "Healthy adult"])
        & (obs["category"] == "Epithelial")
        & (obs["Region"] != "lymph node")
        & ~obs["predicted_doublets"]
    )
    counts = obs.loc[base].groupby(["Diagnosis", "Sample name"], observed=True).size()
    eligible = set(counts[counts >= args.min_cells].index.get_level_values("Sample name"))
    primary_mask = base & obs["Sample name"].isin(eligible)

    primary, primary_manifest = build_pseudobulk(adata, obs, primary_mask)
    primary.to_csv(args.output_dir / "HGCA_primary_pseudobulk_counts.tsv.gz", sep="\t")
    primary_manifest.to_csv(args.output_dir / "HGCA_primary_sample_manifest.csv", index=False)

    proliferative_mask = (
        primary_mask
        & (obs["Diagnosis"] == "Healthy adult")
        & obs["annotation"].isin(PROLIFERATIVE_LABELS)
    )
    proliferative_counts = obs.loc[proliferative_mask].groupby("Sample name", observed=True).size()
    eligible_proliferative = set(proliferative_counts[proliferative_counts >= args.min_cells].index)
    control_mask = (primary_mask & (obs["Diagnosis"] == "fetal")) | (
        proliferative_mask & obs["Sample name"].isin(eligible_proliferative)
    )
    control, control_manifest = build_pseudobulk(adata, obs, control_mask)
    control.to_csv(args.output_dir / "HGCA_proliferative_pseudobulk_counts.tsv.gz", sep="\t")
    control_manifest.to_csv(args.output_dir / "HGCA_proliferative_sample_manifest.csv", index=False)

    obs.reset_index().to_csv(args.output_dir / "HGCA_cell_metadata.csv.gz", index=False)
    gene_manifest = adata.var.reset_index().rename(columns={"index": "gene"})
    gene_manifest.to_csv(args.output_dir / "HGCA_gene_manifest.csv", index=False)
    print(
        f"Primary: {sum(primary_manifest.stage == 'fetal')} fetal and "
        f"{sum(primary_manifest.stage == 'adult')} adult donors. "
        f"Control: {sum(control_manifest.stage == 'adult')} adult proliferative donors."
    )


if __name__ == "__main__":
    main()
