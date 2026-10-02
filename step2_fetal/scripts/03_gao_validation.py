#!/usr/bin/env python3
"""Gao 2018 effect-only validation for fetal-high large-intestinal genes."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_expression(path: Path, cells: list[str] | None = None) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    usecols = ["Gene"] + cells if cells is not None else header
    return pd.read_csv(path, sep="\t", usecols=usecols).set_index("Gene")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    annotation_path = args.data_dir / "41556_2018_105_MOESM4_ESM.xlsx"
    fetal_matrix_path = args.data_dir / "GSE95630_Digestion_TPM_new.txt.gz"
    adult_matrix_path = args.data_dir / "GSE103154_All_Merge_umi_tpm_gene.txt.gz"

    annotation = pd.read_excel(annotation_path)
    fetal_meta = annotation.loc[
        (annotation["Tissue"] == "LI") & (annotation["CellType"] == "Epithelial")
    ].copy()
    # Suffixes such as embryo1.2 denote an additional library from embryo1.
    parsed = fetal_meta["Sample"].str.extract(r"_(\d+W)_embryo(\d+)(?:\.\d+)?_")
    if parsed.isna().any().any():
        raise ValueError("Could not parse fetal age/embryo from all selected sample names")
    fetal_meta["age"] = parsed[0]
    fetal_meta["embryo"] = parsed[1]
    fetal_meta["donor_id"] = fetal_meta["age"] + "_embryo" + fetal_meta["embryo"]
    # GEO matrix headers use lower-case "w" while the corrected supplement uses "W".
    fetal_meta["matrix_sample"] = fetal_meta["Sample"].str.replace(
        r"_(\d+)W_", r"_\1w_", regex=True
    )

    fetal = read_expression(fetal_matrix_path, fetal_meta["matrix_sample"].tolist())
    missing = sorted(set(fetal_meta["matrix_sample"]) - set(fetal.columns))
    if missing:
        raise ValueError(f"Selected fetal cells missing from matrix: {missing[:5]}")
    fetal_donor = pd.DataFrame(
        {
            donor: fetal.loc[:, group["matrix_sample"]].mean(axis=1)
            for donor, group in fetal_meta.groupby("donor_id", sort=True)
        }
    )

    adult = read_expression(adult_matrix_path)
    adult_columns = [column for column in adult.columns if column.startswith(("P1_", "P2_"))]
    if len(adult_columns) != adult.shape[1]:
        raise ValueError("Adult matrix contains columns that cannot be assigned to P1 or P2")
    adult_donor = pd.DataFrame(
        {
            "Adult1": adult.loc[:, [c for c in adult.columns if c.startswith("P1_")]].mean(axis=1),
            "Adult2": adult.loc[:, [c for c in adult.columns if c.startswith("P2_")]].mean(axis=1),
        }
    )

    genes = fetal_donor.index.intersection(adult_donor.index)
    fetal_donor = fetal_donor.loc[genes]
    adult_donor = adult_donor.loc[genes]
    fetal_median = fetal_donor.median(axis=1)
    adult_median = adult_donor.median(axis=1)

    result = pd.DataFrame(
        {
            "gene": genes,
            "Gao_log2FC": np.log2((fetal_median.to_numpy() + 1) / (adult_median.to_numpy() + 1)),
            "Fetal_mean_expression": fetal_donor.mean(axis=1).to_numpy(),
            "Fetal_median_expression": fetal_median.to_numpy(),
            "Adult1_expression": adult_donor["Adult1"].to_numpy(),
            "Adult2_expression": adult_donor["Adult2"].to_numpy(),
            "n_fetal_donors": fetal_donor.shape[1],
            "n_adult_donors": adult_donor.shape[1],
        }
    )
    result["Gao_direction_consistent"] = (
        (result["Fetal_median_expression"] > result["Adult1_expression"])
        & (result["Fetal_median_expression"] > result["Adult2_expression"])
    )
    result["Gao_pass"] = (result["Gao_log2FC"] >= 0.5) & result["Gao_direction_consistent"]
    result = result.sort_values(["Gao_pass", "Gao_log2FC"], ascending=[False, False])
    result.to_csv(args.output_dir / "Gao_validation.csv", index=False)

    sample_manifest = fetal_meta[
        ["Sample", "matrix_sample", "donor_id", "age", "embryo", "Tissue", "CellType"]
    ]
    sample_manifest.to_csv(args.output_dir / "Gao_fetal_cell_manifest.csv", index=False)
    print(
        f"Wrote {len(result):,} genes; {result['Gao_pass'].sum():,} pass H2; "
        f"{fetal_donor.shape[1]} fetal and {adult_donor.shape[1]} adult donors."
    )


if __name__ == "__main__":
    main()
