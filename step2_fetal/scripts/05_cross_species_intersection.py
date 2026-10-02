#!/usr/bin/env python3
"""Join H1, H2, and M1 through frozen Ensembl one-to-one orthologues."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--orthologues", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hgca = pd.read_csv(args.tables_dir / "HGCA_DE.csv")
    gao = pd.read_csv(args.tables_dir / "Gao_validation.csv").rename(
        columns={"gene": "human_gene"}
    )
    mouse = pd.read_csv(args.tables_dir / "Mouse_DE.csv")
    mouse = mouse.loc[mouse["mouse_gene"].notna()].copy()
    orth = pd.read_csv(args.orthologues, sep="\t")
    orth = orth.loc[
        (orth["Mouse homology type"] == "ortholog_one2one")
        & orth["Gene name"].notna()
        & orth["Mouse gene name"].notna()
    ].copy()
    orth = orth.rename(
        columns={
            "Gene stable ID": "human_ensembl_gene_id",
            "Gene name": "human_gene",
            "Mouse gene stable ID": "mouse_ensembl_gene_id",
            "Mouse gene name": "mouse_gene",
            "Mouse homology type": "orthology_type",
            "Mouse orthology confidence [0 low, 1 high]": "orthology_confidence",
        }
    )
    orth = orth.drop_duplicates(subset=["human_ensembl_gene_id", "mouse_ensembl_gene_id"])

    evidence = orth.merge(hgca, on="human_ensembl_gene_id", how="inner", validate="one_to_one")
    evidence = evidence.rename(columns={"human_gene": "human_gene_ensembl116", "gene": "human_gene"})
    evidence["human_symbol_matches_ensembl116"] = (
        evidence["human_gene"] == evidence["human_gene_ensembl116"]
    )
    evidence = evidence.merge(gao, on="human_gene", how="left", validate="one_to_one")
    evidence = evidence.merge(mouse, on="mouse_gene", how="left", validate="many_to_one")

    for column in [
        "HGCA_primary_pass",
        "HGCA_proliferative_direction_pass",
        "H1_pass",
        "Gao_pass",
        "Mouse_pass",
    ]:
        evidence[column] = evidence[column].fillna(False).astype(bool)
    evidence["Conserved_Fetal_High"] = (
        evidence["H1_pass"] & evidence["Gao_pass"] & evidence["Mouse_pass"]
    )
    evidence["evidence_label"] = "other"
    evidence.loc[
        evidence["HGCA_primary_pass"]
        & ~evidence["HGCA_proliferative_direction_pass"],
        "evidence_label",
    ] = "Proliferation-associated reject"
    evidence.loc[
        evidence["H1_pass"] & evidence["Mouse_pass"] & ~evidence["Gao_pass"],
        "evidence_label",
    ] = "Conserved discovery holdout"
    evidence.loc[
        evidence["H1_pass"] & evidence["Gao_pass"] & ~evidence["Mouse_pass"],
        "evidence_label",
    ] = "Human fetal-high only"
    evidence.loc[evidence["Conserved_Fetal_High"], "evidence_label"] = (
        "Conserved Fetal-High Candidate"
    )

    evidence = evidence.sort_values(
        ["Conserved_Fetal_High", "HGCA_log2FC", "Gao_log2FC", "Mouse_log2FC"],
        ascending=[False, False, False, False],
    )
    evidence.to_csv(args.tables_dir / "Cross_species_evidence.csv", index=False)
    conserved = evidence.loc[evidence["Conserved_Fetal_High"]].copy()
    conserved.to_csv(args.tables_dir / "Conserved_Fetal_High.csv", index=False)
    print(f"Wrote {len(evidence):,} mapped genes; {len(conserved):,} strict candidates.")


if __name__ == "__main__":
    main()
