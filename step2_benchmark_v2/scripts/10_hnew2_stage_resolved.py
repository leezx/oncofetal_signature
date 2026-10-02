#!/usr/bin/env python3
"""Addendum v1.3: stage-resolved H-new2 (Fawkner fetal vs Burclaff adult).

Fetal units are the H-new2 units unchanged; PCW comes from the authors'
sample key (step2_fetal_state/config/sample_key_fawkner.tsv). Breakpoint 9 PCW
(frozen in v1.2). Writes unit tables and runs the generic edgeR module:
  Hnew2_early (early vs adult), Hnew2_late (mid/late vs adult),
  Fawkner_late_vs_early (within-fetal). Minimum group size 2.
"""
import argparse
import pathlib
import subprocess

import pandas as pd

BREAK = 9
QC_FAILED = {"ABF1": "25.6% sex-discordant cells", "AAU2": "pool 4 hashtagging failed"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--tables-dir", required=True)
    a = ap.parse_args()
    repo, work, de, tab = map(pathlib.Path, (a.repo, a.work_dir, a.de_dir, a.tables_dir))
    u = pd.read_csv(work / "Hnew2_units.csv")
    key = pd.read_csv(repo / "step2_fetal_state/config/sample_key_fawkner.tsv", sep="\t")
    fetal = u.group == "Fetal"
    tag = u.hashtag.astype(str).str.split("-").str[0]
    tag = tag.where(u.pool != "EPI_run2", tag.str[0])
    u["hashtag_totalseq"] = tag.where(fetal)
    u = u.merge(key[["pool", "hashtag_totalseq", "sample", "donor", "pcw", "region"]],
                on=["pool", "hashtag_totalseq"], how="left")
    if u.loc[fetal.to_numpy(), "pcw"].isna().any():
        raise SystemExit("Fetal unit without PCW in the sample key")
    u["stage_group"] = "Adult"
    u.loc[u.group == "Fetal", "stage_group"] = (u.pcw < BREAK).map({True: "Fetal_early", False: "Fetal_midlate"})
    u["later_hashtag_qc"] = u["sample"].map(QC_FAILED).fillna("")
    u.to_csv(tab / "Hnew2_stage_resolved_units.csv", index=False)
    print(u[u.eligible].groupby("stage_group").size().to_string())
    counts = work / "Hnew2_Fawkner_fetal_Burclaff_adult_counts.csv.gz"
    edger = repo / "step3_cancer/scripts/04_pseudobulk_edger.R"
    for name, case, ref in (("Hnew2_early_Fawkner_vs_Burclaff", "Fetal_early", "Adult"),
                            ("Hnew2_late_Fawkner_vs_Burclaff", "Fetal_midlate", "Adult"),
                            ("Fawkner_late_vs_early", "Fetal_midlate", "Fetal_early")):
        s = u.assign(group=u.stage_group)
        f = work / f"{name}_units.csv"
        s.to_csv(f, index=False)
        subprocess.run(["Rscript", str(edger), str(counts), str(f), case, ref, "none", str(de / f"{name}.csv"), "2"],
                       check=True)
        r = pd.read_csv(de / f"{name}.csv").set_index("symbol").loc["TNFRSF12A"]
        print(f"{name}: TNFRSF12A {r.log2FC:+.2f} P {r.PValue:.3g}")


if __name__ == "__main__":
    main()
