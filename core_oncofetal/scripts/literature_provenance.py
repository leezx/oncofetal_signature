#!/usr/bin/env python3
"""Axis L (Core method v2.0): deterministic literature provenance for the 31 candidates.

Curator code = final parenthetical class of the annotation (A fetal intestine,
B regenerative/fetal-like, C CRC). Named studies = citations in References
matching "Author et al., Journal, Year" or "Author & Author, Journal, Year";
generic phrases are not studies.
  L-A: code contains A (hedged forms included)
  L-B: code contains B and >= 2 distinct named studies
  pass = L-A or L-B
Usage (repo root): python3 core_oncofetal/scripts/literature_provenance.py
"""
import pathlib
import re

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "step1_Literature_curated_intestinal_oncofetal_marker_candidates_merged_v1_v2.csv"
OUT = ROOT / "core_oncofetal/config/literature_provenance_31.tsv"
NAMED = re.compile(r"^[^,]+(?: et al\.| & [^,]+)?, [^,]+, (19|20)\d{2}$")


def main():
    d = pd.read_csv(SRC, encoding="utf-8-sig")
    cand = pd.read_csv(ROOT / "step2_fetal/config/literature_candidates_31.tsv", sep="\t")
    d["gene"] = cand.set_index("literature_label").gene.reindex(d.Oncofetal_marker).to_numpy()
    if d.gene.isna().any():
        raise SystemExit("Unmatched literature labels: " + str(d[d.gene.isna()].Oncofetal_marker.tolist()))
    rows = []
    for _, r in d.iterrows():
        code = re.findall(r"\(([^()]*)\)\.?\s*$", r.Annotation.strip())[0]
        classes = "".join(c for c in "ABC" if c in re.sub(r"[^ABC]", "", code.split(",")[0]))
        refs = [x.strip() for x in str(r.References).split(";")]
        named = [x for x in refs if NAMED.match(x)]
        generic = [x for x in refs if x not in named]
        la = "A" in classes
        lb = "B" in classes and len(set(named)) >= 2
        rows.append(dict(gene=r.gene, literature_label=r.Oncofetal_marker, curator_code=code, classes=classes,
                         n_named_studies=len(set(named)), named_studies="; ".join(named),
                         generic_provenance="; ".join(generic), L_A_direct_fetal=la, L_B_replicated_revival=lb,
                         L_call="pass" if (la or lb) else "fail",
                         L_route="L-A" if la else ("L-B" if lb else "none")))
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(out[["gene", "curator_code", "classes", "n_named_studies", "generic_provenance", "L_call", "L_route"]].to_string(index=False))


if __name__ == "__main__":
    main()
