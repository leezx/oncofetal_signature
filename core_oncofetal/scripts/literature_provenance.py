#!/usr/bin/env python3
"""Axis L (Core method v3.0): literature provenance from the full-text audit.

Input: config/literature_audit_31.tsv, one record per gene x study, classified
from the study's own data after reading the full text (audit corpus: the 55
Step 1 papers plus Mustata 2013, Pikkupeura 2023, Elmentaite 2021,
Fernandez-Vallone 2016, Karo-Atar 2022, Vaquero-Siguero 2026).
  A  fetal intestine: in vivo, or fetal-derived spheroid/organoid vs adult
  B  injury/infection/YAP-driven regenerative fetal-like or revival state
  C  CRC oncofetal / fetal-like tumour state
  mention / not_found / not_read / other_tissue: recorded, never counted
Rule (unchanged from v2.0, now applied to audited evidence):
  L-A: >= 1 audited primary study with A evidence
  L-B: >= 2 distinct audited primary studies with B evidence
  pass = L-A or L-B
Usage (repo root): python3 core_oncofetal/scripts/literature_provenance.py
"""
import pathlib
import re

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
AUDIT = ROOT / "core_oncofetal/config/literature_audit_31.tsv"
CURATED = ROOT / "step1_Literature_curated_intestinal_oncofetal_marker_candidates_merged_v1_v2.csv"
OUT = ROOT / "core_oncofetal/config/literature_provenance_31.tsv"


def main():
    cand = pd.read_csv(ROOT / "step2_fetal/config/literature_candidates_31.tsv", sep="\t")
    au = pd.read_csv(AUDIT, sep="\t").fillna("")
    cur = pd.read_csv(CURATED, encoding="utf-8-sig")
    cur["gene"] = cand.set_index("literature_label").gene.reindex(cur.Oncofetal_marker).to_numpy()
    cur = cur.set_index("gene")
    rows = []
    for g, label in zip(cand.gene, cand.literature_label):
        a = au[au.gene == g]
        studies = lambda c: sorted(set(a[a.evidence_class == c].study))
        A, B, C = studies("A"), studies("B"), studies("C")
        other = a[~a.evidence_class.isin(["A", "B", "C"])]
        code = re.findall(r"\(([^()]*)\)\.?\s*$", cur.loc[g, "Annotation"].strip())[0]
        la, lb = len(A) >= 1, len(B) >= 2
        rows.append(dict(gene=g, literature_label=label, step1_curator_code=code,
                         audited_classes="".join(c for c, s in (("A", A), ("B", B), ("C", C)) if s),
                         n_A_studies=len(A), A_studies="; ".join(A), n_B_studies=len(B), B_studies="; ".join(B),
                         n_C_studies=len(C), C_studies="; ".join(C),
                         unsupported_or_unread="; ".join(f"{r.study} [{r.evidence_class}]" for r in other.itertuples()),
                         step1_references=cur.loc[g, "References"],
                         L_A_direct_fetal=la, L_B_replicated_revival=lb,
                         L_call="pass" if (la or lb) else "fail", L_route="L-A" if la else ("L-B" if lb else "none")))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(out[["gene", "step1_curator_code", "audited_classes", "n_A_studies", "n_B_studies", "n_C_studies",
               "L_call", "L_route"]].to_string(index=False))


if __name__ == "__main__":
    main()
