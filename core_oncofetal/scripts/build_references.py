#!/usr/bin/env python3
"""Methods reference list generated from provenance files (no hand-typed citations).

Sources:
  literature  core_oncofetal/config/literature_audit_31.tsv: every study with
              A/B/C primary evidence for >= 1 of the 31 candidates; DOI from
              docs/step1_literature_candidates/paper_inventory.tsv, else
              core_oncofetal/config/reference_doi_supplement.tsv
  datasets    docs/DATASET_REGISTRY.csv rows used by the CIOC / ECOS-38 /
              Extended CIOC Methods (PMID and DOI from the registry)
  methods     core_oncofetal/config/reference_methods.tsv
Bibliographic metadata is fetched from Crossref by DOI and cached
(core_oncofetal/config/crossref_cache.json); a failed lookup is a hard error.
Outputs: core_oncofetal/results/Methods_references.{md,tsv}
Usage (repo root): python3 core_oncofetal/scripts/build_references.py
"""
import json
import pathlib
import re
import unicodedata
import urllib.request

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
CFG = ROOT / "core_oncofetal/config"
OUT = ROOT / "core_oncofetal/results"
CACHE = CFG / "crossref_cache.json"
DATASETS = [  # (registry Dataset, role in Methods)
    ("Human Gut Cell Atlas (HGCA)", "Gate H: HGCA ≥ 9 PCW vs adult"),
    ("Fawkner-Corbett fetal atlas", "Gate H: H-new2 fetal arm (≥ 9 PCW)"),
    ("Burclaff adult epithelium", "Gate H: H-new2 adult arm"),
    ("Gao 2018 fetal digestive tract", "Gate H: Gao fetal LI (GSE95630) vs GSE103154 adult LI"),
    ("Pikkupeura in vivo (sub-series IV)", "Gate M: GSE230581"),
    ("Hemmerling LCM ileum", "Supportive mouse replication: GSE44433"),
    ("Joanito CRC epithelial scRNA (5 cohorts: CRC-SG1, CRC-SG2, KUL3, KUL5, SMC)", "Gate C: Joanito"),
    ("Pelka CRC immune-hub atlas", "Gate C: Pelka (GSE178341)"),
    ("TCGA COAD + READ", "Supportive bulk tumour vs normal"),
    ("recount3", "TCGA/GTEx uniform processing"),
    ("Khaliq CRC single-cell atlas", "Gate E and coherence, primary (GSE200997)"),
    ("Che CRC + liver-metastasis atlas", "Gate E and coherence, replication (GSE178318)"),
    ("Tabula Sapiens large intestine", "Compartment annotation only"),
    ("Cao fetal atlas (sci-RNA-seq3)", "Breadth annotation: fetal tissue breadth (GSE156793)"),
    ("UCSC Toil TCGA + GTEx", "Breadth annotation: cancer reactivation breadth"),
]


def norm(s):
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]", "", s.split(" et al")[0]) + s.strip()[-4:]


def crossref(doi, cache):
    doi = doi.lower()
    if doi not in cache:
        req = urllib.request.Request(f"https://api.crossref.org/works/{doi}",
                                     headers={"User-Agent": "oncofetal-signature-references (mailto:zxli.hku@gmail.com)"})
        cache[doi] = json.load(urllib.request.urlopen(req, timeout=30))["message"]
    return cache[doi]


def fmt(m):
    au = m.get("author", [])
    names = [f"{a.get('family', '')} {''.join(p[0] for p in re.split(r'[ -]', a.get('given', '')) if p)}".strip()
             if "family" in a else a.get("name", "") for a in au]
    authors = ", ".join(names[:6]) + (", et al" if len(names) > 6 else "")
    if not names:
        authors = "Tabula Sapiens Consortium" if "Tabula" in m["title"][0] else ""
    year = next(((m.get(k) or {}).get("date-parts") or [[None]])[0][0] for k in ("published-print", "issued")
                if ((m.get(k) or {}).get("date-parts") or [[None]])[0][0])
    jr = (m.get("short-container-title") or m.get("container-title") or [""])[0].rstrip(".")
    vol, iss, pg = m.get("volume", ""), m.get("issue", ""), m.get("page") or m.get("article-number", "")
    loc = f"{vol}" + (f"({iss})" if iss else "") + (f":{pg}" if pg else "")
    title = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", m["title"][0])).strip().rstrip(".").replace(" :", ":")
    return f"{authors}. {title}. {jr}. {year}" + (f";{loc}" if loc else "") + f". doi:{m['DOI']}"


def main():
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    rows = []
    a = pd.read_csv(CFG / "literature_audit_31.tsv", sep="\t")
    a = a[a.evidence_class.isin(list("ABC"))]
    inv = pd.read_csv(ROOT / "docs/step1_literature_candidates/paper_inventory.tsv", sep="\t").dropna(subset=["doi"])
    inv_doi = {norm(c): d for c, d in zip(inv.citation, inv.doi) if c != "—"}
    sup = pd.read_csv(CFG / "reference_doi_supplement.tsv", sep="\t").set_index("study")
    for study, g in a.groupby("study"):
        doi = inv_doi.get(norm(study))
        src = "paper_inventory.tsv"
        if study in sup.index:
            doi, src = sup.loc[study, "doi"], "reference_doi_supplement.tsv: " + sup.loc[study, "resolution"]
        assert doi, f"no DOI for {study}"
        cited = "; ".join(f"{cl}: {', '.join(sorted(set(x.gene)))}" for cl, x in g.groupby("evidence_class"))
        rows.append(dict(section="Literature (Step 1 nomination; A/B/C primary evidence)", key=study, doi=doi,
                         cited_for=cited, doi_source=src))
    reg = pd.read_csv(ROOT / "docs/DATASET_REGISTRY.csv").set_index("Dataset")
    for ds, role in DATASETS:
        r = reg.loc[ds]
        rows.append(dict(section="Datasets", key=ds, doi=r.DOI, cited_for=f"{role}; accession {r.Accession}",
                         doi_source=f"DATASET_REGISTRY (PMID {int(r.PMID)})"))
    for _, r in pd.read_csv(CFG / "reference_methods.tsv", sep="\t").iterrows():
        rows.append(dict(section="Methods and software", key=r.key, doi=r.doi, cited_for=r.cited_for,
                         doi_source="reference_methods.tsv"))
    d = pd.DataFrame(rows)
    d["reference"] = [fmt(crossref(x, cache)) for x in d.doi]
    d["crossref_title"] = [crossref(x, cache)["title"][0] for x in d.doi]
    CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    # one entry per DOI; a paper cited both as literature and as a dataset keeps both roles
    extra = d[d.duplicated("doi", keep="first")].groupby("doi").cited_for.apply(" | ".join)
    d = d.drop_duplicates("doi", keep="first").copy()
    d["cited_for"] = [c + (" | " + extra[x] if x in extra.index else "") for c, x in zip(d.cited_for, d.doi)]
    d.to_csv(OUT / "Methods_references.tsv", sep="\t", index=False)
    md = ["# Methods references", "",
          "Generated by `core_oncofetal/scripts/build_references.py` from `literature_audit_31.tsv`, "
          "`paper_inventory.tsv`, `DATASET_REGISTRY.csv` and `reference_methods.tsv`; metadata from Crossref by DOI. "
          "Do not edit by hand.", ""]
    n = 0
    for sec, x in d.groupby("section", sort=False):
        md += [f"## {sec}", ""]
        for _, r in x.iterrows():
            n += 1
            md.append(f"{n}. {r.reference}  ")
            md.append(f"   *Cited for:* {r.cited_for}")
        md.append("")
    (OUT / "Methods_references.md").write_text("\n".join(md))
    print(d.groupby("section", sort=False).size().to_string())
    print(d[["key", "crossref_title"]].to_string(index=False))


if __name__ == "__main__":
    main()
