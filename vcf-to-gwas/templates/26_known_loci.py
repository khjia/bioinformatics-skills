#!/usr/bin/env python3
"""26_known_loci.py ---- plant known-locus check via Europe PMC literature search.

No human GWAS Catalog: the query is a gene (typically the Arabidopsis/rice
ortholog of a candidate) plus trait keywords, and the output is the most-cited
matching papers. This is literature evidence, not a curated association database.

Usage:
  python3 26_known_loci.py --gene AT5G43760 --trait "fatty acid seed oil" \
      --out known.tsv --meta known_meta.tsv --max 20
"""
import argparse
import json
import sys
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "plant-known-loci/1.0"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True, help="gene symbol or locus ID, e.g. AT5G43760")
    ap.add_argument("--trait", default="", help="trait keywords, e.g. 'fatty acid seed oil'")
    ap.add_argument("--max", type=int, default=20)
    ap.add_argument("--out", default="known_loci.tsv")
    ap.add_argument("--meta", default="", help="optional metadata TSV")
    a = ap.parse_args()

    query = a.gene if not a.trait else f"{a.gene} AND ({a.trait})"
    url = (f"{BASE}?query={urllib.parse.quote(query)}&format=json"
           f"&pageSize={a.max}&resultType=core&sort=CITED%20desc")
    payload = get(url)
    hits = payload.get("resultList", {}).get("result", [])

    cols = ["pmid", "doi", "year", "title", "journal", "first_author", "cited_by"]
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for h in hits:
            j = (h.get("journalInfo") or {}).get("journal") or {}
            row = {
                "pmid": h.get("pmid", ""),
                "doi": h.get("doi", ""),
                "year": h.get("pubYear", ""),
                "title": (h.get("title") or "").replace("\t", " "),
                "journal": j.get("title", ""),
                "first_author": h.get("authorString", "").split(",")[0],
                "cited_by": h.get("citedByCount", ""),
            }
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    if a.meta:
        with open(a.meta, "w", encoding="utf-8") as fh:
            fh.write("gene\ttrait\tn_hits\tsource\n")
            fh.write(f"{a.gene}\t{a.trait}\t{payload.get('hitCount', len(hits))}\tEuropePMC\n")

    print(f"{len(hits)} papers (hitCount {payload.get('hitCount')}) -> {a.out}")


if __name__ == "__main__":
    main()
