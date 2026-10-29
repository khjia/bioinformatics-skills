#!/usr/bin/env python3
"""26_known_loci.py ---- look up known associations in the NHGRI-EBI GWAS Catalog (v2 API).

Human-centric: use rsIDs or human gene symbols. For non-model plants, query the
orthologous Arabidopsis/rice gene symbol and treat hits as hypothesis-generating.

Usage:
  python3 26_known_loci.py --rsid rs7903146 --out known.tsv
  python3 26_known_loci.py --gene ELOVL6 --out known.tsv --max 50
"""
import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
UA = {"User-Agent": "gwas-known-loci/1.1"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 410:
            sys.exit("GWAS Catalog v1 API is retired; this script uses v2 "
                     "(https://www.ebi.ac.uk/gwas/rest/api/v2).")
        raise


def rows_from(payload):
    out = []
    for a in payload.get("_embedded", {}).get("associations", []):
        rsids = [x.get("rs_id", "") for x in a.get("snp_allele", []) or []]
        traits = [t.get("efo_trait", "") for t in a.get("efo_traits", []) or []]
        if not traits:
            traits = a.get("reported_trait", []) or []
        out.append({
            "rsid": ";".join([x for x in rsids if x]),
            "traits": ";".join([t for t in traits if t]),
            "genes": ";".join(a.get("mapped_genes", []) or []),
            "pvalue": a.get("p_value", ""),
            "beta": a.get("beta", ""),
            "risk_allele": ";".join(a.get("snp_effect_allele", []) or []),
            "study": a.get("accession_id", ""),
            "pmid": a.get("pubmed_id", ""),
            "author": a.get("first_author", ""),
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rsid")
    ap.add_argument("--gene")
    ap.add_argument("--max", type=int, default=200)
    ap.add_argument("--out", default="known_loci.tsv")
    a = ap.parse_args()
    if not (a.rsid or a.gene):
        sys.exit("provide --rsid or --gene")

    if a.rsid:
        url = f"{BASE}/associations?rs_id={urllib.parse.quote(a.rsid)}&size={a.max}"
    else:
        url = f"{BASE}/associations?mapped_gene={urllib.parse.quote(a.gene)}&size={a.max}"
    rows = rows_from(get(url))

    cols = ["rsid", "traits", "genes", "pvalue", "beta", "risk_allele", "study", "pmid", "author"]
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows[: a.max]:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"{len(rows[:a.max])} associations -> {a.out}")


if __name__ == "__main__":
    main()
