#!/usr/bin/env python3
"""26_known_loci.py ---- look up known associations in the NHGRI-EBI GWAS Catalog.

Human-centric: use rsIDs or human gene symbols. For non-model plants, query the
orthologous Arabidopsis/rice gene symbol and treat hits as hypothesis-generating.

Usage:
  python3 26_known_loci.py --rsid rs7903146 --out known.tsv
  python3 26_known_loci.py --gene APOE --out known.tsv
  python3 26_known_loci.py --gene KCS19 --out known.tsv --max 50
"""
import argparse
import json
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/gwas/rest/api"
UA = {"User-Agent": "gwas-known-loci/1.0"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def assoc_rows(a):
    rows = []
    for ass in a.get("associations", []):
        loci = ass.get("loci", [])
        genes = []
        for l in loci:
            for g in l.get("authorReportedGenes", []) or []:
                genes.append(g.get("geneName", ""))
        traits = [t.get("trait", "") for t in ass.get("efoTraits", []) or []]
        rows.append({
            "rsid": (ass.get("snps", [{}])[0].get("rsId", "") if ass.get("snps") else ""),
            "traits": ";".join([t for t in traits if t]),
            "genes": ";".join([g for g in genes if g]),
            "pvalue": ass.get("pvalue", ""),
            "or_beta": ass.get("orPerCopyNum") or ass.get("betaNum") or "",
            "risk_allele": ass.get("riskAlleleName", ""),
            "study": ass.get("study", {}).get("accessionId", ""),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rsid")
    ap.add_argument("--gene")
    ap.add_argument("--max", type=int, default=200)
    ap.add_argument("--out", default="known_loci.tsv")
    a = ap.parse_args()
    if not (a.rsid or a.gene):
        sys.exit("provide --rsid or --gene")

    rows = []
    if a.rsid:
        url = f"{BASE}/singleNucleotidePolymorphisms/{urllib.parse.quote(a.rsid)}/associations"
        rows = assoc_rows(get(url))
    else:
        url = (f"{BASE}/singleNucleotidePolymorphisms/search/findByGene"
               f"?geneName={urllib.parse.quote(a.gene)}&size={a.max}")
        d = get(url)
        for snp in d.get("_embedded", {}).get("singleNucleotidePolymorphisms", []):
            time.sleep(0.2)
            try:
                rows += assoc_rows(get(f"{BASE}/singleNucleotidePolymorphisms/{snp['rsId']}/associations"))
            except Exception as e:
                print("skip", snp.get("rsId"), e, file=sys.stderr)

    cols = ["rsid", "traits", "genes", "pvalue", "or_beta", "risk_allele", "study"]
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows[: a.max]:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"{len(rows[:a.max])} associations -> {a.out}")


if __name__ == "__main__":
    main()
