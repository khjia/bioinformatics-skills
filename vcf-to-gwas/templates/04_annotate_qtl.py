#!/usr/bin/env python3
"""04_annotate_qtl.py ---- annotate Fast3VmrMLM QTN with the nearest gene.

Env:
  PROJECT   analysis root (required)
  GFF       gene GFF3 (default <PROJECT>/final.gene.mRNA.longest.gff3)
  WINDOW    max distance to gene, bp (default 50000)
  OUT       output TSV (default <PROJECT>/06.Fast3VmrMLM/report/qtl_annotated.tsv)
  TRAITS    comma-separated subset (default: all traits in trait.csv)

Output columns: trait, marker, chr, pos, lod, significance, gene_id, gene_name, dist
QTN farther than WINDOW from any gene are kept with gene_id = "intergenic".
"""
import os
import re
import sys
import numpy as np
import pandas as pd

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("set PROJECT")
GFF    = os.environ.get("GFF", os.path.join(PROJECT, "final.gene.mRNA.longest.gff3"))
WINDOW = int(os.environ.get("WINDOW", "50000"))
OUT    = os.environ.get("OUT", os.path.join(PROJECT, "06.Fast3VmrMLM/report/qtl_annotated.tsv"))
RESDIR = os.path.join(PROJECT, "06.Fast3VmrMLM/results")
TRAIT_CSV = os.path.join(PROJECT, "06.Fast3VmrMLM/inputs/trait.csv")

traits = list(pd.read_csv(TRAIT_CSV, nrows=0).columns)[1:]
sel = os.environ.get("TRAITS", "").strip()
if sel:
    traits = [t.strip() for t in sel.split(",") if t.strip()]

# ---- genes ----
gff = pd.read_csv(GFF, sep="\t", header=None, comment="#",
                  names=["chr", "src", "type", "start", "end", "score", "strand", "phase", "attrs"])
gff = gff[gff["type"] == "gene"].copy()

def attr(s, key):
    m = re.search(key + r'=([^;]+)', str(s))
    return m.group(1) if m else ""

gff["gid"] = gff["attrs"].apply(lambda s: attr(s, "ID"))
gff["gname"] = gff["attrs"].apply(lambda s: attr(s, "Name"))
genes = gff[["chr", "start", "end", "gid", "gname"]].reset_index(drop=True)
print(f"genes: {len(genes)}")

def pick(cols, *keys):
    for c in cols:
        low = str(c).lower()
        if any(k.lower() in low for k in keys):
            return c
    return None

rows = []
for tr in traits:
    xf = os.path.join(RESDIR, f"{tr}_result.xlsx")
    if not os.path.exists(xf):
        print(f"  skip (no xlsx): {tr}")
        continue
    xl = pd.read_excel(xf)
    if len(xl) == 0:
        print(f"  skip (empty): {tr}")
        continue
    cm, cc = pick(xl.columns, "marker"), pick(xl.columns, "chromosome", "chr")
    cp, cl = pick(xl.columns, "position", "pos"), pick(xl.columns, "lod")
    cs = pick(xl.columns, "significance")
    if not all([cm, cc, cp, cl]):
        print(f"  skip (columns not found): {tr} -> {list(xl.columns)}")
        continue

    for _, r in xl.iterrows():
        ch = str(r[cc]); pos = int(float(r[cp])); lod = float(r[cl])
        sig = str(r[cs]) if cs else "SIG"
        cg = genes[genes["chr"] == ch]
        if len(cg) == 0:
            rows.append(dict(trait=tr, marker=r[cm], chr=ch, pos=pos, lod=lod,
                             significance=sig, gene_id="intergenic", gene_name="", dist=-1))
            continue
        s = cg["start"].to_numpy()
        e = cg["end"].to_numpy()
        d = np.minimum(np.abs(s - pos), np.abs(e - pos))
        d = np.where((s <= pos) & (e >= pos), 0, d)
        j = int(np.argmin(d))
        if d[j] <= WINDOW:
            n = cg.iloc[j]
            rows.append(dict(trait=tr, marker=r[cm], chr=ch, pos=pos, lod=lod,
                             significance=sig, gene_id=n["gid"], gene_name=n["gname"],
                             dist=int(d[j])))
        else:
            rows.append(dict(trait=tr, marker=r[cm], chr=ch, pos=pos, lod=lod,
                             significance=sig, gene_id="intergenic", gene_name="",
                             dist=int(d[j])))
    print(f"  {tr}: {len(xl)} QTN")

adf = pd.DataFrame(rows)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
adf.to_csv(OUT, sep="\t", index=False)
genic = (adf["gene_id"] != "intergenic").sum()
print(f"wrote {OUT}: {len(adf)} QTN ({genic} genic)")
