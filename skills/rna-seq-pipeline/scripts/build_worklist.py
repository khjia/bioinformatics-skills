#!/usr/bin/env python3
"""
Query ENA filereport API for a list of BioProjects or accessions and build
a parallel-download worklist for download_parallel.sh.

Usage:
    python3 build_worklist.py \
        --projects PRJNA267198 PRJNA887921 \
        --dst /path/to/01.rawdata \
        --out worklist.txt \
        [--layout PAIRED] \
        [--platform ILLUMINA BGISEQ]

Output: worklist.txt with lines "URL|DST|MD5"
Also writes runinfo.tsv alongside (full metadata for audit trail).
"""
import argparse, sys, urllib.request, urllib.parse, os, ssl

PROXY = os.environ.get("HTTPS_PROXY", "http://127.0.0.1:10809")
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,instrument_platform,scientific_name"


def fetch(accession: str) -> str:
    params = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": FIELDS,
        "format": "tsv",
    })
    url = f"{ENA}?{params}"
    proxy_handler = urllib.request.ProxyHandler({"http": PROXY, "https": PROXY})
    opener = urllib.request.build_opener(proxy_handler)
    with opener.open(url, timeout=120) as r:
        return r.read().decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--projects", nargs="+", required=True,
                    help="Project / sample / run accessions (PRJ*, SRR*, ERR*, ...)")
    ap.add_argument("--dst", required=True, help="Destination directory for fastq.gz")
    ap.add_argument("--out", default="worklist.txt")
    ap.add_argument("--runinfo", default="runinfo.tsv")
    ap.add_argument("--layout", default=None, help="Filter library_layout (PAIRED|SINGLE)")
    ap.add_argument("--platform", nargs="+", default=None,
                    help="Filter instrument_platform (e.g. ILLUMINA BGISEQ)")
    args = ap.parse_args()

    header_written = False
    n_runs = 0
    n_files = 0
    with open(args.runinfo, "w") as rinfo, open(args.out, "w") as wl:
        for proj in args.projects:
            try:
                tsv = fetch(proj)
            except Exception as e:
                print(f"FETCH FAIL {proj}: {e}", file=sys.stderr)
                continue
            lines = tsv.strip().split("\n")
            if not lines or len(lines) < 2:
                print(f"EMPTY {proj}", file=sys.stderr)
                continue
            header = lines[0]
            if not header_written:
                rinfo.write(header + "\n")
                header_written = True
            cols = header.split("\t")
            idx = {c: i for i, c in enumerate(cols)}
            for ln in lines[1:]:
                f = ln.split("\t")
                if len(f) < len(cols):
                    f += [""] * (len(cols) - len(f))
                if args.layout and f[idx["library_layout"]] != args.layout:
                    continue
                if args.platform and f[idx["instrument_platform"]] not in args.platform:
                    continue
                urls = f[idx["fastq_ftp"]].split(";")
                md5s = f[idx["fastq_md5"]].split(";")
                if not urls or not urls[0]:
                    continue
                n_runs += 1
                rinfo.write(ln + "\n")
                for u, m in zip(urls, md5s):
                    wl.write(f"https://{u}|{args.dst}|{m}\n")
                    n_files += 1

    print(f"runs={n_runs} files={n_files} → {args.out}")
    print(f"metadata        → {args.runinfo}")


if __name__ == "__main__":
    main()
