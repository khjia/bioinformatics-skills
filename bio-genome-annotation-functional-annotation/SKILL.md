---
name: bio-genome-annotation-functional-annotation
description: "Complete protein functional annotation pipeline: eggNOG-mapper (GO/KEGG/COG), Diamond blastp (SwissProt/TrEMBL/NR/Arabidopsis), InterProScan (Pfam/CDD/SMART/PRINTS/PANTHER). Merge all results into unified annotation table with per-DB coverage. Includes CNS cluster pitfalls and verification scripts."
---

# Functional Annotation — Complete Pipeline

Three-strategy functional annotation for predicted proteins:
1. **eggNOG-mapper** → GO, KEGG, COG, functional description
2. **Diamond blastp** → best hit against SwissProt/TrEMBL/NR/Arabidopsis (id≥30%, E-value≤1e-5)
3. **InterProScan** → domain/motif annotations (Pfam, CDD, SMART, PRINTS, PANTHER, etc.)

Plus: merge script + coverage statistics + verification.

## When to use

- New genome: annotate predicted protein-coding genes
- Input: protein FASTA file (from BRAKER4/EviAnn/etc.)
- Output: merged annotation TSV with GO, KEGG, domain info per gene

## CNS Cluster Pitfalls (MUST READ)

### eggNOG-mapper

| Issue | Fix |
|-------|-----|
| `emapper.py: command not found` | Use full path: `/media/nfs1/hermes/miniforge3/envs/eggnog/bin/emapper.py` |
| diamond not found in eggnogmapper/bin | `ln -sf /media/nfs1/hermes/miniforge3/envs/eggnog/bin/diamond /media/nfs1/hermes/miniforge3/envs/eggnog/lib/python3.11/site-packages/eggnogmapper/bin/diamond` |
| HMMER server crash during pfam_realign denovo | Use `--pfam_realign none` instead (GO/KEGG from orthologs, not Pfam denovo) |
| annotations file empty (only header) | Check log for Python traceback — usually HMMER server issue |

### InterProScan

| Issue | Fix |
|-------|-----|
| FunFam crashes: "Cannot run program python3" | Exclude FunFam + Gene3D with `-app` flag |
| Multi-format output with `-o` | Run twice: once for TSV, once for GFF3 |
| PANTHER data missing | Extract `panther-data-14.1.tar.gz` (~15GB) into `data/panther/19.0/` |

### Diamond

| Issue | Fix |
|-------|-----|
| `--salltitles` not supported in makedb | Remove `--salltitles` from `diamond makedb` |
| `--temp-dir` not recognized | Remove `--temp-dir` from `diamond blastp` |

### Protein databases

| DB | Size | Path |
|----|------|------|
| SwissProt | ~90MB | `/media/nfs1/hermes/db/protein_db/uniprot_sprot.fasta.gz` |
| TrEMBL | ~38GB | `/media/nfs1/hermes/db/protein_db/uniprot_trembl.fasta.gz` |
| NR | ~187GB | `/media/nfs1/hermes/db/protein_db/nr.gz` |
| Arabidopsis | ~9MB | `/media/nfs1/hermes/db/protein_db/arabidopsis_thaliana_reviewed.fa` |
| eggNOG | ~44GB | `/media/nfs1/hermes/db/eggnog/` (symlinked from /genome/eggnog_data/) |

---

## Strategy 1: eggNOG-mapper

### Script

```bash
#!/bin/bash
# Usage: bash run_eggnog.sh <sample> <proteins.fa> <output_dir>
SAMPLE=${1:-Pv25_08}
PROTEIN=${2:-proteins.fa}
OUTDIR=${3:-./eggNOG}
EMAPPER=/media/nfs1/hermes/miniforge3/envs/eggnog/bin/emapper.py

export http_proxy=http://127.0.0.1:10809
export https_proxy=http://127.0.0.1:10809

$EMAPPER \
    -i "$PROTEIN" \
    -o "$OUTDIR/${SAMPLE}_eggnog" \
    --data_dir /media/nfs1/hermes/db/eggnog \
    --cpu 40 \
    --tax_scope 33090 \
    --go_evidence non-electronic \
    --pfam_realign none \
    --temp_dir /tmp \
    --override
```

### Output columns (0-indexed after gene name)

| idx | Column |
|-----|--------|
| 0 | seed_ortholog |
| 5 | COG_category |
| 6 | Description |
| 7 | Preferred_name |
| 8 | GOs |
| 10 | KEGG_ko |
| 11 | KEGG_Pathway |

---

## Strategy 2: Diamond blastp

### Script

```bash
#!/bin/bash
# Usage: bash run_diamond.sh <sample> <proteins.fa> <output_dir>
SAMPLE=${1:-Pv25_08}
PROTEIN=${2:-proteins.fa}
OUTDIR=${3:-./diamond}
DB_DIR=/media/nfs1/hermes/db/protein_db

export PATH=/media/nfs1/hermes/miniforge3/bin:$PATH
export http_proxy=http://127.0.0.1:10809
export https_proxy=http://127.0.0.1:10809

# Build databases (skip if already built)
for DB_NAME in swissprot trembl nr arabidopsis; do
    [ -f "$OUTDIR/${DB_NAME}.dmnd" ] && continue
    case $DB_NAME in
        swissprot) IN=$DB_DIR/uniprot_sprot.fasta.gz;;
        trembl)    IN=$DB_DIR/uniprot_trembl.fasta.gz;;
        nr)        IN=$DB_DIR/nr.gz;;
        arabidopsis) IN=$DB_DIR/arabidopsis_thaliana_reviewed.fa;;
    esac
    diamond makedb --in $IN --db $OUTDIR/$DB_NAME
done

# Search (id≥30%, E-value≤1e-5, top 5 hits)
for DB_NAME in swissprot trembl nr arabidopsis; do
    diamond blastp \
        -d $OUTDIR/$DB_NAME \
        -q "$PROTEIN" \
        -o $OUTDIR/${SAMPLE}_vs_${DB_NAME}.blastp \
        --id 30 --evalue 1e-5 \
        --max-target-seqs 5 --outfmt 6 --threads 40
done
```

---

## Strategy 3: InterProScan

### Script (exclude FunFam/Gene3D)

```bash
#!/bin/bash
# Usage: bash run_interproscan.sh <sample> <proteins.fa> <output_dir>
SAMPLE=${1:-Pv25_08}
PROTEIN=${2:-proteins.fa}
OUTDIR=${3:-./interproscan}
IPS=/media/nfs2/hermes/soft/interproscan/interproscan-5.78-109.0/interproscan.sh

# TSV output
$IPS -i "$PROTEIN" -f tsv \
    -app AntiFam-8.0,CDD-3.21,Coils-2.2.1,Hamap-2026_01,MobiDBLite-4.0,NCBIfam-19.0,PANTHER-19.0,Pfam-38.2,PIRSF-3.10,PIRSR-2025_05,PRINTS-42.0,ProSitePatterns-2026_01,ProSiteProfiles-2026_01,SFLD-4,SMART-9.0,SUPERFAMILY-1.75 \
    -goterms -pa -dp -cpu 40 \
    -o $OUTDIR/${SAMPLE}_interproscan.tsv -T /tmp

# GFF3 output (separate run)
$IPS -i "$PROTEIN" -f gff3 \
    -app AntiFam-8.0,CDD-3.21,Coils-2.2.1,Hamap-2026_01,MobiDBLite-4.0,NCBIfam-19.0,PANTHER-19.0,Pfam-38.2,PIRSF-3.10,PIRSR-2025_05,PRINTS-42.0,ProSitePatterns-2026_01,ProSiteProfiles-2026_01,SFLD-4,SMART-9.0,SUPERFAMILY-1.75 \
    -goterms -pa -dp -cpu 40 \
    -o $OUTDIR/${SAMPLE}_interproscan.gff3 -T /tmp
```

### InterProScan TSV columns

| Col | Name |
|-----|------|
| 0 | protein_id |
| 3 | analysis (Pfam, CDD, SMART, etc.) |
| 4 | signature_accession |
| 5 | signature_description |
| 8 | evalue |
| 11 | interpro_accession |
| 13 | go_terms |
| 14 | pathways |

---

## Merge Script

Combines all three strategies into one unified annotation table.

### Step 1: Preprocess InterProScan with awk (fast, handles 1.5GB TSV)

```bash
# Skip MobiDBLite and Coils (low info), keep top 50 domains per gene
awk -F'\t' 'NR>4 && $6!="" && $6!="-" && $4!="MobiDBLite" && $4!="Coils" {
    gene=$1; db=$4; desc=$6;
    if (!(gene in cnt) || cnt[gene]<50) {
        key=gene SUBSEP cnt[gene]
        domains[key]=db":"desc
        cnt[gene]++
    }
}
END {
    for (g in cnt) {
        r=""
        for (i=0; i<cnt[g]; i++) {
            if (i>0) r=r";"
            r=r domains[g SUBSEP i]
        }
        print g "\t" r
    }
}' interproscan.tsv > interpro_summary.tsv
```

### Step 2: Python merge (reads eggNOG + Diamond + preprocessed InterPro)

```python
#!/usr/bin/env python3
"""Three-strategy functional annotation merge."""
import sys
from collections import defaultdict

def parse_eggnog(f):
    r = {}
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        p = line.strip().split('\t')
        if len(p) < 21: continue
        r[p[0]] = p[1:]
    return r

def parse_diamond(files):
    r = defaultdict(dict)
    for f in files:
        db = f.split('_vs_')[1].replace('.blastp','')
        for line in open(f):
            if line.startswith('#'): continue
            p = line.strip().split('\t')
            if len(p) < 12: continue
            gene, hit, pid, eval_, bits = p[0], p[1], float(p[2]), float(p[10]), float(p[11])
            if db not in r[gene] or bits > r[gene][db][3]:
                r[gene][db] = (hit, pid, eval_, bits)
    return r

def parse_interpro_summary(f):
    r = {}
    for line in open(f):
        p = line.strip().split('\t')
        if len(p) >= 2: r[p[0]] = p[1]
    return r

def parse_interpro_per_db(filepath):
    """Best InterPro hit per database per gene."""
    result = defaultdict(dict)
    for line in open(filepath):
        if line.startswith('#') or not line.strip(): continue
        p = line.strip().split('\t')
        if len(p) < 10: continue
        gene, db, desc = p[0], p[3], p[5]
        evalue = float(p[8]) if p[8] not in ['-', ''] else 999
        if db in ['MobiDBLite', 'Coils'] or desc in ['', '-']: continue
        if db not in result[gene] or evalue < result[gene][db][1]:
            result[gene][db] = (desc, evalue)
    return result

# Usage:
# python3 merge_annotations.py <sample> <proteins.fa> \
#   <eggnog.annotations> <interpro_raw.tsv> <interpro_summary.tsv> \
#   <diamond_files...> <output.tsv>
```

### Output format (26 columns)

```
Gene_ID | Description | Preferred_name | COG | GOs | KEGG_ko | KEGG_Pathway |
Diamond_best_hit | Diamond_best_pident | Diamond_best_evalue | Diamond_best_DB |
InterPro_Domains (summary) |
InterPro_Pfam | InterPro_CDD | InterPro_SMART | InterPro_PRINTS | InterPro_PANTHER |
InterPro_Gene3D | InterPro_SUPERFAMILY | InterPro_FunFam | InterPro_ProSiteProfiles |
InterPro_ProSitePatterns | InterPro_NCBIfam | InterPro_PIRSF | InterPro_Hamap | InterPro_SFLD
```

Diamond best hit prefers SwissProt > NR > Arabidopsis > TrEMBL.

---

## Coverage Statistics

```bash
# After merge, verify coverage
F=output_functional_annotation.tsv
IN=$(grep -c '^>' input_proteins.fa)
OUT=$(tail -n+2 $F | cut -f1 | sort -u | wc -l)
echo "Input: $IN  Merged: $OUT  Missing: $((IN-OUT))"

# Per-strategy counts
echo "GO:     $(tail -n+2 $F | awk -F'\t' '$5!="-" && $5!=""' | wc -l)"
echo "KEGG:   $(tail -n+2 $F | awk -F'\t' '$7!="-" && $7!=""' | wc -l)"
echo "Diamond: $(tail -n+2 $F | awk -F'\t' '$8!="-" && $8!=""' | wc -l)"
echo "InterPro: $(tail -n+2 $F | awk -F'\t' '$12!="-" && $12!=""' | wc -l)"

# Per-DB InterPro coverage (vs raw InterProScan)
for DB in Pfam CDD SMART PRINTS PANTHER; do
    ORIG=$(awk -F'\t' -v db="$DB" 'NR>4 && $4==db {print $1}' interproscan.tsv | sort -u | wc -l)
    MERGED=$(awk -F'\t' -v db="$DB" 'NR>1 && $12~db {print $1}' $F | sort -u | wc -l)
    echo "$DB: $MERGED/$ORIG"
done

# Unannotated gene length distribution
python3 -c "
import subprocess
lengths = {}
name = None; seq = []
for line in open('proteins.fa'):
    line = line.strip()
    if line.startswith('>'):
        if name: lengths[name] = len(''.join(seq))
        name = line.split()[0][1:]; seq = []
    else: seq.append(line)
if name: lengths[name] = len(''.join(seq))

unannot = [g.strip() for g in subprocess.run(
    ['awk','-F','\t','NR>1 && (\$5==\"-\"||\$5==\"\") && (\$8==\"-\"||\$8==\"\") && (\$12==\"-\"||\$12==\"\")','{print \$1}'],
    input=open('$F').read(), capture_output=True, text=True).stdout.split('\n') if g.strip()]
lens = [lengths.get(g,0) for g in unannot]
print(f'Unannotated: {len(lens)} genes, avg {sum(lens)/len(lens):.0f} aa')
"
```

### Expected coverage (plant genome, ~30k genes)

| Metric | Typical |
|--------|---------|
| eggNOG GO | 40-50% |
| eggNOG KEGG | 20-25% |
| Diamond (any DB) | 85-95% |
| InterPro (any DB) | 70-85% |
| At least one strategy | 88-95% |
| No annotation | 5-12% (mostly short proteins <100aa) |

---

## Verification Checklist

- [ ] All input proteins present in merged output (0 missing)
- [ ] eggNOG annotations file has data rows (not just header)
- [ ] Diamond hit counts match between original blastp and merged
- [ ] InterProScan per-DB coverage: Pfam >80%, CDD >80%, SMART >85%
- [ ] No duplicate Gene_IDs in merged output
- [ ] Unannotated genes are mostly short (<100aa)
- [ ] SwissProt preferred as Diamond best hit when available

---

## Project Layout

```
<sample>/04.function_annotation/
├── eggNOG/
│   └── {sample}_eggnog.emapper.{annotations,hits,seed_orthologs}
├── diamond/
│   ├── {sample}_vs_{swissprot,trembl,nr,arabidopsis}.blastp
│   └── {db}.dmnd
├── interproscan/
│   ├── {sample}_interproscan.tsv
│   └── {sample}_interproscan.gff3
├── scripts/
│   ├── run_eggnog.sh
│   ├── run_diamond.sh
│   └── run_interproscan.sh
├── logs/
│   ├── eggnog_{sample}.log
│   ├── diamond_{sample}.log
│   └── interproscan_{sample}.log
└── results/
    └── {sample}_functional_annotation.tsv
```

---

## Citations

- eggNOG-mapper: Huerta-Cepas et al. (2019) "Fast Genome-Wide Functional Annotation through Orthology Assignment by eggNOG-Mapper." *Molecular Biology and Evolution* 34(8):2115-2122. doi:10.1093/molbev/msz099
- Diamond: Buchfink et al. (2021) "Fast and sensitive protein alignment using DIAMOND." *Nature Methods* 12:59-60. doi:10.1038/nmeth.3176
- InterProScan: Jones et al. (2014) "InterProScan 5: genome-scale protein function classification." *Bioinformatics* 30(9):1236-1240. doi:10.1093/bioinformatics/btu031

### Methods Template (Chinese)

> 蛋白编码基因的功能注释采用三种策略：（1）eggNOG-mapper v2.1.13（Huerta-Cepas et al., 2019）与eggNOG直系同源基因数据库比对，注释GO、KEGG和COG功能分类，分类群范围设为Viridiplantae（tax_scope 33090）；（2）DIAMOND v2.x（Buchfink et al., 2021）将蛋白序列与Swiss-Prot、TrEMBL、NR和拟南芥蛋白数据库进行blastp比对（一致性≥30%，E-value≤1e-5），选取最佳比对结果；（3）InterProScan v5.78（Jones et al., 2014）对蛋白序列进行结构域和motif搜索，涵盖Pfam、CDD、SMART、PRINTS、PANTHER、Gene3D、SUPERFAMILY等子数据库。三种策略的结果合并为统一注释表。
