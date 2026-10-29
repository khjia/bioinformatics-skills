# Known-locus check for plants (ortholog + literature)

Companion to `templates/26_known_loci.py`. There is no plant equivalent of the
human GWAS Catalog, so the practical check is: map the candidate to its
Arabidopsis/rice ortholog, then search the literature for that gene and the
trait.

## What the script does

`26_known_loci.py` queries **Europe PMC** with `<gene> AND (<trait keywords>)`,
sorts by citation count, and writes a TSV of papers plus an optional metadata
row (gene, trait, hit count, source).

```bash
python3 templates/26_known_loci.py \
  --gene AT5G43760 \
  --trait "fatty acid seed oil" \
  --out known.tsv --meta known_meta.tsv --max 20
```

Output columns: `pmid, doi, year, title, journal, first_author, cited_by`.

## Plant sources worth using

| Source | Use | Access |
|---|---|---|
| Europe PMC | literature by ortholog + trait | REST, no key |
| AraPheno | Arabidopsis phenotype/GWAS records | web/API, phenotype-centric |
| 1001 Genomes | Arabidopsis natural variation / GWAS | web portal |
| Ensembl Plants | orthologs, gene models, sequences | BioMart / REST |
| Gramene | comparative grass genomics, QTL | web |
| Brassica / rice QTL literature | seed oil / fatty acid QTL | Europe PMC |

## How to read the result

- A hit means "this gene family has been studied for this trait in a model
  species", **not** that the locus is validated in *Acer truncatum*.
- Ortholog mapping is imperfect: tandem KCS/FAE1 arrays often have
  species-specific copy numbers; one Arabidopsis gene may map to several
  paralogs.
- Trait mismatch matters: seed oil content, fatty-acid composition and VLCFA
  elongation are related but distinct.
- Absence of papers is weak evidence of novelty; negative results are rarely
  published.

## Recommended wording in a manuscript

"The candidate region contains a cluster of KCS/FAE1-like elongases
(orthologous to Arabidopsis AT5G43760); members of this family have been
associated with very-long-chain fatty acid and seed-oil traits in Arabidopsis
and Brassica (refs). No association has been reported for *Acer truncatum*."
