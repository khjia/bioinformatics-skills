# Known-locus lookup (GWAS Catalog and friends)

Companion to `templates/26_known_loci.py`. Distilled from the davila7
`gwas-database` skill (GWAS Catalog REST) and the ClawBio `gwas-lookup`
concept (federated variant query).

## When to use it

Before writing "novel locus", check whether the variant, gene or region has
already been reported. This is a literature/annotation step — it does **not**
run an association test.

## GWAS Catalog REST (human)

```bash
# by rsID
python3 templates/26_known_loci.py --rsid rs7903146 --out known_rs.tsv

# by gene symbol
python3 templates/26_known_loci.py --gene APOE --out known_apoe.tsv
```

Endpoints (no authentication, ~0.2 s between requests):

- `https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={rsid}`
- `https://www.ebi.ac.uk/gwas/rest/api/v2/associations?mapped_gene={SYMBOL}`
- Summary statistics: `https://www.ebi.ac.uk/gwas/summary-statistics/api/...`

Output columns: `rsid, traits, genes, pvalue, or_beta, risk_allele, study`.

## Non-model plants (this project)

The Catalog is human-only. Practical workflow:

1. Take the candidate gene (e.g. `Actnc04aG0088000`).
2. Map to its Arabidopsis ortholog (already in the annotation file:
   `AT5G43760`, KCS19).
3. Query the Catalog with the ortholog symbol, or search literature
   (PubMed/Europe PMC) for the gene family + trait (e.g. "KCS VLCFA seed oil").
4. Treat any hit as a **hypothesis**, not validation: ortholog, trait and
   population differ.

For cross-species variant lookup, Ensembl Plants / Gramene and the
Arabidopsis 1001 Genomes portal are more useful than the human Catalog.

## Federated lookup (optional)

ClawBio's `gwas-lookup` queries gnomAD, ClinVar, Open Targets, GTEx, LDlink
and others for a single rsID. Useful for human data; for plants only GTEx-like
or expression atlases apply, and most do not cover Acer truncatum.

## Interpretation cautions

- **Trait mismatch**: a lipid GWAS hit for LDL cholesterol is not evidence for
  seed fatty-acid composition.
- **LD, not causality**: Catalog associations are lead SNPs; the causal variant
  may be a proxy.
- **Ancestry**: effect sizes and even lead variants differ across ancestries.
- **Publication bias**: absence from the Catalog is not evidence of novelty.
