---
name: plant-functional-annotation
description: Annotate plant protein-coding genes with eggNOG-mapper, DIAMOND homology searches, and InterProScan domain/signature analyses. Use for plant genome annotation, proteome functional annotation, GO/KEGG assignment, or integrating orthology, sequence-similarity, and conserved-domain evidence; supports haplotype-resolved assemblies and Chinese or English requests.
---

# Plant protein functional annotation

Run reproducible functional annotation for a protein-coding gene set using three complementary evidence sources: orthology, sequence similarity, and conserved signatures/domains. Keep the original gene or protein IDs in every output and preserve each method's raw results.

## Inputs and defaults

- Start from the current release's protein FASTA and matching GFF3. Prefer a one-protein-per-gene representative set such as `*.protein.longest.fa`; use the full isoform set only when requested.
- Do not infer that a file named `*.cds.fa` or `*.transcript.fa` contains the expected sequence type. Check the first few sequences and lengths, confirm IDs against the GFF3, and translate representative CDSs against the protein FASTA when sequence labels are uncertain.
- Keep samples separate. For a phased/diploid individual, a combined A+B proteome is suitable for functional annotation when IDs are unique and haplotypes remain identifiable. Split by haplotype only when requested or when the analysis requires haplotype-specific reporting. Retain organellar proteins if they are part of the gene set, and record whether they were included.
- Before running tools, inventory executable versions, local database versions and paths, free space, CPU/memory availability, input record counts, duplicate IDs, and any missing sequences. Prefer databases already installed on the analysis host. Do not download multi-gigabyte databases, install software, or replace existing databases without authorization.
- If user-level Python packages shadow the selected Conda environment, set `PYTHONNOUSERSITE=1` for that command rather than modifying unrelated environments.

## 1. Orthology and functional terms with eggNOG-mapper

Use protein mode and a local eggNOG database directory. The search stage commonly uses DIAMOND; choose sensitivity and CPU settings to match the compute budget.

```bash
emapper.py \
  -i "$PROTEINS" --itype proteins -m diamond --sensmode sensitive \
  --data_dir "$EGGNOG_DATA" --cpu "$THREADS" \
  --output "$SAMPLE" --output_dir "$OUT/eggnog"
```

Keep the `.emapper.annotations` table and run log. Report GO, KEGG KO/pathway, ortholog group, and description fields when present. A blank field means that method did not assign that annotation; it does not establish that the gene lacks the function.

## 2. Sequence similarity with DIAMOND

Search against the databases requested by the user. A common plant panel is Swiss-Prot, TrEMBL, NCBI NR, a model species (for example *Arabidopsis thaliana*), and a relevant taxonomic protein set. Keep each database as a separate result so source and hit counts remain clear.

Request several top alignments, then select the best hit **among alignments that pass the specified filters**. The default thresholds below implement the common requirement of identity strictly greater than 30% and E-value strictly less than 1e-5. No query-coverage filter is implied unless the user specifies one.

```bash
diamond blastp \
  --query "$PROTEINS" --db "$DIAMOND_DB" --out "$RAW_TSV" \
  --evalue 1e-5 --id 30 --max-target-seqs 25 --sensitive \
  --threads "$THREADS" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen \
    qstart qend sstart send evalue bitscore qlen slen stitle
```

Then apply strict `pident > 30` and `evalue < 1e-5`, select one passing hit per query by highest bitscore (ties: lower E-value, then higher identity), and retain the raw candidates. Use the bundled helper:

```bash
python scripts/select_best_diamond_hits.py \
  --input "$RAW_TSV" --output "$BEST_TSV" --database "$DATABASE_NAME"
```

If a different threshold is requested, pass `--min-identity` and `--max-evalue`; the script applies strict comparisons. Record database release/date and sequence count when available. Database titles are evidence labels, not curated functional conclusions.

## 3. Conserved signatures and domains with InterProScan

Use the installed member-database analyses, including the available Pfam, PRINTS, SMART, PANTHER, CDD, and other supported signatures. Check the installed InterProScan application list and data completeness first; do not assume every member database is installed merely because it is named in a methods draft.

```bash
interproscan.sh \
  -i "$PROTEINS" -b "$OUT/interproscan/$SAMPLE" -f TSV \
  -goterms -iprlookup -pa -cpu "$THREADS" -dp
```

`-dp` disables pre-calculated match lookup so results are generated locally. Preserve the TSV and run log. Report InterPro accessions, member database, signature accession, match coordinates, GO terms, and pathway fields when present.

## Run management and evidence integration

- Put outputs under a per-sample directory with separate `eggnog/`, `diamond/raw_hits/`, `diamond/best_hits/`, `interproscan/`, and `logs/` subdirectories. Store the exact commands, software versions, database versions, thresholds, input FASTA checksum, and run date in a README or manifest.
- On shared servers, check the current load and choose a conservative thread count. Prefer one long-running job per method or serialize runs when memory or database I/O is limiting. Capture stdout and stderr to per-sample logs; for jobs longer than a few minutes, use the site's job scheduler or a durable background job and inspect every terminal status.
- Validate that result IDs map to the input proteins, no duplicate query IDs were introduced, and each method completed successfully. Summarize per sample: total input genes, genes with an eggNOG ortholog, genes with GO/KEGG terms, genes with a passing DIAMOND hit per database, and genes with InterPro/domain/GO/pathway assignments.
- Produce an integrated table with one row per gene and separate columns for each evidence source. Preserve raw outputs; do not overwrite the source GFF or protein FASTA. Do not collapse conflicting functions into one unsupported label. Mark evidence source and missing values explicitly.
- In haplotype-resolved genomes, homologous genes can both receive annotations and database-dependent best hits can differ. Do not interpret duplicate annotations as erroneous duplication or as proof that both copies have identical function.
- Distinguish annotation coverage from biological validation. Homology and domain matches support functional hypotheses; they do not prove expression, activity, or causality. State database releases, thresholds, and incomplete or failed analyses in the report.
