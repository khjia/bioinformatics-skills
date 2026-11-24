# Plant Functional Annotation Skill

A Claude skill for reproducible protein-coding gene functional annotation in plant genome projects. It combines three complementary methods:

1. eggNOG-mapper for orthology-based functional annotations, including GO and KEGG when available.
2. DIAMOND for protein similarity searches against curated, broad, and taxonomic databases.
3. InterProScan for conserved signatures, domains, InterPro accessions, GO terms, and pathway evidence.

The workflow favors one representative protein per gene, preserves source IDs and raw outputs, records software/database versions, and reports per-method coverage rather than merging evidence into unqualified claims.

## Install as a Claude skill

Copy this folder into the relevant skills directory, for example:

```text
~/.claude/skills/plant-functional-annotation/
```

Keep `SKILL.md` at the skill folder root and retain the `scripts/` helper directory.

## Included helper

`scripts/select_best_diamond_hits.py` filters a DIAMOND tabular output made with the field order documented in `SKILL.md`, applies strict identity and E-value thresholds, and returns the best passing hit per query.

## Database and compute requirements

The skill does not bundle software or reference databases. Install or locate eggNOG-mapper, DIAMOND, and InterProScan on the analysis host and point each tool at a suitable local database. Database downloads can be large; review size, release, storage, licensing, and network constraints before obtaining them.
