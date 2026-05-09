# bioinformatics-skills

Reusable bioinformatics workflow skills — structured runbooks for Hermes Agent / Claude / any LLM-driven dev environment that supports the "skills" pattern (a `SKILL.md` with frontmatter + optional `templates/`, `scripts/`, `references/`).

Each skill captures an end-to-end pipeline with:
- when to use / when not to use
- canonical project layout
- step-by-step commands (real, tested)
- parameters rationale table
- common pitfalls and their fixes
- verification checklist
- reusable templates

## Skills

### gwas-rmvp-htcondor
Run multi-trait multi-model GWAS (GLM + MLM + FarmCPU) with rMVP on an HTCondor cluster. Handles:
- PLINK BED → rMVP `big.matrix` conversion
- VanRaden kinship with caching
- External PC covariates (from LD-pruned PCA)
- Per-trait HTCondor batching with auto-resubmit on eviction
- **Adaptive λ_GC-driven PC tuning**: sweep PC ∈ {0,1,2,3,5,7,10}, pick the count that lands each (trait, model) inside the λ_GC [0.85, 1.15] band
- **Hardlink-based `final_results/` consolidation** to save hundreds of GB of intermediate sweep output
- Python post-processing: Manhattan + QQ plots (English-only titles to sidestep CJK font issues) + top-hits summary tables

Tested on common-bean (Phaseolus vulgaris) genotype panel, ~16M SNPs × 140–225 samples, 26–29 traits × 3 models per run.

## Layout

```
bioinformatics-skills/
├── README.md
├── LICENSE
└── <skill-name>/
    ├── SKILL.md
    └── templates/
        └── ...
```

## Usage

### With Hermes Agent

Drop the skill directory into your Hermes skills path:
```bash
cp -r gwas-rmvp-htcondor ~/.hermes/skills/bioinformatics/
```

Hermes will auto-index it on next startup. Call `skill_view(name='gwas-rmvp-htcondor')` to load the full content.

### Standalone / manual

Each `SKILL.md` is a standalone runbook — read top-to-bottom, adapt the project paths in `templates/*.py` / `.R` / `.sh`, run.

## License

MIT — see `LICENSE`. Free to use, modify, redistribute.

## Contributing

Found an edge case or a better parameter value? Open an issue or PR. The pitfalls sections grow best when people share what bit them.
