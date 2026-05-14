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
Run multi-trait multi-model GWAS (GLM + MLM + FarmCPU + BLINK) with rMVP on an HTCondor cluster. Handles:
- PLINK BED → rMVP `big.matrix` conversion
- VanRaden kinship with **fingerprint cache validation** (auto-invalidates on genotype change)
- External PC covariates (from LD-pruned PCA)
- Per-trait HTCondor batching with **trait-level retry** (`RMVP_RETRY`) and auto-resubmit on eviction
- **Adaptive λ_GC-driven PC tuning**: sweep PC ∈ {0,1,2,3,5,7,10}, pick the count that lands each (trait, model) inside the λ_GC [0.85, 1.15] band
- **Hardlink-based `final_results/` consolidation** to save hundreds of GB of intermediate sweep output
- **Multi-model consensus**: SNPs significant in ≥ N models → `high_confidence_loci.tsv`
- Python post-processing: per-trait Manhattan + QQ, multi-trait stacked Manhattan, trait × chrom hotspot heatmaps (English-only titles to sidestep CJK font issues)
- **Single-file HTML report** (base64-embedded PNGs, stdlib + pandas only)

Tested on common-bean (Phaseolus vulgaris) genotype panel, ~16M SNPs × 140–225 samples, 26–29 traits × 4 models per run.

**Citing**: any paper using this skill must cite rMVP (Yin et al., 2021) plus the original method paper for each model reported (Yu 2006 for MLM, Liu 2016 for FarmCPU, Huang 2019 for BLINK, VanRaden 2008 for K). Full reference list + BibTeX + a suggested Methods-section sentence are in `gwas-rmvp-htcondor/SKILL.md` under **Citations**.

**中文使用说明**: see `gwas-rmvp-htcondor/README_zh.md` for a Chinese quickstart covering inputs, the 5-step pipeline, adaptive PC tuning, multi-model consensus, common pitfalls, and the verification checklist.

### wgcna-tf-network
WGCNA module TF-structural gene co-expression network visualization. Pipeline:
- iTAK-based transcription factor annotation (hmmscan + TF_Rule classification)
- Per-pathway ring-layout network plots: inner ring = TF, outer ring = structural genes
- Edge width proportional to WGCNA co-expression weight; TF→structural edges bolded
- NPG (Nature Publishing Group) color palette

Includes reusable scripts: `itak_classify.py` (standalone TF classifier) and `plot_tf_network.py` (configurable network plotter).

**Citing**: iTAK (Zheng et al., 2016), WGCNA (Langfelder & Horvath, 2008), HMMER (Eddy, 2011). Full BibTeX in `wgcna-tf-network/SKILL.md`.

**中文使用说明**: 见 `wgcna-tf-network/SKILL.md` 底部中文章节。

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
