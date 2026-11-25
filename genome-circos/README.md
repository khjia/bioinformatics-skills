# Genome Circos

Compare two chromosome-level assemblies in one circular plot. Chromosome pairs are interleaved clockwise: the first assembly uses darker fills and solid outlines, the second lighter fills and dashed outlines. Rings show gene and repeat density, telomere status, and rDNA markers. Optional ribbons show alignments between paired chromosomes, not structural variant calls.

## Install as a Claude Code skill

Copy this entire `genome-circos` directory to `~/.claude/skills/` (personal) or `.claude/skills/` (project), including its scripts and references. Invoke `/genome-circos` with the path to your manifest and a new output directory.

## Run the synthetic example

From this directory, with Python 3.10 or later:

```powershell
python -I scripts/make_demo.py demo-data
```

```powershell
python -I scripts/draw_circos.py demo-data/manifest.json demo-output
```

The outputs are `circos.svg`, `legend.en.md` (English figure legend), and `window_values.tsv` (raw counts and coverage). Add `--export` to produce a vector PDF and a 300-dpi PNG if [PyMuPDF](https://pymupdf.readthedocs.io/) is installed. The SVG, legend, and TSV need only the Python standard library. Re-runs require a new output directory.

For your own assemblies, follow the [manifest and coordinate specification](references/inputs.md). Bar heights are normalized separately for each track and chromosome sector; use the raw TSV for between-assembly comparisons.

Run the checks with:

```powershell
python -I -m unittest discover -s tests
```
