---
name: genome-circos
description: Build a publication-style, two-assembly Circos overview with paired chromosome sectors, gene and repeat density, telomeres, rDNA, optional alignment ribbons, and an English figure legend.
argument-hint: <manifest.json> <new-output-directory>
disable-model-invocation: true
---

# Paired genome Circos

Build an interleaved circular comparison of exactly two chromosome-level assemblies. The first assembly is darker with a solid outline; the second is lighter with a dashed outline. Both retain the chromosome-pair order supplied by the user.

1. Read [the input contract](references/inputs.md). Use the matching `.fai`, final filtered gene GFF3, and classified repeat GFF3 for each assembly. Ask for the exact versions or chromosome pairing if they are ambiguous; do not silently pick an older annotation or pair chromosomes by similarity of names.
2. Verify the assembly order, nuclear chromosome list, and PAF direction. Use reviewed marker tables if available. Treat unsupplied telomere calls as not assessed, and omit unsupplied rDNA marks or PAF ribbons. Do not infer telomere status from a chromosome's name, length, or the presence of a row in an upstream report.
3. Run `python -I "${CLAUDE_SKILL_DIR}/scripts/draw_circos.py" "<manifest.json>" "<new-output-directory>"`. Add `--export` when PDF and PNG are requested and PyMuPDF is available; SVG and the English legend need only Python 3.10+. Do not install dependencies or overwrite an existing output directory without the user's authorization.
4. Inspect the exported figure, English legend (`legend.en.md`), and raw window values (`window_values.tsv`). Check sector order, label placement, optional inputs, and the reported number of ribbons. The heights are scaled independently per chromosome sector and track; they do not compare absolute density between assemblies.
5. Deliver the figure and English legend together. Ribbons are input sequence alignments, not SyRI calls. Do not infer inversions, translocations, or causal phenotypes from this figure alone. Only share data or outputs externally if the user explicitly asks to share those items.

For a local smoke test, generate synthetic data with `python -I "${CLAUDE_SKILL_DIR}/scripts/make_demo.py" "<new-demo-directory>"` and use its `manifest.json`. No real genomes or manuscript figures are bundled.
