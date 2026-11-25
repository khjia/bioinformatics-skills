# Input contract and methods

Create a UTF-8 JSON manifest. Paths are relative to the manifest unless absolute. The first assembly is darker and drawn first clockwise within each pair.

```json
{
  "assemblies": [
    {"label":"Assembly A", "fai":"A.fa.fai", "genes":"A.final.gff3", "repeats":"A.EDTA.gff3", "telomeres":"A.ends.tsv", "rdna":"A.rdna.tsv"},
    {"label":"Assembly B", "fai":"B.fa.fai", "genes":"B.final.gff3", "repeats":"B.EDTA.gff3"}
  ],
  "pairs":"pairs.tsv",
  "paf":"B_to_A.paf",
  "window_bp":1000000,
  "min_link_bp":10000
}
```

- `fai`: contig IDs and lengths in the first two FAI columns. Extra contigs are not displayed unless selected in `pairs`.
- `genes`: final, filtered GFF3 from the same assembly. Only `gene` records count. The lower of the two central positions of an even-length interval is its midpoint. GFF3 is 1-based inclusive and is converted before binning. End windows are truncated to the chromosome length; values are per window, not divided by its size.
- `repeats`: classified GFF3 (e.g. EDTA), with `Classification` or `classification` attributes. The class rules are LTR/Gypsy, LTR/Copia, DNA/* or MITE/* or RC/Helitron, and LINE/* or LINE or Unknown or LTR/unknown or SINE/unknown. Other labels are not shown. Repeat values are union-covered base pairs per bin **within each class**, not feature counts; overlapping classes can cover the same bases. This implementation does not mask rDNA or reproduce project-specific repeat subfamily filters from the originating figure. If masking is needed, document the preprocessing and supply the reviewed annotations; update the English legend accordingly.
- `pairs`: tab-separated header `label`, `first`, `second`. One row per chromosome pair in display order; `first` and `second` must match their respective FAI/GFF3 IDs. `label` is the shared outer label. Include only the intended chromosomes (e.g. omit organelles from a nuclear genome overview). Gray is used as a neutral pair accent after the first 11 positions; labels still identify those pairs.
- `telomeres` (optional): TSV header `chromosome`, `end`, `status`. One row per assessed end. End is `left` (low coordinates) or `right` (high coordinates). Status is `detected`, `not_detected`, or `not_assessed`. Omitted ends are gray (not assessed). Caps are schematic and do not encode length.
- `rdna` (optional): TSV header `chromosome`, `start`, `end`, `type`; type is `5S` or `18S`. Supply reviewed intervals in 0-based, half-open coordinates, not raw BLAST coordinates. Only chromosomes present in the pairing table may appear in marker inputs. Markers are widened to a minimum visible angular size when needed.
- `paf` (optional): standard PAF, with the **first assembly as target** and **second as query**, e.g. from minimap2 with that reference/query order. Show only explicitly paired chromosomes, discard records explicitly tagged non-primary, accept untagged records, and filter on alignment span in **both** genomes. Reverse-strand ribbons are twisted but are not inversion calls. Record the upstream aligner and filtering settings in the manuscript methods.
- Files above may be gzip-compressed using a `.gz` suffix, except the JSON manifest. A missing optional file key (or a `null` value) means that input was not supplied.

Output `window_values.tsv` retains raw values and 0-based, half-open coordinates. Retain the manifest and provenance of all input datasets with your analysis; do not publish them by default.
