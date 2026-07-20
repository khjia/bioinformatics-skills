---
name: bio-alternative-splicing
description: "Complete alternative splicing analysis: PSI quantification, QC metrics, differential splicing (rMATS, leafcutter, MAJIQ, SUPPA2), isoform switching with functional consequences, and splice variant prediction. Covers SE/A5SS/A3SS/MXE/RI event types."
---


## Algorithmic Taxonomy


| Family | Unit | Reference tools | Fails when |
|--------|------|-----------------|------------|
| Event-based | Pre-defined SE/A5SS/A3SS/MXE/RI events from annotation | rMATS-turbo, SUPPA2, VAST-TOOLS | Event isn't in annotation; complex multi-junction events split arbitrarily; AFE/ALE confounded with splicing |
| LSV-based | Local Splice Variations at single source/target nodes | MAJIQ V3 | Memory-constrained environments; cohorts smaller than ~3 reps; non-academic users (license) |
| Junction-cluster | Annotation-free intron clusters by shared splice sites | leafcutter, leafcutter2 | Undersampled clusters lose power; topology biologically uninterpretable for novel events |
| Splice-graph | Graph nodes (non-overlapping exonic regions) | Whippet, Shiba | Whippet maintenance status uncertain since ~2022; complex multi-exon graphs |
| Coverage-aware (IR) | Intron body coverage + flanking junctions | IRFinder-S, S-IRFindeR, iREAD | Confounded by overlapping exons, repeats, low mappability regions |
| Isoform-based | Transcript abundance via EM | Salmon/kallisto + tximport | Salmon EM uncertainty propagates; many similar isoforms (TTN, MAPT) become indistinguishable |

The agent's first decision is **which family** the question requires, not which tool. Switching family is the response to a tool failing within a family — switching tool within a family rarely fixes systematic blind spots.


## Event Taxonomy (Beyond Standard SE/A5SS/A3SS/MXE/RI)


| Class | Code | Biology | Detection caveat |
|-------|------|---------|-------------------|
| Skipped exon (cassette) | SE | Default cassette-exon AS | Most common; well-handled by all tools |
| Alternative 5' splice site | A5SS | Alternative donor; intron 5' end varies | Sign convention tool-specific (see below) |
| Alternative 3' splice site | A3SS | Alternative acceptor; intron 3' end varies | Sensitive to BPS cancer mutations (SF3B1) — cryptic 3'ss ~10-30nt upstream |
| Mutually exclusive exons | MXE | Two exons paired, one included | Tool implementations vary; verify which "form 1" is yours |
| Retained intron | RI | Whole intron retained in mature mRNA | Junction-only quant systematically underdetects; needs IRFinder-S |
| **Microexon** | (sub-SE) | 3-27 nt; neural-enriched, SRRM4-regulated | Missed by default aligner anchor lengths (>=20-30 nt); needs VAST-TOOLS, MicroExonator, or long-read |
| **Exitron** | n/a | Intronic region within an annotated CDS exon | Mis-classified as A5SS/A3SS by most tools; use ExitronFinder, ScanExitron |
| **Alternative first exon (AFE)** | AFE | Alternative TSS / promoter use | Promoter-driven, NOT spliceosomal; confirm with FANTOM CAGE before reporting as splicing |
| **Alternative last exon (ALE)** | ALE | Alternative cleavage/polyadenylation | APA-driven, NOT spliceosomal; confirm with 3'-end-seq |
| **Detained intron (DI)** | (RI subtype) | Nuclear-retained on mature mRNA, regulated by Clk kinases | Distinct from cytoplasmic NMD-targeted RI (Boutz 2015 *Genes Dev*); requires fractionation to confirm |
| **Recursive splicing** | (intron subtype) | Long introns >50kb spliced via internal ratchet points | Sibley 2015 *Nature*; only detectable with long-read or nascent RNA-seq |

Tool-agnostic taxonomy reference: Wang 2008 *Nature*; Vaquero-Garcia 2016 *eLife* (LSV); Tapial 2017 *Genome Res* (VastDB).


## Tool Selection Matrix


| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| rMATS-turbo | Standard SE/A5SS/A3SS/MXE/RI in annotated organism, n>=3 | BAM + GTF | Fast, well-calibrated at n>=3, novel SS support | Junction read imbalance; novel multi-junction events; underdetects RI and microexons |
| SUPPA2 | Quick PSI from existing TPM, pilot analysis | Salmon/kallisto TPM + GTF | No alignment; fastest | Annotation-bound; high FDR (15-30%) at n<=2 vs n<=2 |
| MAJIQ V3 | Complex events, heterogeneous cohorts (HET) | BAM + GFF3 | Bayesian posterior PSI, complete LSV semantics | High memory (~50+ GB on cohorts); academic license; complex LSV interpretation needs care |
| leafcutter | Novel junction discovery, sQTL, low-memory environments | BAM (regtools junctions) | Annotation-free; ~400 MB memory; SOTA for unannotated organisms | Sensitive to read depth; cluster topology arbitrary for complex multi-junction events |
| Shiba (2025) | Low-coverage / few-replicate designs; junction imbalance correction | BAM + GTF | Best calibration in own benchmark at n=2 vs n=2 | New (2025); limited community calibration |
| VAST-TOOLS | Cross-species comparative AS, microexons | FASTQ | VastDB orthology, ExOrthist co-tool | Limited to species in VastDB |
| Whippet | Laptop-scale exploratory | FASTQ | Fast splice-graph-based PSI | Reduced active development since ~2022; underperforms on complex topologies |
| IRFinder-S | Intron retention specifically | FASTQ | Coverage + junction integration; CNN-based artifact filtering | IR-only; not for cassette events |
| S-IRFindeR | Replicate-stable IR ratio | BAM | Stable IR ratio metric | IR-only; less integrated than IRFinder-S |

Methodology evolves; verify benchmarks (Olofsson 2023 *Brief Bioinform*; Kubota 2025 *NAR*; Tran 2025 *WIREs RNA*) and tool docs before committing. Default 2026 recommendation: run rMATS-turbo + leafcutter and reconcile; add MAJIQ V3 for complex events / heterogeneous cohorts; switch to Shiba for n=2 vs n=2.


## PSI Definition and Effective Length Normalization


For a cassette exon, naive PSI ignores that the inclusion isoform contains more positions where a junction read can map than the skipping isoform. rMATS reports `IncFormLen` (= 2*(read_length - anchor) for the two flanking junctions, plus exon body bases) and `SkipFormLen` (= read_length - anchor) and computes:

```
PSI = (IJC / IncFormLen) / (IJC / IncFormLen + SJC / SkipFormLen)
```

where IJC = inclusion junction counts, SJC = skipping junction counts. **Skipping this normalization biases PSI by ~10-30% depending on read length and exon size.** SUPPA2 derives PSI from transcript TPMs, which the upstream Salmon/kallisto already accounts for. leafcutter operates on intron usage proportions within a cluster (a different statistic).

For long-read data, every read carries full isoform identity — effective-length normalization becomes unnecessary because each read counts as one isoform.


## Sign Conventions for Alternative Splice Sites


| Tool | A5SS interpretation | A3SS interpretation | "Inclusion" direction |
|------|----------------------|----------------------|------------------------|
| rMATS | "long" form = donor downstream of alternative donor | "long" form = acceptor upstream of alternative acceptor | PSI > 0 = more long form |
| SUPPA2 | Same as rMATS (long = inclusion of additional exon body) | Same | PSI > 0 = more long form |
| VAST-TOOLS | Encoded in event ID (`_D1` vs `_D2`) | Encoded in event ID (`_A1` vs `_A2`) | Document the chosen reference |
| MAJIQ | Per-junction within LSV; explicit donor/acceptor naming in VOILA | Same | PSI per junction in the LSV |

Always record which alternative form ΔPSI > 0 corresponds to in publication-grade reporting. Confusion is the most common reviewer comment for AS papers.


## rMATS-turbo Workflow


**Goal:** Quantify SE/A5SS/A3SS/MXE/RI events from BAMs aligned with STAR 2-pass.

**Approach:** Group BAMs by condition, run rMATS with `--statoff` for quantification only, then parse JC.txt files for per-replicate PSI.

```bash
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --statoff
```

Key flags: `--novelSS` discovers junctions absent from the GTF (recommended with STAR 2-pass output). `--variable-read-length` allows mixed read lengths in the cohort. `--libType fr-firststrand` matches Illumina TruSeq stranded; verify with RSeQC `infer_experiment.py`. `--statoff` is for quantification-only runs; omit for differential testing.

```python
import pandas as pd

se_jc = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

inc_cols = [c for c in se_jc.columns if c.startswith('IncLevel')]
se_jc['mean_PSI'] = se_jc[inc_cols].mean(axis=1)

per_rep_inc = se_jc['IJC_SAMPLE_1'].str.split(',').apply(lambda x: list(map(int, x)))
per_rep_skip = se_jc['SJC_SAMPLE_1'].str.split(',').apply(lambda x: list(map(int, x)))
min_inc = per_rep_inc.apply(min)
min_skip = per_rep_skip.apply(min)

reliable = se_jc[(min_inc + min_skip) >= 20]
```

### JC vs JCEC files

`SE.MATS.JC.txt` uses **only junction-spanning reads**. `SE.MATS.JCEC.txt` adds **reads contained within the alternative exon body** as inclusion evidence.

- Prefer **JC** for clean cassette-exon analysis when reads spanning junctions are sufficient.
- Use **JCEC** when alternative exons are short (<50nt) and junction-spanning reads are scarce.
- Avoid **JCEC** when intron retention overlaps the alternative exon body — exon-body reads may come from retained introns, not inclusion isoform.


## SUPPA2 Workflow


**Goal:** Compute event PSI from transcript TPM without alignment; useful when Salmon/kallisto TPMs already exist.

**Approach:** Generate IOE event definitions from GTF, then aggregate TPMs of transcripts including/excluding each event.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI AF AL

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e transcript_tpm.tsv -o psi_${ev}
done
```

SUPPA2 is annotation-bound: events absent from the GTF cannot be quantified. Whether an event is detected depends entirely on which transcripts the upstream Salmon/kallisto index contains. Use GENCODE comprehensive over basic when SUPPA2 detection sensitivity matters.


## MAJIQ V3 Workflow


**Goal:** Quantify LSVs with Bayesian posterior PSI distributions; ideal for complex multi-junction events that don't fit canonical event types.

**Approach:** Build a splice graph from BAMs + GFF3, compute per-junction coverage with bootstrap, then run `majiq psi` for posterior PSI per LSV.

```bash
majiq build annotation.gff3 -c settings.ini -j 8 -o build_output
majiq psi build_output/sample1.majiq build_output/sample2.majiq -j 4 -o psi_output -n condition_psi
voila view -p 5000 -j 8 build_output/splicegraph.zarr psi_output/condition_psi.psi.voila -o voila_output
```

MAJIQ V3 (Slaff et al *bioRxiv* 2024; public release 2025) replaced V2's SQLite splicegraph (`splicegraph.sql`) with **Zarr storage** (`splicegraph.zarr`); the `.sql` is deprecated. V3 is ~3.2x faster than V2 via xarray/zarr/Dask parallelization. LSV output includes posterior mean PSI plus the full posterior distribution; this enables threshold-based testing (e.g. P(|ΔPSI| > 0.2)) rather than frequentist p-values.


## leafcutter Junction Quantification


**Goal:** Detect junctions and intron clusters annotation-free for downstream cluster-level usage.

**Approach:** Extract junctions per BAM with regtools, write filenames into a list, then cluster introns sharing splice sites.

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

`-a 8` = 8nt anchor minimum (raise to 12 for stricter; lower to 6 for microexon-friendly). `-m 50` = minimum junction reads per cluster. `-l 500000` = max intron length (relevant for long brain-gene introns; raise for genes like DSCAM, ROBO2, ANK3).


## Per-Tool Failure Modes


### rMATS-turbo: Junction Read Imbalance

**Trigger:** A cassette exon's flanking exons have unequal read mapping opportunity (very short upstream exon, repeat-overlapping flanks, or low-mappability regions).

**Mechanism:** rMATS' binomial model treats inclusion vs skipping junctions as having equal mappability. When mappability differs between the two junction types, the PSI estimate is biased.

**Symptom:** "Significant" rMATS calls with no concordant change in leafcutter or MAJIQ at the same locus; ΔPSI direction inconsistent with sashimi-plot intuition.

**Fix:** Run Shiba (Kubota 2025 *NAR*) which corrects junction-imbalance, or filter rMATS hits requiring concordant detection in leafcutter.

### SUPPA2: Sparse Empirical Null at Low Replicate Count

**Trigger:** n=2 vs n=2 (or n=3 vs n=2) design with `--method empirical`.

**Mechanism:** SUPPA2's empirical null is constructed from between-replicate ΔPSI distributions binned by transcript expression. With few replicates, the binned null is sparse and conservative-looking but actually under-calibrated.

**Symptom:** Inflated FDR (15-30% in benchmarks); many "significant" hits don't replicate or validate.

**Fix:** Switch to leafcutter or Shiba for n<=3 designs; or use `--method classical` (Wilcoxon) for very low replicate count; reconcile against orthogonal tool.

### MAJIQ V3: Complex LSV Interpretation

**Trigger:** A gene with 4+ alternative splice sites at one node (e.g. one source, multiple acceptors).

**Mechanism:** A complete LSV at a single node lists all observed junctions; PSI is per-junction within the LSV, not "PSI of one event."

**Symptom:** Reporting "PSI of the gene" doesn't make sense; per-junction PSIs sum to 1 across the LSV but no single number represents the gene.

**Fix:** Use VOILA to visualize the LSV graph and identify which junction(s) shifted; for cassette-style reporting, derive equivalent PSI from sum of inclusion-junctions / total junctions in the LSV.

### leafcutter: Cluster Topology Arbitrariness

**Trigger:** A cluster has 4+ introns sharing splice sites with non-canonical topology (e.g. mixed cassette + alternative donor + IR).

**Mechanism:** leafcutter clusters introns by shared splice sites; complex topologies don't map onto SE/A5SS/A3SS taxonomy and cluster-level "ΔPSI" hides which intron drove the change.

**Symptom:** Significant cluster-level p-value but multiple introns showing different effect-size directions.

**Fix:** Inspect the cluster in leafviz; report per-intron effect sizes (`effect_sizes.txt`); for canonical event reporting, map to SE/A5SS/A3SS via flanking exon coordinates manually.


## Reconciliation: When rMATS and leafcutter Disagree


The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction read imbalance OR rMATS event hits an annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Check `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | Event class mismatch (e.g. rMATS calls SE positive but leafcutter sees A5SS shift in same cluster) | Manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based).


## Intron Retention: Canonical vs Detained vs Co-Transcriptional Unspliced


Three biologically distinct states all called "IR" by generic tools:

1. **Canonical RI (cytoplasmic, NMD-substrate often)**: mature polyadenylated mRNA carries the intron; usually PTC-bearing and NMD-targeted, sometimes encoding an alternative protein.
2. **Detained intron (DI)** (Boutz 2015 *Genes Dev*): nuclear-localized, mature transcripts retaining a specific intron; a regulated reservoir released into translation upon signaling.
3. **Co-transcriptional unspliced**: nascent pre-mRNA captured before splicing complete; not a regulated state.

**Library prep determines which state(s) you see:**
- Poly(A) selection: enriches (1), depletes (2)/(3)
- rRNA depletion (cytoplasmic): captures (1)
- rRNA depletion (whole cell or nuclear): captures all three

**To distinguish DI from canonical RI:** subcellular fractionation (nuclear vs cytoplasmic RNA-seq), or NMD inhibitor (cycloheximide, NMDi-14) treatment — canonical RI mRNA increases under NMD inhibition; DI does not.

```bash
IRFinder -m FullAuto -r REF/ -d ir_output sample.fastq.gz
```

IRFinder-S (Lorenzi 2021 *Genome Biol*) uses CNN-based filtering of true IR vs noise; current SOTA for IR analysis. iREAD and S-IRFindeR (Broseus & Ritchie 2020 *bioRxiv*) are alternatives.


## Microexon Detection


Microexons (3-27 nt, neural-enriched, SRRM4-regulated; Irimia 2014 *Cell*) are missed by default short-read aligners requiring 20-30 nt anchors. Options:

| Approach | Tool | Notes |
|----------|------|-------|
| Curated database lookup | VAST-TOOLS + VastDB | Cross-species, microexon-aware (Tapial 2017 *Genome Res*) |
| De novo discovery | MicroExonator (Parada 2021 *Genome Biol*) | Snakemake pipeline |
| Tune the upstream aligner | `STAR --alignSJoverhangMin 6 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04` | rMATS itself cannot recover microexons that STAR didn't pass through; lower DB-junction overhang to 1 (trusts annotated microexon coords) and combine with strict mismatch filter. Typical AS pipelines use STAR 8/3 which is too strict for microexons |
| Long-read sequencing | PacBio Iso-Seq, ONT | Solves the problem entirely; reads span microexons fully |

For brain / neural tissue or autism-spectrum studies, **microexon analysis must be explicit** — default short-read pipelines underdetect them by ~70%.


## Quality Thresholds


| Metric | Threshold | Source / Rationale |
|--------|-----------|---------------------|
| Junction reads per replicate | >=10-20 (per-replicate minimum) | Empirical PSI variance becomes <0.05 above this; below, PSI becomes a coin flip |
| PSI dynamic range | mean PSI 0.05-0.95 | Outside is near-constitutive; rMATS, SUPPA2 default filters drop these |
| Missing values | <50% of samples | Higher missingness indicates low expression — re-test with subset |
| Read length | >=75nt PE preferred; >=100nt for microexons | 50nt SE biases toward shorter exons (Wang 2008 *Nature*) |
| Library | rRNA depletion for IR analysis; poly(A) acceptable for cassette | Sims 2014 *Genome Res*; poly(A) loses pre-mRNA |
| STAR 2-pass | Cohort-style preferred over per-sample basic | Veeneman 2016 *Bioinformatics*: >=94% novel junction recovery |
| MAJIQ minreads / minpos | --minreads 10 --minpos 3 | Default; lower for low-coverage |
| leafcutter -m | 50 reads per cluster | Higher for rare events; lower for sQTL discovery |
| Anchor length | >=8 nt for short-read | Below this, false-positive junctions dominate (CIGAR-N noise) |


## Decision Tree by Scenario


| Scenario | Recommended tool(s) | Why |
|----------|----------------------|-----|
| Standard cassette analysis, n>=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Default workflow; complementary algorithmic families |
| Non-model organism, no GENCODE-grade annotation | leafcutter + de novo discovery | Annotation-free |
| Heterogeneous cohort, n>=10 vs n>=10 (clinical, GTEx-style) | MAJIQ V3 with HET module | HET designed for between-sample variability dominance |
| Low coverage / few replicates (n=2 vs n=2) | Shiba | Junction-imbalance correction; SOTA at low coverage in 2025 benchmarks |
| Cross-species comparative (vertebrate panel) | VAST-TOOLS + VastDB | Orthology-aware events; ExOrthist co-tool |
| TPM-only available (no BAMs) | SUPPA2 | Annotation-bound but fast |
| Microexon focus (neural / ASD) | VAST-TOOLS or MicroExonator | Default tools systematically miss microexons |
| Intron retention focus | IRFinder-S (rRNA-depleted library) | Coverage-aware; CNN artifact filter |
| Detained introns specifically | IRFinder-S + nuclear/cytoplasmic fractionation | Required to separate DI from cytoplasmic RI |
| Long reads available | rMATS-long, FLAIR, IsoQuant | Full-isoform resolution; see long-read-splicing |
| Single-cell (full-length plate) | MARVEL, BRIE2 | See single-cell-splicing |
| Single-cell (10X 3') | Likely don't attempt; consider Sierra for APA | 10X 3' chemistry insufficient for AS |


## Common Errors


| Error | Cause | Solution |
|-------|-------|----------|
| `error: GTF gene_id parsing` (rMATS) | rMATS expects GENCODE-style gene_id; some Ensembl GTFs use different attribute order | `gffread input.gff3 -T -o standardized.gtf` |
| `KeyError: 'IJC_SAMPLE_1'` (rMATS parsing) | Output column missing; sometimes occurs when --statoff combined with novel events on older versions | Update rMATS-turbo to >=4.3.x; re-run |
| `MAJIQ: too few reads at junction` | Default `--minreads 10 --minpos 3` filters out the locus | Lower thresholds for low-coverage data; document filtering |
| `leafcutter: dispersion estimation failed` | Cluster has all-zero counts in one group | Pre-filter clusters with `--min-samps-feature-prop 3` |
| `SUPPA2: empirical p computed on N=4 nulls` | Insufficient replicates for empirical mode | Switch `--method classical` (Wilcoxon) for very low replicate count |
| `regtools: invalid CIGAR` | Non-BAM-spec read in input | Filter with `samtools view -h -F 0x100 -F 0x800` (drop secondary/supplementary) |
| `STAR: too many SJs` (in pass 2) | Cohort SJ.out.tab too large | Filter to junctions seen in >=3 samples or with >=3 unique reads before merging |


## Output Interpretation


PSI ranges 0 to 1: 1 = always included, 0 = always skipped, 0.5 = balanced. Sign of `IncLevelDifference` matches `--b1 minus --b2` group order — always document which is which in publications.

**NMD direction matters:** an increase in PSI of a poison exon (PTC-introducing) **decreases** functional protein due to NMD. Always check whether the alternative form is PTC-bearing using ORF-aware annotation (IsoformSwitchAnalyzeR consequences, or manual stop-codon distance check vs last exon-exon junction).

**Disease signatures:**
- **SF3B1** mutations (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30 nt upstream of canonical (Darman 2015 *Cell Rep*). Look for clustered A3SS hits.
- **U2AF1** mutations (lung adeno, MDS): altered preferences at 3'ss -3 position; cassette-exon shifts.
- **TDP-43 loss** (ALS/FTD): de novo cryptic exons in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*) — annotation-free tools required (leafcutter denovo).


## Common Pitfalls


- Treating AFE/ALE as splicing — these are typically promoter-driven (AFE) or APA-driven (ALE), not spliceosomal. Confirm with FANTOM CAGE or 3'-end-seq.
- Confusing detained introns with NMD-targeted RI — both call as "IR" but have opposite biological fates.
- Using poly(A) libraries for IR analysis — biases toward mature transcripts, depletes pre-mRNA.
- Single-end short reads — junction-spanning reads need >=8nt overhang on both sides; biases toward shorter exons.
- Quoting "PSI of the gene" from MAJIQ LSV output — only per-junction PSI within an LSV is meaningful.
- Skipping STAR 2-pass — loses ~14% of novel junctions; matters for any non-canonical organism or condition.
- Trusting rMATS calls without `--novelSS` when STAR 2-pass found new junctions — rMATS will only quantify pre-annotated events.


## Related Skills


- differential-splicing - Compare PSI between conditions; use the same upstream alignment but switch to with-stat tools
- splicing-qc - Run BEFORE quantification to verify library, depth, strandedness, alignment quality
- isoform-switching - DTU framework with NMD/ORF/domain consequences; complementary to event-level PSI
- sashimi-plots - Visualize specific events for QC and reporting; concordance check across tools
- splice-variant-prediction - SpliceAI/Pangolin for variant impact predictions to test against PSI changes
- long-read-splicing - Full-isoform PSI without anchor-length limits; preferred for microexons and complex isoforms
- read-alignment/star-alignment - STAR 2-pass cohort-style alignment is required upstream
- rna-quantification/alignment-free-quant - Salmon/kallisto TPM is required for SUPPA2


## References


- Wang et al 2008 *Nature* - AS event taxonomy
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ LSV framework
- Trincado et al 2018 *Genome Biol* - SUPPA2
- Li et al 2018 *Nat Genet* - leafcutter
- Tapial et al 2017 *Genome Res* - VAST-TOOLS / VastDB
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Slaff et al 2024 *bioRxiv* - MAJIQ V3
- Kubota et al 2025 *NAR* - Shiba
- Lorenzi et al 2021 *Genome Biol* - IRFinder-S
- Boutz et al 2015 *Genes Dev* - detained introns
- Irimia et al 2014 *Cell* - SRRM4 microexons
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Olofsson et al 2023 *Brief Bioinform* - benchmark
- Tran et al 2025 *WIREs RNA* - methodology review
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing
- Veeneman et al 2016 *Bioinformatics* - STAR 2-pass benchmark


## QC Layer Taxonomy


| Layer | Target | Tool | Fails when |
|-------|--------|------|------------|
| Experimental design | Read length, depth, replicates, library type | Pre-sequencing review | <PE 75nt; n<3 vs n<3; <30M reads/sample |
| Library prep | poly(A) vs rRNA depletion | Pre-sequencing review | poly(A) library used for IR analysis |
| Alignment | STAR 2-pass cohort-style | STAR | 1-pass loses 14% novel junctions; per-sample 2-pass introduces inconsistency |
| Junction discovery | Saturation, novelty | RSeQC `junction_saturation`, `junction_annotation` | Curve still rising = under-sequenced; novel% >40% suggests biology or artifact |
| Strand specificity | Library protocol consistency | RSeQC `infer_experiment` | Wrong `--libType` halves usable junctions |
| Splice site strength | Cryptic vs canonical | MaxEntScan, SpliceAI | Weak splice sites (MaxEnt<5) may indicate cryptic, regulated, or annotation error |
| Junction overhang | Read-junction support quality | pysam CIGAR parsing | Overhang <8nt = high false-positive rate |
| Contamination | rRNA, adapters | fastq_screen | >20% rRNA in "depleted" library = failed depletion |
| Annotation | GENCODE basic vs comprehensive | Annotation choice | Basic for canonical events; comprehensive for DTU |


## Decision Tree by Question


| Question | Recommended QC |
|----------|-----------------|
| Will my planned RNA-seq design support AS analysis? | Pre-sequencing audit: library type, read length, depth, replicates |
| Is my data suitable for cassette exon analysis? | Junction saturation + known/novel ratio + read length |
| Why does my AS analysis call so few events? | Saturation curve, depth, library type, alignment 2-pass |
| Why does my AS analysis call so many novel junctions? | Annotation completeness + novel% + biology check (TDP-43, SF3B1) |
| Are my SpliceAI predictions calibrated for my tissue? | MaxEntScan + SpliceAI concordance for known sites |
| Did STAR 2-pass actually run cohort-style? | Verify SJ.out.tab merging across samples |
| Is intron retention detectable in my data? | Library type (must be rRNA-depleted); strand-specific |
| Are my microexons detectable? | Read length >=100; aligner anchor settings; consider VAST-TOOLS |


## Experimental Design Audit (Before Sequencing)


| Decision | For splicing analysis | Rationale |
|----------|------------------------|-----------|
| **Library prep** | rRNA depletion (Ribo-Zero, RiboCop) | poly(A) selection loses pre-mRNA, nascent transcripts, and detained introns; for IR analysis rRNA depletion is mandatory (Sims 2014 *Genome Res*) |
| **Read length** | PE 100-150 nt (PE 150 strongly preferred) | Junction-spanning reads need >=8 nt overhang on each exon; 50 nt SE biases toward shorter exons (Wang 2008 *Nature*) |
| **Pairing** | Paired-end | Single-end loses fragment-level disambiguation of junctions |
| **Depth** | 50-100M reads/sample | DGE-grade 30M misses low-PSI events; 100M for low-abundance event discovery |
| **Strandedness** | Stranded library (Illumina TruSeq stranded) | Distinguishes overlapping antisense; some tools double-count unstranded junctions |
| **Replicates** | n>=3 per condition | n=2 vs n=2 has poor calibration in most tools (especially SUPPA2) |
| **Annotation** | GENCODE basic for canonical, comprehensive for DTU/discovery | basic = high-confidence; comprehensive includes putative — affects FDR control |
| **Microexons** | PE 100+ with `--alignSJoverhangMin 8`; VAST-TOOLS | Default aligners miss 3-27nt exons |
| **Long-intron genes (TTN, brain)** | Increased `--alignIntronMax` | Default 1Mb may miss >1Mb introns |


## STAR 2-Pass Alignment


**Goal:** Maximize novel-junction sensitivity for downstream AS analysis.

**Approach:** Run STAR once per sample to discover novel junctions (pass 1), merge novel junctions across cohort, then re-align with the augmented junction set (pass 2). Cohort-style 2-pass beats per-sample basic 2-pass for differential splicing because all samples use the same junction reference.

```bash
# Pass 1: per-sample
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix pass1_${sample}_ \
    --outSJtype Standard \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1
```

```bash
# Cohort-style 2-pass: collect all SJ.out.tab from pass 1
cat pass1_*_SJ.out.tab | awk '$5 > 0 && $7 >= 3' | sort -u > cohort_novel_SJ.tab

# Pass 2: re-align with augmented junctions
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbFileChrStartEnd cohort_novel_SJ.tab \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix pass2_${sample}_ \
    --outSJtype Standard \
    --twopassMode None \
    --quantMode GeneCounts \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3
```

| Approach | Novel-junction recovery | Cohort consistency |
|----------|-------------------------|--------------------|
| 1-pass with annotation | ~80-86% (depends on GENCODE completeness) | High (annotation-based) |
| Per-sample basic 2-pass (`--twopassMode Basic`) | >=94% | Variable (each sample has its own junction set) |
| Cohort-style 2-pass (manual merge) | >=94% | High (shared junction reference) |

Per-sample 2-pass (`--twopassMode Basic`) is simpler but produces inconsistent junction sets across samples; for differential splicing the **cohort-style** version is preferred (Veeneman 2016 *Bioinformatics*).

The pass-1 filter `awk '$5 > 0 && $7 >= 3'` keeps junctions with strand info AND >=3 unique reads — adjust threshold to balance discovery vs noise.


## Junction Saturation


**Goal:** Determine whether sequencing depth is sufficient for comprehensive splicing detection.

**Approach:** Run RSeQC junction saturation; check whether the discovery curve plateaus.

```bash
junction_saturation.py \
    -i sample.bam \
    -r gencode_v45.bed \
    -o sample_junc_sat \
    -m 50 \
    -v 100000
```

```python
import subprocess
import pandas as pd

samples = ['s1.bam', 's2.bam', 's3.bam']
for sample in samples:
    subprocess.run([
        'junction_saturation.py',
        '-i', sample,
        '-r', 'gencode_v45.bed',
        '-o', sample.replace('.bam', '_junc_sat')
    ], check=True)
```

The output `*.junctionSaturation_plot.r` plots known + novel junctions vs subsampled reads.

**Plateau detection rule:** if from 80% to 100% of reads, the junction count rises by <2%, consider it plateaued. Still rising means more sequencing would yield more junctions.

For AS analysis, **plateau on the known junction curve** is the requirement; novel-junction curves often don't plateau even at deep coverage (which is biologically informative — novel junctions are inherently rarer events).


## Novel-vs-Known Junction Ratio


**Goal:** Detect annotation/mapping issues or biologically interesting cryptic splicing.

**Approach:** Classify junctions with RSeQC and compute the novel:known ratio.

```bash
junction_annotation.py -i sample.bam -r gencode_v45.bed -o sample_junc_annot
```

```python
import pandas as pd

junc = pd.read_csv('sample_junc_annot.junction.xls', sep='\t')
total = junc['total_splicing_events'].sum()

by_class = junc.groupby('annotation')['total_splicing_events'].sum()
known_frac = by_class.get('known', 0) / total
novel_frac = (by_class.get('partial_novel', 0) + by_class.get('novel', 0)) / total

print(f'known: {known_frac:.1%}, novel: {novel_frac:.1%}')
```

| Known fraction | Status | Interpretation |
|----------------|--------|----------------|
| >=80% | Healthy | Comprehensive annotation, good alignment |
| 60-80% | Acceptable | Check annotation completeness or organism |
| <60% | Suspect or interesting | Mapping artifacts, contamination, OR biologically informative |

**High novel-junction rate may be biology, not artifact:**
- **TDP-43 loss** (ALS/FTD post-mortem brain): cryptic exon de-repression in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*)
- **SF3B1-mutant** cancer (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30nt upstream of canonical (Darman 2015 *Cell Rep*)
- **Non-model organism**: GENCODE-grade annotation unavailable; novel junctions reflect annotation gaps not biology
- **Microbial / viral contamination**: reads aligning to host but with unusual junctions

If novel% >40%, drill down: check organism, check spliceosomal mutation status, check known disease signatures.


## Junction Read Overhang and Coverage


**Goal:** Profile per-junction read counts and overhang distribution to identify weakly-supported events.

**Approach:** Parse CIGAR for N (intron) operations; tally per-junction reads and minimum exon overhangs.

```python
import pysam
from collections import defaultdict

def junction_stats(bam_path):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    counts = defaultdict(int)
    min_overhang = defaultdict(lambda: float('inf'))

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary:
            continue
        ref_pos = read.reference_start
        cumulative_query = 0
        cigar = read.cigartuples
        for i, (op, length) in enumerate(cigar):
            if op == 3:
                left_match = sum(l for o, l in cigar[:i] if o in (0, 7, 8))
                right_match = sum(l for o, l in cigar[i+1:] if o in (0, 7, 8))
                overhang = min(left_match, right_match)
                key = (read.reference_name, ref_pos, ref_pos + length)
                counts[key] += 1
                min_overhang[key] = min(min_overhang[key], overhang)
            if op in (0, 2, 3, 7, 8):
                ref_pos += length

    bam.close()
    return counts, dict(min_overhang)

counts, overhang = junction_stats('sample.bam')
print(f'total junctions: {len(counts)}')
print(f'>= 10 reads: {sum(1 for c in counts.values() if c >= 10)}')
print(f'overhang >= 8 nt: {sum(1 for k, c in counts.items() if overhang[k] >= 8)}')
```

Junction reads with overhang <8 nt are common false positives, especially for novel sites. Most callers default to >=8 nt anchor for this reason. Microexon-aware aligners use overhang as low as 6 nt with explicit configuration.


## Splice Site Strength (MaxEntScan and SpliceAI)


**Goal:** Score donor and acceptor splice sites to flag weak / cryptic sites and to predict variant impact on splicing.

**Approach:** Use MaxEntScan (sequence information content) and SpliceAI (context-aware deep-learning) — they answer different questions.

```python
from maxentpy.maxent import score5, score3

donor = 'CAGGTAAGT'
acceptor = 'TTTTTTTTTTTTTTTTTTTTCAG'
print(f"5'ss MaxEnt: {score5(donor):.2f}")
print(f"3'ss MaxEnt: {score3(acceptor):.2f}")
```

| Score | Interpretation | Source |
|-------|----------------|--------|
| 5'ss MaxEnt > 8 | Strong donor | Yeo & Burge 2004 *J Comput Biol* |
| 5'ss MaxEnt 5-8 | Moderate | |
| 5'ss MaxEnt < 5 | Weak / cryptic | |
| 3'ss MaxEnt > 8 | Strong acceptor | |
| 3'ss MaxEnt < 5 | Weak / cryptic | |
| SpliceAI delta > 0.2 | PP3 supporting (ClinGen SVI 2023) | Walker 2023 *AJHG* |
| SpliceAI delta > 0.5 | PP3 moderate | |
| SpliceAI delta > 0.8 | PP3 strong | |

**MaxEntScan vs SpliceAI:**
- **MaxEntScan** scores sequence information content (intrinsic strength). Captures position-wise dependencies at the consensus.
- **SpliceAI** predicts in-vivo usage probability given full pre-mRNA context (10 kb window).
- A position with **high MaxEnt but low SpliceAI** is intrinsically strong but contextually silenced (chromatin, trans factors).
- A position with **low MaxEnt but high SpliceAI** is intrinsically weak but contextually used (enhancer-driven, e.g. weak donors stabilized by ESEs).
- Report both for variant interpretation; for variant impact see `splice-variant-prediction`.


## Picard CollectRnaSeqMetrics and Gene-Body Coverage


**Goal:** Get integrated RNA-seq QC including intronic / exonic / intergenic mapping rates and gene-body coverage uniformity.

**Approach:** Run picard CollectRnaSeqMetrics for mapping distribution; RSeQC `geneBody_coverage.py` for 5'-3' bias.

```bash
picard CollectRnaSeqMetrics \
    I=sample.bam \
    O=sample.rna_metrics.txt \
    REF_FLAT=refFlat.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
    RIBOSOMAL_INTERVALS=rRNA_intervals.interval_list

# Strandedness conversion (foot-gun):
# Reverse-stranded (Illumina TruSeq Stranded; NEB Ultra II Directional — both dUTP):
#   rMATS  --libType fr-firststrand
#   featureCounts -s 2
#   Picard STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
# Forward-stranded (Lexogen QuantSeq FWD, certain ligation-based kits):
#   rMATS  --libType fr-secondstrand
#   featureCounts -s 1
#   Picard STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND
# STAR has no library-strand flag; pass --outSAMstrandField intronMotif
# (works for any library) so downstream tools can read XS tags.

geneBody_coverage.py \
    -i sample.bam \
    -r gencode_v45.bed \
    -o sample_geneBody
```

| Metric | Healthy | Concerning |
|--------|---------|------------|
| PCT_CODING_BASES | >=50% | <30% (suggests degradation or mis-priming) |
| PCT_UTR_BASES | 20-40% | >>50% (3' bias) |
| PCT_INTRONIC_BASES | <30% (poly(A)); <60% (rRNA-depleted) | >50% (poly(A)) suggests pre-mRNA contamination |
| PCT_INTERGENIC_BASES | <10% | >20% (genomic DNA contamination) |
| MEDIAN_5PRIME_TO_3PRIME_BIAS | 0.7-1.3 | >2 or <0.5 (severe degradation) |
| Gene body coverage curve | Flat | Strong 3' skew = RIN low or library mis-prep |

3' bias (degraded RNA) directly reduces splicing-event detection because junction reads scatter across the gene body; with 3' bias they concentrate near the 3' end and miss CDS junctions.


## Strandedness Verification


```bash
infer_experiment.py -i sample.bam -r gencode_v45.bed -s 200000
```

Output reports the fraction of reads consistent with each library type:

| Output pattern | Library type | rMATS `--libType` |
|----------------|---------------|---------------------|
| ~50% / ~50% | Unstranded | `fr-unstranded` |
| >=90% "++ , --" | Forward-stranded | `fr-secondstrand` |
| >=90% "+- , -+" | Reverse-stranded (Illumina TruSeq stranded) | `fr-firststrand` |

**Wrong strand setting halves usable junction reads** — always verify before quantification. RSeQC `infer_experiment.py` is fast and authoritative.


## Annotation Choice


| GENCODE level | Contents | Use for |
|---------------|----------|---------|
| Basic | High-confidence canonical isoforms | Standard rMATS, leafcutter, SUPPA2 |
| Comprehensive | All transcripts including putative/predicted | DTU pipelines (DRIMSeq+DEXSeq, satuRn), isoform discovery |
| RefSeq | NCBI curated | Less complete than GENCODE; legacy use |
| Ensembl | Same content as GENCODE in vertebrates | Different attribute conventions |

Comprehensive captures more biology but inflates DTU multiple-testing burden and includes annotation noise. For event-level (rMATS) AS, basic is usually adequate; for transcript-level DTU (DRIMSeq, satuRn), comprehensive may be necessary to capture rare isoforms.


## rRNA Contamination Check


```bash
fastq_screen --conf fastq_screen.conf --threads 8 sample_R1.fq.gz
```

Or post-alignment:

```bash
samtools view -c sample.bam | awk '{print "total:",$0}'
samtools view -c -L rRNA_intervals.bed sample.bam | awk '{print "rRNA:",$0}'
```

| rRNA fraction | Library type | Status |
|----------------|---------------|--------|
| >=20% | "depleted" | Failed depletion; redo |
| 5-20% | "depleted" | Acceptable; some rRNA leakage |
| <5% | poly(A) | Healthy |
| <5% | "depleted" | Excellent depletion |
| 1-3% | poly(A) | Suggests RNA degradation |

>5% rRNA in a poly(A) library suggests degraded RNA; >20% in a "depleted" library indicates failed depletion.


## Troubleshooting Low Event Detection


| Issue | Possible causes | Solutions |
|-------|-----------------|-----------|
| Few events called | Low depth; short reads; SE; wrong strand | Increase depth; use PE150; verify libType |
| High novel junctions | Annotation gaps; mapping artifacts; biology (TDP-43, SF3B1) | Update annotation; check 2-pass; consider biology |
| Low IR detection | poly(A) library | Use rRNA depletion |
| Microexons missing | Default aligner anchors too long | VAST-TOOLS, MicroExonator, or long-read |
| Many weak splice sites | Cryptic splicing | Validate with MaxEnt + SpliceAI; consider RNA-seq from secondary tissue |
| FDR uncalibrated at low n | n=2 vs n=2 | Use leafcutter or Shiba; avoid SUPPA2 alone |
| PSI variance high across replicates | Library prep / RIN inconsistency | Check RIN; consider RNA degradation |
| Sashimi plot mismatch with PSI | Junction-imbalance bias in rMATS | Run Shiba; or filter by overhang distribution |


## Statistical Model Taxonomy


| Tool | Model | Test statistic | Min reps per group | Calibration regime | Fails when |
|------|-------|-----------------|---------------------|---------------------|------------|
| rMATS-turbo | Binomial counts with hierarchical PSI variance | LRT on \|ΔPSI\| > `cutoff` (default 0.0001) | n>=3 | Well-calibrated at n>=3 with adequate junction reads | Junction read imbalance; very low coverage; uncorrected for confounders |
| leafcutter | Dirichlet-multinomial GLM at cluster level | LRT on group factor | n>=2 (n>=3 preferred) | Strong at n>=3; novel-junction-friendly | Undersampled clusters (DM dispersion unstable); cluster topology arbitrariness |
| MAJIQ deltapsi | Beta-binomial bootstrap → posterior over PSI per LSV | P(\|ΔPSI\| > T) threshold (T=0.2) | n>=3 | Replicate-structured n=3 vs n=3 | Cohorts where between-sample variability dominates between-group |
| MAJIQ HET | Same model, heterogeneity-aware | Per-LSV permutation-based test | n>=10 | n>=10 vs n>=10 cohort designs | Tightly-controlled small replicate experiments |
| SUPPA2 (empirical) | Empirical null from between-replicate ΔPSI | ECDF on \|ΔPSI\| conditioned on TPM | n>=4 | n>=4 vs n>=4 with paired-end deep sequencing | n<=3 vs n<=3 (sparse null collapses) |
| SUPPA2 (classical) | Wilcoxon rank-sum on PSI distributions | Wilcoxon p-value | n>=2 | Small samples; non-parametric backup | Cassette events with tight PSI distributions |
| Shiba (2025) | Beta-binomial with explicit junction-imbalance correction | LRT | n>=2 | n=2-3 vs n=2-3 | Established benchmarks limited (new tool) |
| LeafcutterMD | Dirichlet-multinomial outlier mode | Per-sample p-value | n=1 vs cohort >=20 | Single-patient vs cohort | Too few controls (<20) |
| FRASER 2.0 | Beta-binomial autoencoder on Intron Jaccard Index | Per-sample p-value with delta cutoff | n=1 vs cohort >=20 | n>=20 control cohort, single-patient query | See `outlier-splicing-detection` for this regime |

The first decision is which **regime** the design falls into: between-group with replicates, heterogeneous cohort, or single-sample-vs-cohort. Within each regime, tool choice is much smaller (1-2 options).

Comprehensive 2023-2026 benchmarks: Olofsson 2023 *Brief Bioinform*; Tran 2025 *WIREs RNA*; Kubota 2025 *NAR*. Methodology evolves — verify benchmarks and tool docs before reporting. Default 2026 recommendation: run **two complementary tools** (rMATS + leafcutter) and require concordance for high-confidence calls.


## Decision Tree by Experimental Design


| Scenario | Recommended tool | Why | Threshold |
|----------|------------------|-----|-----------|
| Standard n=3 vs n=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Two algorithmic families; concordant hits = high-confidence | FDR<0.05, \|ΔPSI\|>0.10 |
| n=2 vs n=2 small pilot | Shiba | Junction-imbalance correction matters most at low coverage | FDR<0.10, \|ΔPSI\|>0.10 |
| n=10+ vs n=10+ heterogeneous (clinical, GTEx-style) | MAJIQ V3 HET | HET designed for between-sample heterogeneity | P(\|ΔPSI\|>0.2)>0.95 |
| Single rare-disease patient vs panel of n>=20 | FRASER 2.0 (see outlier-splicing-detection) | Outlier detection statistical model is fundamentally different | padj<0.05, \|delta-jaccard\|>=0.1 |
| Time-course / multi-condition design | Custom DEXSeq or limma on PSI matrix | rMATS/leafcutter primarily 2-group | FDR<0.05 on time:group interaction |
| Paired tumor-normal | rMATS with `--paired-stats` | Paired test reduces inter-patient variance | FDR<0.05, paired \|ΔPSI\|>0.10 |
| Cancer with spliceosomal mutation (SF3B1, U2AF1) | leafcutter or MAJIQ denovo | Cryptic events not in annotation | FDR<0.05; check 3'ss shifts in IGV |
| TDP-43 loss / ALS post-mortem | leafcutter denovo | Cryptic exons not in annotation | FDR<0.05; expect UNC13A, STMN2 |
| Non-model organism without GENCODE-grade annotation | leafcutter | Annotation-free | FDR<0.05, \|ΔPSI\|>0.10 |
| Long-read available | rMATS-long, FLAIR diffSplice | See long-read-splicing | Tool-specific |


## rMATS-turbo Differential Analysis


**Goal:** Detect statistically significant differential splicing between two groups from BAMs.

**Approach:** Run rMATS-turbo without `--statoff`, then filter by FDR + ΔPSI + per-replicate coverage.

```bash
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --cstat 0.05
```

`--cstat 0.05` tests `|ΔPSI| > 0.05`; raise to 0.10 for stricter discovery. `--novelSS` enables novel-junction discovery (recommended with STAR 2-pass). For paired designs, add `--paired-stats`.

```python
import pandas as pd
import numpy as np

se = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

def min_per_rep(s):
    return s.str.split(',').apply(lambda x: min(int(v) for v in x))

se['min_inc'] = min_per_rep(se['IJC_SAMPLE_1']).combine(min_per_rep(se['IJC_SAMPLE_2']), min)
se['min_skip'] = min_per_rep(se['SJC_SAMPLE_1']).combine(min_per_rep(se['SJC_SAMPLE_2']), min)

significant = se[
    (se['FDR'] < 0.05) &
    (se['IncLevelDifference'].abs() > 0.10) &
    ((se['min_inc'] + se['min_skip']) >= 10)
].copy()

significant['score'] = -np.log10(significant['FDR']) * significant['IncLevelDifference'].abs()
top = significant.nlargest(50, 'score')
```


## leafcutter Differential Intron Usage


**Goal:** Detect differential intron-cluster usage annotation-free, capturing novel junctions and complex multi-junction events.

**Approach:** Extract junctions with regtools, cluster introns by shared splice sites, run cluster-level Dirichlet-multinomial test.

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

```r
library(leafcutter)

groups <- data.frame(
    sample = c('s1', 's2', 's3', 's4', 's5', 's6'),
    group = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
)
write.table(groups, 'groups.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

system('leafcutter_ds.R --num_threads 4 --exon_file gencode_exons.txt.gz \
    leafcutter_perind_numers.counts.gz groups.txt -o ds_results')

cluster_sig <- read.table('ds_results_cluster_significance.txt', header = TRUE, sep = '\t')
intron_effects <- read.table('ds_results_effect_sizes.txt', header = TRUE, sep = '\t')

sig_clusters <- subset(cluster_sig, p.adjust < 0.05)
```

**LeafCutter2** (Quan 2025 *bioRxiv*) extends leafcutter with NMD-aware classification of unproductive splicing — useful when AS-NMD coupling is the question.


## MAJIQ V3 Differential Analysis


**Goal:** Detect differential LSVs with full posterior distributions over ΔPSI; ideal for complex multi-junction events and heterogeneous cohorts.

**Approach:** Build splice graph → compute coverage per group → run deltapsi (replicate-structured) or heterogen (cohort-style).

```bash
majiq build annotation.gff3 -c settings.ini -j 8 -o build_output

majiq deltapsi \
    -grp1 build_output/ctrl1.majiq build_output/ctrl2.majiq build_output/ctrl3.majiq \
    -grp2 build_output/trt1.majiq build_output/trt2.majiq build_output/trt3.majiq \
    -n control treatment \
    -o deltapsi_output \
    --minreads 10 --minpos 3 \
    -j 8

majiq heterogen \
    -grp1 build_output/het_ctrl{1..20}.majiq \
    -grp2 build_output/het_trt{1..20}.majiq \
    -n control treatment \
    -o heterogen_output \
    -j 8

voila view -p 5000 -j 8 build_output/splicegraph.zarr deltapsi_output/control_treatment.deltapsi.voila -o voila_html
```

MAJIQ V3 (Slaff et al *bioRxiv* 2024; public release 2025) uses Zarr storage (`splicegraph.zarr`); V2's SQLite splicegraph is deprecated. MAJIQ reports posterior probability `P(|ΔPSI| > 0.2)`; thresholds are interpreted differently from FDR. Use HET for n>=10 vs n>=10 cohort designs (clinical, GTEx-style); deltapsi for tightly controlled n=3 vs n=3.


## SUPPA2 Differential Analysis


**Goal:** Quick differential splicing from existing transcript quantifications, useful as a sanity check or pilot.

**Approach:** Generate per-condition PSI files from Salmon TPM, then run `diffSplice` with empirical or classical p-values.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e ctrl_tpm.tsv -o ctrl_${ev}
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e trt_tpm.tsv -o trt_${ev}

    suppa.py diffSplice \
        -m empirical \
        -gc \
        -i events_${ev}_strict.ioe \
        -p ctrl_${ev}.psi trt_${ev}.psi \
        -e ctrl_tpm.tsv trt_tpm.tsv \
        -o diff_${ev}
done
```

For n<=3 designs, switch `-m classical` (Wilcoxon). Empirical null requires sufficient between-replicate observations to construct.


## Shiba for Low-Coverage / Few-Replicate Designs


**Goal:** Detect differential splicing with explicit junction-imbalance correction — addresses a known false-positive source for rMATS-style methods.

**Approach:** Shiba is a Snakemake-based pipeline configured via YAML. Install via bioconda, write a config file describing groups + BAMs, then run with snakemake.

```bash
conda install -c bioconda shiba

# Edit config.yaml with reference GTF, BAM groups, output dir, thresholds
# Then run the Snakemake workflow:
snakemake -s snakeshiba.smk \
    --configfile config.yaml \
    --cores 8 \
    --use-singularity \
    --singularity-args "--bind $HOME:$HOME"
```

Shiba (Kubota 2025 *NAR*) reportedly outperforms rMATS at n=2 vs n=2 by correcting differential mappability between inclusion and skipping junctions; community calibration still emerging. See https://sika-zheng-lab.github.io/Shiba/ for the full config.yaml schema.


## Reconciliation: When Tools Disagree


The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction imbalance OR rMATS event hits annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba on the same locus |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Verify `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | Event class mismatch (rMATS calls SE positive, leafcutter sees A5SS shift in same cluster) | Manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |
| All tools null but biology suggests change | Underpowered design or wrong regime | Increase replicates; check whether outlier-splicing-detection regime applies |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based). Document both calls and any explainable disagreements.


## rMATS Output Columns Reference


| Column | Meaning |
|--------|---------|
| IJC_SAMPLE_1 / SJC_SAMPLE_1 | Comma-delimited inclusion / skipping junction counts per replicate, group 1 |
| IJC_SAMPLE_2 / SJC_SAMPLE_2 | Same for group 2 |
| IncFormLen / SkipFormLen | Effective lengths normalizing PSI for differential mapping opportunity |
| upstreamES/EE, downstreamES/EE | Flanking exon coordinates (genomic order; strand-agnostic in column meaning) |
| exonStart_0base / exonEnd | Cassette exon coordinates (0-based half-open) |
| PValue | LRT p-value of \|ΔPSI\| > cutoff |
| FDR | BH-adjusted PValue within event class |
| IncLevel1, IncLevel2 | Comma-delimited per-replicate PSI values |
| IncLevelDifference | mean(IncLevel1) - mean(IncLevel2); sign matches --b1 - --b2 order |


## Replicate Count and Power


| Design | Recommended tools | Expected power for ΔPSI=0.2 |
|--------|-------------------|------------------------------|
| n=2 vs n=2 | leafcutter or Shiba; **avoid SUPPA2** | Marginal; many real effects missed |
| n=3 vs n=3 | rMATS-turbo + leafcutter | Adequate at moderate coverage; standard |
| n=5 vs n=5 | rMATS or leafcutter, MAJIQ deltapsi | Good; recommended for publication |
| n=10+ vs n=10+ heterogeneous | MAJIQ-HET | Designed for this scale |
| Single patient vs n=20+ controls | leafcutterMD or FRASER2 | Outlier regime; see outlier-splicing-detection |

For an effect-size of |ΔPSI|=0.10 (typical biological signal), power generally requires n>=4 and >=20 junction reads per replicate. Below this, expect to miss most real changes.


## Significance and Effect-Size Thresholds


| Stringency | \|ΔPSI\| | FDR | Use case |
|------------|----------|-----|----------|
| Lenient | > 0.05 | < 0.10 | Discovery, exploratory, hypothesis generation |
| Standard | > 0.10 | < 0.05 | Publication; default reporting threshold |
| Stringent | > 0.20 | < 0.01 | Validation cohort, follow-up targets |

For MAJIQ: posterior probability `P(|ΔPSI| > 0.2) >= 0.95` is roughly equivalent to standard stringency. Always document tool, threshold, and rationale.

**Biologically meaningful ΔPSI varies by context:**
- A poison exon shift of |ΔPSI|=0.10 can halve functional protein (huge biology, modest number).
- A stoichiometric isoform shift of |ΔPSI|=0.10 may be physiologically silent.
- Therapeutic ASO target: SMA nusinersen aims for ΔPSI~+0.30 in SMN2 exon 7.


## Confounder Handling


**rMATS** does not natively accept arbitrary covariates. Workarounds:
1. **Stratification**: run rMATS within each batch separately and meta-analyze.
2. **PSI residuals (logit-transformed)**: PSI is bounded [0,1]; raw linear regression near the boundaries is biased. Logit-transform first, regress on confounders, then test residuals.
3. **Switch to leafcutter** (R function accepts `confounders` matrix; CLI accepts confounders as additional columns in the groups file).

```python
import numpy as np
import statsmodels.formula.api as smf

# logit-transform PSI before residualization (PSI is bounded [0,1])
eps = 1e-3
psi['logit_psi'] = np.log((psi['psi'].clip(eps, 1 - eps)) / (1 - psi['psi'].clip(eps, 1 - eps)))
psi['psi_resid'] = smf.ols('logit_psi ~ batch + RIN', data=psi).fit().resid
# then test psi_resid by group via Wilcoxon
```

**leafcutter** accepts confounders two ways:
- **R function**: `differential_splicing(counts, x, confounders=numeric_matrix)` accepts a numeric covariate matrix
- **CLI script**: `leafcutter_ds.R` reads confounders from **additional columns in the groups file** (3rd, 4th, ... columns), NOT from a `--confounders` flag

**MAJIQ** does not accept arbitrary confounders; use stratification or switch tool.

**Always check confounding before reporting:** PCA on PSI matrix; if PC1 separates by batch rather than group, the comparison is confounded.


## Multi-Group / Multi-Factor Designs


| Design | Approach |
|--------|----------|
| 3 groups (e.g. drug A, drug B, control) | Pairwise rMATS or leafcutter; OR limma/DESeq2 on logit-PSI matrix |
| Time-course (e.g. 0h, 6h, 24h) | DEXSeq on event counts with time as factor; or limma::lmFit on PSI matrix |
| 2x2 factorial (genotype × treatment) | DEXSeq with interaction term; rMATS pairwise on interaction subsets |
| Continuous covariate (dose, age) | limma::lmFit on logit-PSI ~ covariate |

For complex designs, custom regression on the PSI matrix is more flexible than rMATS/leafcutter pairwise.


## Result Prioritization


**Goal:** Rank events by combined statistical and biological significance for follow-up.

**Approach:** Composite score combining FDR and effect size, then enrich for biology (RBP binding, NMD sensitivity, conservation, disease relevance).

```python
import pandas as pd
import numpy as np

sig['score'] = -np.log10(sig['FDR']) * sig['IncLevelDifference'].abs()
sig['exon_length'] = sig['exonEnd'] - sig['exonStart_0base']
sig['nmd_likely'] = (sig['exon_length'] % 3 != 0)
top_events = sig.nlargest(50, 'score')
```

Cross-reference top hits with:
- **eCLIP/ENCODE RBP target databases** (POSTAR3, oRNAment, RBP2GO) → candidate trans-regulators
- **Disease-specific signatures**: SF3B1 cryptic 3'ss for MDS/CLL/UM; TDP-43 cryptic exons (UNC13A, STMN2) for ALS/FTD
- **Conservation**: VastDB cross-species PSI for evolutionary support
- **Splice-site predictions**: SpliceAI scores for the involved sites (see splice-variant-prediction)


## Decision Tree by Goal


| Goal | Recommended tool |
|------|-------------------|
| Validate a specific rMATS hit | rmats2sashimiplot (one-line) or ggsashimi (custom) |
| Validate a leafcutter cluster | leafviz (interactive) or ggsashimi with cluster coordinates |
| Validate a MAJIQ LSV (complex topology) | MAJIQ-VOILA (only tool that shows full LSV graph) |
| Publication-quality two-condition comparison | ggsashimi `-O 3 -A mean_j` for grouped overlay |
| Multi-track figure (RNA-seq + H3K4me3 + ATAC) | pyGenomeTracks |
| Quick ad-hoc browsing during development | IGV sashimi |
| Tool-agnostic batch heatmap of significant events | Jutils |
| Interactive cohort-level filtering of leafcutter results | leafviz Shiny |


## ggsashimi for Publication Overlays


**Goal:** Generate publication-quality sashimi plot for a region with samples grouped by condition and per-sample tracks aggregated.

**Approach:** Define samples + groups + colors in a TSV (no header), then call ggsashimi with coordinates, GTF, and visual flags.

```python
import subprocess
import pandas as pd

groups = pd.DataFrame({
    'bam': ['ctrl1.bam', 'ctrl2.bam', 'ctrl3.bam', 'trt1.bam', 'trt2.bam', 'trt3.bam'],
    'group': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment'],
    'color': ['#1f77b4'] * 3 + ['#ff7f0e'] * 3
})
groups.to_csv('sashimi_groups.tsv', sep='\t', index=False, header=False)

subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr17:43094000-43125000',
    '-o', 'BRCA1_sashimi',
    '-M', '10',
    '--alpha', '0.25',
    '--height', '3',
    '--width', '10',
    '--shrink',
    '--fix-y-scale',
    '--ann-height', '4',
    '-g', 'gencode_v45.gtf',
    '--base-size', '14',
    '-O', '3',
    '-A', 'mean_j',
    '-F', 'pdf'
], check=True)
```

Key ggsashimi flags (Garrido-Martin 2018 *PLoS Comput Biol*):
- `--overlay 3` (or `-O 3`): aggregate multiple samples within a group into a single overlay track with summary statistics — its signature feature
- `-A mean_j`: junction aggregation method (`mean`, `median`, `mean_j` accounts for sample-wise normalization); use `mean_j` for biological replicates
- `--shrink`: rescale long introns (>2x flanking exons) for compact display
- `--fix-y-scale`: identical y-axis across groups (essential for visual comparison)
- `--alpha 0.25`: transparency for per-sample coverage in overlay mode
- `-M 10`: minimum junction reads to display (lower = noisier; 5-10 typical; raise to 20+ for crowded plots)
- `--ann-height`: gene annotation track height
- `-F pdf`: output format (pdf, png, svg, eps)


## Batch Plotting from rMATS Hits


**Goal:** Auto-generate sashimi plots for all significant rMATS differential events.

**Approach:** Parse SE.MATS.JC.txt, expand coordinates to flanking exons + 500nt context, iterate ggsashimi.

```python
import subprocess
import pandas as pd
from pathlib import Path

diff = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')
sig = diff[(diff['FDR'] < 0.05) & (diff['IncLevelDifference'].abs() > 0.10)]

Path('sashimi_plots').mkdir(exist_ok=True)
for idx, ev in sig.head(25).iterrows():
    region = f'{ev["chr"]}:{ev["upstreamES"] - 500}-{ev["downstreamEE"] + 500}'
    safe_name = f'{ev["geneSymbol"]}_{ev["chr"]}_{ev["upstreamES"]}'
    subprocess.run([
        'ggsashimi.py',
        '-b', 'sashimi_groups.tsv',
        '-c', region,
        '-o', f'sashimi_plots/{safe_name}',
        '-M', '5',
        '--shrink',
        '--fix-y-scale',
        '-O', '3',
        '-A', 'mean_j',
        '-g', 'annotation.gtf',
        '-F', 'pdf'
    ], check=True)
```

For MXE events, plot from upstreamES of exon 1 to downstreamEE of exon 2 to show both alternative exons in the same figure.


## rmats2sashimiplot


**Goal:** Plot directly from rMATS event coordinates without manual region calculation.

**Approach:** Pass rMATS event file + BAM lists + event type; rmats2sashimiplot extracts coordinates and produces per-event PDFs.

```bash
rmats2sashimiplot \
    --b1 ctrl1.bam,ctrl2.bam,ctrl3.bam \
    --b2 trt1.bam,trt2.bam,trt3.bam \
    -t SE \
    -e rmats_output/SE.MATS.JC.txt \
    --l1 Control \
    --l2 Treatment \
    -o sashimi_rmats \
    --exon_s 1 \
    --intron_s 5 \
    --color '#1f77b4,#ff7f0e' \
    --group-info group_def.txt
```

`--exon_s 1 --intron_s 5` shrinks intron-to-exon visual ratio 5:1 (introns drawn 1/5 their actual length). The `--group-info` flag (newer versions) allows custom replicate groupings.


## MAJIQ-VOILA Interactive HTML


**Goal:** Browse LSV posterior PSI distributions interactively with splice-graph topology.

**Approach:** Run `voila` on MAJIQ output to generate self-contained HTML.

```bash
# MAJIQ V3 (June 2025+) uses Zarr-format splicegraph (V2's .sql is deprecated)
voila view -p 5000 -j 8 build/splicegraph.zarr psi_output/sample.psi.voila -o voila_psi_html

voila view -p 5000 -j 8 build/splicegraph.zarr deltapsi_output/group1_group2.deltapsi.voila -o voila_dpsi_html
```

VOILA shows:
- Complete LSV graphs (single source / single target nodes)
- Per-junction posterior PSI violin plots
- ΔPSI distributions across all conditions
- Confidence by junction within an LSV

**The only tool that visualizes complex multi-junction LSVs intuitively.** For events that don't fit canonical SE/A5SS/A3SS, VOILA is the visualization of choice.


## leafviz Shiny App


**Goal:** Browse leafcutter clusters with intron-level effects, sashimi-like plots, and NMD annotation.

**Approach:** Prepare leafviz input from leafcutter differential output, then launch Shiny.

```bash
prepare_results.R \
    -o leafviz \
    -m groups.txt \
    leafcutter_perind_numers.counts.gz \
    ds_results_cluster_significance.txt \
    ds_results_effect_sizes.txt \
    annotation_codes
```

```r
library(leafviz)
run_leafviz('leafviz.RData')
```

Standalone alternative: `jackhump/leafviz` GitHub repo for the lightweight installable subset. Useful for cohort-level interactive filtering.


## Jutils for Tool-Agnostic Output


**Goal:** Visualize differential splicing output uniformly across rMATS, leafcutter, SUPPA2, and MAJIQ.

**Approach:** Convert tool output to Jutils' standard format, then plot.

```bash
jutils convert -t rmats -i SE.MATS.JC.txt -o rmats_jutils.tsv
jutils heatmap -i rmats_jutils.tsv -o heatmap.pdf --top 50
jutils sashimi -i rmats_jutils.tsv -b sashimi_groups.tsv -g annotation.gtf -o sashimi_jutils/
jutils venn -i rmats_jutils.tsv leafcutter_jutils.tsv -o overlap_venn.pdf
```

(Yang 2021 *Bioinformatics*) Useful when comparing multiple tools' outputs across publications or doing meta-analysis.


## pyGenomeTracks for Multi-Track Figures


**Goal:** Combine splicing with chromatin or coverage tracks for publication figures.

**Approach:** Define tracks in an INI file (genes, BAM, BigWig, BED), then run `pyGenomeTracks --tracks tracks.ini --region ... -o figure.pdf`.

```ini
[gene_models]
file = annotation.gtf
height = 3
title = GENCODE v45
fontsize = 10
file_type = gtf

[ctrl_coverage]
file = ctrl_merged.bw
title = Control
color = #1f77b4
height = 3
file_type = bigwig

[trt_coverage]
file = trt_merged.bw
title = Treatment
color = #ff7f0e
height = 3
file_type = bigwig

[junctions]
file = junctions.bedpe
title = Junctions
height = 2
file_type = links
links_type = arcs
```

The `junctions.bedpe` file must be in **BEDPE format** (6 columns: chr1 start1 end1 chr2 start2 end2 [+ optional score]). Convert from regtools .bed12 junctions:

```bash
# Convert regtools junctions BED12 to BEDPE for pyGenomeTracks.
# regtools BED12 column 11 is blockSizes (anchor_left, anchor_right);
# column 12 is blockStarts (0, intron_length + anchor_left).
# Intron start = chromStart + anchor_left = $2 + a[1]
# Intron end   = chromStart + blockStarts[2] = $2 + b[2]
awk 'BEGIN{OFS="\t"} {split($11,a,","); split($12,b,","); s=$2+a[1]; e=$2+b[2]; print $1, s, s+1, $1, e-1, e, $5}' \
    regtools_junctions.bed > junctions.bedpe
```

```bash
pyGenomeTracks --tracks tracks.ini --region chr17:43094000-43125000 -o figure.pdf
```


## Reading Sashimi Plots (Interpretation Guide)


| Visual element | What it represents |
|----------------|--------------------|
| Filled coverage track | Read coverage at each genomic position (depth-normalized in `-A` mode) |
| Arc / curve between exons | Junction-spanning reads; arc connects donor to acceptor |
| Number on arc | Count of junction-spanning reads (raw, not normalized, unless `-A` set) |
| Arc thickness | Often proportional to read count (tool-dependent) |
| Gene model below | Exons (boxes) and introns (lines) from GTF |
| Multiple parallel tracks | Per-sample (default) or per-group (with `-O`) |

**Junction count interpretation:** the number on an arc is the absolute count of reads whose CIGAR string contained an `N` operation matching that intron coordinate. Higher = more usage. Compare counts on inclusion vs skipping arcs to estimate PSI visually.

**Color convention:** by convention, control = blue (`#1f77b4`), treatment = orange (`#ff7f0e`); always document. Use ColorBrewer or matplotlib defaults for >2 groups.


## Customization Reference


| Visual goal | ggsashimi flag |
|-------------|-----------------|
| Reduce intron whitespace | `--shrink` |
| Identical y-axis across groups | `--fix-y-scale` |
| Per-group overlay aggregation | `-O 3 -A mean_j` |
| Larger figure | `--width 12 --height 4` |
| Bigger fonts | `--base-size 16` |
| Vector output | `-F pdf` or `-F svg` |
| Custom palette | Edit colors in groups TSV |
| Filter junction noise | `-M 10` (raise to 20+) |
| Transparency | `--alpha 0.25` |
| GTF feature filter | `--gtf-filter protein_coding` |


## Best Practices


| Tip | Rationale |
|-----|-----------|
| Use `--shrink` for genes with large introns | Keeps exons visible (TTN, brain genes with multi-kb introns) |
| `--fix-y-scale` for cross-group comparisons | Otherwise auto-rescaling visually exaggerates differences |
| Aggregate replicates with `-O 3 -A mean_j` | Reduces clutter; per-sample variance still shown via alpha |
| Limit to 3-4 groups per figure | More becomes hard to read |
| Include 200-500 nt flanking exons | Show full splicing context |
| For MXE events, plot both alternative exons | Otherwise only half of the event is visible |
| Check accessibility colors | Use ColorBrewer-safe palettes for color-blind readers |
| Always include a legend | Sashimi figures without legends are uninformative for non-experts |
| Specify output format explicitly | PDF for publication; PNG for slides; SVG for editing |


## Troubleshooting


| Issue | Cause | Solution |
|-------|-------|----------|
| No junctions shown | Default `-M 10` too strict | Lower to `-M 3` or `-M 5` |
| Plot too crowded | Many samples without aggregation | Use `-O 3` to overlay groups |
| Annotation missing or wrong gene | GTF lacks gene_name attribute or wrong build | Verify GTF version vs BAM reference; switch to `--gtf-filter protein_coding` |
| Memory issues on large regions | >100 kb regions with many samples | Plot smaller windows or pre-extract reads with samtools view |
| Y-axis dominated by one peak | Outlier sample | Use `-A mean_j` to aggregate; or filter outlier |


## Predictor Taxonomy


| Family | Architecture | Output | Fails when |
|--------|--------------|--------|------------|
| Context-aware CNN | 10 kb dilated ResNet | Per-position donor/acceptor probability | Long-range (>5 kb) regulatory effects; tissue-specific events |
| Tissue-aware CNN/transformer | Same arch + multi-tissue training | Per-tissue ΔPSI | Tissue not in training set; novel cell types |
| Modular per-region CNN | Separate sub-models for 5'ss/3'ss/exon/intron | Calibrated quantitative ΔPSI | Atypical events; complex multi-junction effects |
| Foundation transformer | Pretrained on broad genomic context | Splice probability or ΔPSI | New tools; less battle-tested |
| Empirical lookup | Public RNA-seq event database | Top-N most likely mis-splicing outcomes | Variant types not represented in training cohorts |
| Composite score | Blend of multiple predictors | Single scaled score | When component predictors disagree internally |


## Decision Tree by Use Case


| Use case | Recommended approach |
|----------|----------------------|
| Clinical variant report (single variant, ACMG classification) | SpliceAI default 50nt + ClinGen SVI 2023 thresholds |
| Tissue-specific clinical question (brain disease, cardiomyopathy) | SpliceAI + Pangolin (tissue-matched) |
| Unsolved Mendelian case (suspect deep-intronic) | SpliceAI extended window (-D 500-2000) + SpliceVault |
| VUS panel screening | SpliceAI + Pangolin + MMSplice concordance scoring |
| Predict consequence of canonical-disrupting variant | SpliceVault top-N empirical events |
| Branchpoint variant suspected | BPHunter (branchpoint screen) — SpliceAI is weak here |
| Splice-switching ASO design (target ESE/ESS occlusion) | SpliceAI on masked sequence + RNAfold accessibility |
| Validate predicted splice change in patient | RNA-seq + FRASER2 (see outlier-splicing-detection) |
| Pseudoexon prediction in deep intron | SpliceAI extended window + CI-SpliceAI; require RNA validation |


## ClinGen SVI 2023 Framework


The ClinGen Sequence Variant Interpretation (SVI) splicing subgroup (Walker 2023 *Am J Hum Genet*; Riepe 2024 *Genet Med*) extended the ACMG/AMP 2015 framework with explicit splice-prediction rules.

| Evidence code | Threshold | Notes |
|----------------|-----------|-------|
| **PP3** (supporting pathogenic) | SpliceAI delta >= 0.20 | Computational evidence supporting pathogenicity |
| **PP3 moderate** | SpliceAI delta >= 0.50 | Or concordance across multiple predictors |
| **PP3 strong** | SpliceAI delta >= 0.80 | Typically requires concordance + canonical site |
| **BP4** (supporting benign) | SpliceAI delta <= 0.10 | Computational evidence against pathogenicity |
| **PVS1** (very strong null) | Canonical +/-1, +/-2 site disruption with predicted LoF + NMD | Requires gene where LoF is established mechanism (Abou Tayoun 2018 *Hum Mutat* PVS1 decision tree) |
| **PS3 / BS3** (functional) | RNA evidence (RT-PCR, RNA-seq, minigene) | Supersedes computational evidence |

**Operational rules:** Computational evidence (PP3/BP4) is *supporting*, not standalone. Splicing variants benefit from concordance across SpliceAI + Pangolin + MMSplice. RNA validation supersedes prediction. Always log SpliceAI version, distance window, and reference transcript. SpliceAI alone is **not sufficient** for PVS1; canonical site disruption requires gene-level LoF context.


## SpliceAI Workflow


**Goal:** Annotate VCF variants with per-variant delta scores for splice-site change.

**Approach:** Run `spliceai` CLI with reference genome and annotation; parse INFO field for delta scores. **SpliceAI is human-only** (`-A grch37` or `-A grch38`); the model was trained on GENCODE human and does not directly transfer to mouse, fly, or other species. For mouse, retrained variants exist (e.g. mouseSpliceAI); for other species, use Pangolin (4 species: human, mouse, rat, rhesus macaque) or accept that prediction will be unreliable.

```bash
spliceai \
    -I input.vcf \
    -O output.vcf \
    -R GRCh38.primary_assembly.genome.fa \
    -A grch38 \
    -D 50 \
    -M 0
```

`-D 50` = distance window in nt around variant (default 50). For deep-intronic variants suspected of creating pseudoexons, raise to **500-2000**:

```bash
spliceai -I input.vcf -O output_extended.vcf -R genome.fa -A grch38 -D 500 -M 1
```

`-M 0` (default) returns raw scores; `-M 1` masks splice gains at annotated sites and losses at unannotated sites (cleaner for clinical use). Output INFO format: `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. Delta score = max(DS_AG, DS_AL, DS_DG, DS_DL).

```python
import pandas as pd
import re

def parse_spliceai_vcf(vcf_path):
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info = fields[7]
            m = re.search(r'SpliceAI=([^;]+)', info)
            if not m:
                continue
            for ann in m.group(1).split(','):
                parts = ann.split('|')
                allele, symbol = parts[0], parts[1]
                ds = [float(p) if p != '.' else 0 for p in parts[2:6]]
                dp = parts[6:10]
                rows.append({
                    'chrom': fields[0], 'pos': int(fields[1]),
                    'ref': fields[3], 'alt': allele,
                    'gene': symbol,
                    'DS_AG': ds[0], 'DS_AL': ds[1],
                    'DS_DG': ds[2], 'DS_DL': ds[3],
                    'delta_max': max(ds),
                })
    return pd.DataFrame(rows)

df = parse_spliceai_vcf('output.vcf')
df['acmg_evidence'] = pd.cut(
    df['delta_max'],
    bins=[-0.01, 0.10, 0.20, 0.50, 0.80, 1.01],
    labels=['BP4', 'inconclusive', 'PP3_supporting', 'PP3_moderate', 'PP3_strong']
)
```

DS labels: AG = acceptor gain, AL = acceptor loss, DG = donor gain, DL = donor loss.


## Pangolin for Tissue-Specific Prediction


**Goal:** Get tissue-specific splice impact predictions when disease tissue is known.

**Approach:** Run Pangolin CLI with VCF + reference + gffutils annotation database.

```bash
python -c "import gffutils; gffutils.create_db('gencode.v45.annotation.gff3', 'gencode.db', force=True)"

pangolin \
    input.vcf \
    GRCh38.primary_assembly.genome.fa \
    gencode.db \
    pangolin_output \
    -d 500 \
    -m True \
    -s 0.2
```

`-m True` masks splice gains at annotated sites and losses at unannotated sites (recommended for clinical use). `-s 0.2` outputs all sites with predicted change >= cutoff.

Pangolin output is a VCF with per-tissue predictions across the **4 tissues used at training: brain, heart, liver, testis** (Zeng & Li 2022 *Genome Biol*). The model outputs per-species per-tissue predictions but extrapolates poorly to tissues outside this set. Use the tissue closest to disease-relevant context. **For tissues not in the 4-tissue training set, fall back to SpliceAI** — Pangolin extrapolates poorly to unseen tissues.


## SpliceVault for Empirical Mis-Splicing Outcomes


**Goal:** Predict the *type* of mis-splicing (exon skipping vs cryptic site activation) given a canonical-disrupting variant.

**Approach:** Query SpliceVault's database of empirical mis-splicing events from public RNA-seq.

```python
import requests

# Web API: https://kidsneuro.shinyapps.io/splicevault/
# Or use the R/Python package at github.com/kidsneuro-lab/SpliceVault

# Example: NM_000546.6:c.673-2A>G (TP53)
# Returns top-N most likely mis-splicing events: exon skipping, cryptic 3'ss usage, etc.
```

SpliceVault (Dawes 2023 *Nat Genet*) showed that the **Top-4 events** at any splice site explain >95% of empirical mis-splicing — a striking regularity that makes consequence prediction tractable. Use SpliceVault when the question is not "will splicing change?" but "what specific aberrant splicing will occur?".


## MMSplice for Calibrated ΔPSI


**Goal:** Predict quantitative ΔPSI (not just probability of disruption) for cassette exons.

**Approach:** Score variant impact on each splicing region (5'ss, 3'ss, exon, intron-3'/5') and combine.

```python
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save

dl = SplicingVCFDataloader(
    gtf='gencode.v45.basic.gtf',
    fasta_file='GRCh38.fa',
    vcf_file='input.vcf'
)

model = MMSplice()
predict_save(model, dl, 'mmsplice_predictions.csv', pathogenicity=True)
```

MMSplice (Cheng 2019 *Genome Biol*) reports Δlogit_psi per variant. Useful when calibrated effect sizes matter (research) more than probability of disruption (clinical screening). Companion **MTSplice** (Cheng 2021 *Genome Biol*) adds tissue-specific Δψ predictions.


## HGVS Splicing Nomenclature


Following den Dunnen 2016 *Hum Mutat*:

| Notation | Meaning |
|----------|---------|
| `c.123+1G>A` | +1 of intron downstream of exon ending at cDNA position 123 (canonical 5'ss G) |
| `c.123+5G>A` | +5 position of donor (consensus region) |
| `c.124-1G>A` | -1 of acceptor (canonical AG) |
| `c.124-3T>G` | -3 of acceptor (Py-tract / BPS region) |
| `c.124-50A>G` | Deep-intronic; may activate cryptic site |
| `r.123_456del` | RNA-level deletion (predicted exon skipping) |
| `r.spl?` | Unknown splice consequence |
| `r.0?` | No detectable RNA |
| `p.0?` | Unknown protein consequence |
| `p.(=)` | No predicted protein change (silent) |

Validation tools: VariantValidator (Freeman 2018 *Hum Mutat*), Mutalyzer 3 (Lefter 2021 *Hum Mutat*).


## Extended-Window Scoring for Deep-Intronic Variants


SpliceAI's default precomputed scores use a **50-nt window**, missing variants that create pseudoexons in deep intronic regions. For unsolved Mendelian cases:

```bash
# Recompute with extended window
spliceai -I input.vcf -O output_2kb.vcf -R genome.fa -A grch38 -D 2000

# Or use CI-SpliceAI (Strauch 2022 Bioinformatics) optimized for distal effects
```

| Window | Tradeoff |
|--------|----------|
| -D 50 (default) | Fast; captures canonical-site disruption; misses deep-intronic |
| -D 500 | Captures most pseudoexon-creating deep-intronic variants |
| -D 2000 | Maximum sensitivity; some false positives at large distances |

Pseudoexon creation in deep introns explains ~5-15% of unsolved Mendelian disease alleles in current cohorts (Smith 2024 *Nat Commun*). Disease examples: CFTR 3849+10kbC>T, USH2A c.7595-2144A>G, CEP290 c.2991+1655A>G (LCA10).


## Concordance Across Predictors


```python
import pandas as pd

merged = (spliceai_df
    .merge(pangolin_df, on=['chrom', 'pos', 'alt'], suffixes=('_sai', '_pang'))
    .merge(mmsplice_df, on=['chrom', 'pos', 'alt'])
)

merged['concordance'] = (
    (merged['delta_max_sai'] >= 0.2).astype(int) +
    (merged['pangolin_score'].abs() >= 0.2).astype(int) +
    (merged['delta_logit_psi'].abs() >= 1.0).astype(int)
)

merged['interpretation'] = merged['concordance'].map({
    0: 'concordant_benign',
    1: 'discordant_low_evidence',
    2: 'concordant_evidence',
    3: 'high_concordance_pathogenic'
})
```

| Concordance | Interpretation | Action |
|-------------|----------------|--------|
| 3/3 above threshold | High confidence | Report PP3 strong |
| 2/3 above | Concordant evidence | Report PP3 moderate |
| 1/3 above | Discordant | Report inconclusive; flag for RNA validation |
| 0/3 above | Concordant benign | BP4 supporting |

Discordance is the most informative pattern — variants where one model sees impact and others don't are high priority for RNA validation.


## Branchpoint Variant Detection


All current tools are **weak at branchpoint variants** because the BPS motif (yUnAy) has low information content. Specific branchpoint tools:

| Tool | Method | Notes |
|------|--------|-------|
| BPP | Position-weight matrix | Zhang 2017 *NAR* |
| LaBranchoR | Bidirectional LSTM | Paggi & Bejerano 2018 *Genome Biol* |
| SVM-BPfinder | SVM on conservation+sequence | Corvelo 2010 *PLoS Comput Biol* |
| BPHunter | Genome-wide branchpoint screen using GTEx-derived BP database | Zhang 2022 *PNAS* |

Branchpoint variants are under-recognized in clinical pipelines; SpliceAI captures only some because branchpoint motifs have low information content. **Recommendation:** when SpliceAI delta is borderline (0.1-0.3) for a variant in the BPS region (-18 to -40 from 3'ss), run BPHunter as supplement.


## Splice-Switching ASO Design


**Goal:** Design antisense oligonucleotides to modulate splicing therapeutically (e.g. SMA ISS-N1, DMD exon skipping).

**Approach:** Use SpliceAI to predict impact of binding-site occlusion; check accessibility (RNAfold); avoid SR/hnRNP off-target motifs.

```python
# Conceptual workflow - actual design uses ASO synthesis platforms
# 1. Identify target ESE/ESS/ISE/ISS region from MaxEntScan + SpliceAI scan
# 2. Design candidate 18-22 nt ASOs spanning the regulatory element
# 3. For each ASO, simulate splice-site occlusion impact via SpliceAI on the masked sequence
# 4. Filter for RNA accessibility (avoid stable hairpins) using RNAfold
# 5. Whole-transcriptome SpliceAI scan for off-target binding (>=17/20 nt match)
# 6. Avoid TLR9 immunostimulatory CpG motifs

# Chemistry choices:
# - 2'-MOE-PS: nusinersen-like (CNS, intrathecal)
# - PMO: DMD ASOs (systemic IV)
# - GalNAc-conjugated: hepatic targeting
```

Approved precedents: **nusinersen** (SMA ISS-N1 occlusion, exon 7 inclusion); **risdiplam** (small-molecule SMN2 splicing modulator); **eteplirsen/golodirsen/casimersen/viltolarsen** (DMD exon skipping). Design references: Hua 2008 *AJHG*; Aartsma-Rus 2023 *Nat Rev Drug Discov*.


## Population Database Lookup


| Database | Use for |
|----------|---------|
| gnomAD v4 | Allele frequency; SpliceAI annotations integrated |
| ClinVar | Existing classifications; SpliceAI integrated since 2020 |
| SpliceVarDB | Curated splice variants with experimental RNA validation |
| dbNSFP4 | Pre-computed splice scores aggregated |
| Recount3 | Tissue-specific PSI lookups from public RNA-seq |
| GTEx sQTL v8 | Tissue-specific splicing QTLs across 49 tissues |
| MaveDB | Splice MAVE results (e.g. BRCA1 saturation; Findlay 2018 *Nature*) |

Always check ClinVar first for existing classifications; cross-reference with gnomAD for population frequency before committing to PP3/PP4.

