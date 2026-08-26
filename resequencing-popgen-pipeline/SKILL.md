---
name: resequencing-popgen-pipeline
description: "End-to-end plant resequencing population genetics pipeline: fastp QC, minibwa mapping, freebayes joint calling, VCF norm/filter, plink2 LD+PCA, ADMIXTURE, FastTree2, ANGSD pi/TajD/Fst, ggtree visualization (5-figure template), HTML report. Okabe-Ito colorblind palette."
category: bioinformatics
triggers:
  - "重测序分析"
  - "resequencing popgen"
  - "群体遗传全流程"
  - "群体结构三剑客"
  - "freebayes joint calling pipeline"
---

# Plant Resequencing Population Genetics Pipeline

From raw FASTQ → VCF → population structure (5-figure template) + Fst/π/TajD/LD decay → HTML report.

## Environment

- **fb_calling** conda env: freebayes, bcftools, vcftools
- **pop_genetics** conda env: plink2, admixture, fasttree, angsd
- **R** packages: ggplot2, ggtree, ggtreeExtra, ggnewscale, ape, patchwork
- **minibwa**: fast BWA-MEM successor (~/bin/minibwa)

## Directory Structure

```
project/req/
├── 01.rawdata/          # Raw FASTQ files
├── 02.cleandata/        # fastp-trimmed reads
├── 03.Mapping/          # minibwa alignment + samtools sort
│   ├── ref/             # Reference genome (genome.fasta symlink)
│   ├── *.bam            # Sorted BAM files
│   └── summary/         # Mapping statistics
├── 04.Variants/         # freebayes joint calling + VCF filtering
├── 05.population/       # Main analysis (all samples)
├── 06.pv25_only/        # Subgroup analysis (optional)
├── 07.selection/        # Selection scans (Fst, π, TajD, LD decay)
└── 00.report/           # HTML report + figures
```

## Pipeline Stages

### Stage 1: QC (fastp)
```bash
fastp -i R1.fq.gz -I R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
  --thread 8 --html sample_fastp.html --json sample_fastp.json
```

### Stage 2: Mapping (minibwa + samtools)
```bash
minibwa map -t 8 -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
  ref/genome.fasta clean_R1.fq.gz clean_R2.fq.gz | \
  samtools sort -@ 4 -o sample.bam -
samtools index -c sample.bam   # Use -c for CSI (large genomes >500Mb)
```

**Note**: Use CSI index (not BAI) for genomes with chromosomes >500Mb (e.g., pea chr03=553Mb). BAI limit is 2^29 = 536Mb.

### Stage 3: Variant Calling (freebayes joint calling)
```bash
# Generate BAM list
ls *.bam > bamlist.txt

# Joint calling
freebayes -f ref/genome.fasta -L bamlist.txt \
  -p 1 -P 0 -E 3 --min-coverage 4 \
  --use-best-n-alleles 5 \
  > raw_variants.vcf.gz
```

### Stage 4: VCF Filtering

| Parameter | Threshold | Rationale |
|-----------|-----------|-----------|
| QUAL | ≥20 | Base quality filter |
| GQ | ≥20 | Genotype quality |
| DP | [4, 1000] | Remove low/high coverage artifacts |
| F_MISSING | ≤0.2 | Max 20% missing genotypes |
| MAF | ≥0.05 | Minor allele frequency |
| Alleles | bi-allelic SNP only | Remove indels/multiallelic |

```bash
bcftools view -i 'QUAL>=20' raw_variants.vcf.gz | \
bcftools view -i 'F_MISSING<=0.2' | \
bcftools view -m2 -M2 -v snps | \
bcftools view -i 'FORMAT/DP>=4 && FORMAT/DP<=1000 && FORMAT/GQ>=20' | \
bcftools view -i 'MAF>=0.05' > filtered.vcf
```

### Stage 5: Population Structure

**PCA (plink2):**
```bash
plink2 --vcf filtered.vcf --set-all-var-ids @:# --rm-dup \
  --pca 10 --out pca_results
```

**ADMIXTURE:**
```bash
# Convert to bed/bim/fam first
plink2 --vcf filtered.vcf --set-all-var-ids @:# --rm-dup --make-bed --out pop

# Run K=2 to K=6
for K in 2 3 4 5 6; do
  admixture --cv pop.bed $K -j8 | tee log${K}.out
done
```

**Neighbor-Joining Tree (FastTree2):**
```bash
# Generate SNP alignment from VCF
vcf2phylip --input filtered.vcf --output aligned.phy
FastTree2 -gtr -nt aligned.phy > tree.nwk
```

### Stage 6: Selection Scans

**Fst (vcftools, window-based):**
```bash
vcftools --vcf filtered.vcf --weir-fst-pop popA.txt --weir-fst-pop popB.txt \
  --fst-window-size 50000 --fst-window-step 10000 --out fst_results
```

**Nucleotide Diversity π (vcftools):**
```bash
vcftools --vcf filtered.vcf --keep popA.txt \
  --window-pi 50000 --window-pi-step 10000 --out pi_popA
```

**Tajima's D (vcftools):**
```bash
vcftools --vcf filtered.vcf --keep popA.txt \
  --TajimaD 50000 --out tajd_popA
```

**LD Decay (PopLDdecay):**
```bash
PopLDdecay -InVCF filtered.vcf -OutStat ld_stat \
  -SubPop popA.txt -MAF 0.05 -Miss 0.2
perl Plot_OnePop.pl -inFile ld_stat.stat.gz -output ld_decay_plot
```

## Visualization (5-Figure Template)

**Fig1:** Circular tree (clade-colored, upward propagation)
**Fig2:** Circular tree + K=2 Admixture outer ring (Q matrix aligned to tree tip order)
**Fig3:** PCA (PC1/2 + PC1/3, population-colored)
**Fig4:** Admixture K=2-4 (sorted by ancestry proportion within groups)
**FigCV:** Cross-validation error curve

### Key R Code Pattern
```r
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(patchwork)

# Circular tree with clade coloring
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tiplab2(size = 2.5, align = TRUE, linesize = 0.3)

# Add Admixture ring
p <- p + new_scale_fill() +
  geom_fruit(data = q_df, geom = geom_bar,
             mapping = aes(y = tip, fill = factor(K), x = 1),
             stat = "identity", width = 0.8, offset = 0.1)
```

## Critical Pitfalls

1. **plink2 LD**: Must use `--set-all-var-ids @:# --rm-dup` or get duplicate ID errors
2. **ADMIXTURE**: Integer chromosome names; .Q files must be in current working directory
3. **ANGSD `-fold`**: Only use with `realSFS`, not standalone
4. **ANGSD 2D SFS**: NOT feasible with >140 samples — use vcftools `--weir-fst-pop` instead (minutes vs hangs)
5. **PCA eigenvalues**: Divide by sum of all eigenvalues for variance explained
6. **Q matrix alignment**: `Q[match(tree$tip.label, fam_order), ]` — must match tree tip order
7. **Clade coloring**: Propagate upward (child → parent), not majority vote
8. **Long tasks**: Use `nohup cmd > log 2>&1 &` to avoid wrapper timeouts
9. **PopLDdecay**: Full path `/media/nfs1/hermes/miniforge3/bin/PopLDdecay`; use `comment.char=""` to read `.stat.gz`; avoid `-SubPop` (segfault)
10. **bcftools/PopLDdecay**: Not in pop_genetics env — use `conda activate fb_calling` or full path

## Citations

If you use this pipeline in your research, please cite the following tools:

### Core Tools

- **fastp**: Chen S, Zhou Y, Chen Y, Gu J. "fastp: an ultra-fast all-in-one FASTQ preprocessor." *Bioinformatics*, 2018, 34(17):884-890. doi:10.1093/bioinformatics/bty560
  ```bibtex
  @article{chen2018fastp,
    title={fastp: an ultra-fast all-in-one FASTQ preprocessor},
    author={Chen, Shifu and Zhou, Yanqing and Chen, Yaru and Gu, Jia},
    journal={Bioinformatics},
    volume={34},
    number={17},
    pages={884--890},
    year={2018},
    publisher={Oxford University Press}
  }
  ```

- **minibwa**: Li H. "New alignment algorithms for long reads and assembly." *bioRxiv*, 2024. doi:10.1101/2024.08.29.611371 (successor to BWA-MEM)
  ```bibtex
  @article{li2024minibwa,
    title={New alignment algorithms for long reads and assembly},
    author={Li, Heng},
    journal={bioRxiv},
    year={2024},
    publisher={Cold Spring Harbor Laboratory}
  }
  ```

- **samtools**: Danecek P, Bonfield JK, Liddle J, et al. "Twelve years of SAMtools and BCFtools." *GigaScience*, 2021, 10(2):giab008. doi:10.1093/gigascience/giab008
  ```bibtex
  @article{danecek2021samtools,
    title={Twelve years of SAMtools and BCFtools},
    author={Danecek, Petr and Bonfield, James K and Liddle, John and Marshall, John and Ohan, Valeriu and Pollard, Martin O and Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A and Davies, Robert M and Li, Heng},
    journal={GigaScience},
    volume={10},
    number={2},
    pages={giab008},
    year={2021},
    publisher={Oxford University Press}
  }
  ```

- **freebayes**: Garrison E, Marth G. "Haplotype-based variant calling from next-generation sequencing datasets." *arXiv:1207.4127*, 2012.
  ```bibtex
  @article{garrison2012freebayes,
    title={Haplotype-based variant calling from next-generation sequencing datasets},
    author={Garrison, Erik and Marth, Gabor},
    journal={arXiv preprint arXiv:1207.4127},
    year={2012}
  }
  ```

- **bcftools**: See samtools citation above (Danecek et al. 2021).

- **vcftools**: Danecek P, Auton A, Abecasis G, et al. "The variant call format and VCFtools." *Bioinformatics*, 2011, 27(15):2156-2158. doi:10.1093/bioinformatics/btr330
  ```bibtex
  @article{danecek2011vcftools,
    title={The variant call format and VCFtools},
    author={Danecek, Petr and Auton, Adam and Abecasis, Goncalo and others},
    journal={Bioinformatics},
    volume={27},
    number={15},
    pages={2156--2158},
    year={2011},
    publisher={Oxford University Press}
  }
  ```

### Population Genetics Tools

- **plink2**: Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. "Second-generation PLINK: rising to the challenge of larger and richer datasets." *GigaScience*, 2015, 4(1):s13742-015-0047-8. doi:10.1186/s13742-015-0047-8
  ```bibtex
  @article{chang2015plink2,
    title={Second-generation PLINK: rising to the challenge of larger and richer datasets},
    author={Chang, Christopher C and Chow, Chow Chong and Tellier, Laurent CAM and Vattikuti, Shashank and Purcell, Shaun M and Lee, James J},
    journal={GigaScience},
    volume={4},
    number={1},
    pages={s13742-015-0047-8},
    year={2015},
    publisher={BioMed Central}
  }
  ```

- **ADMIXTURE**: Alexander DH, Novembre J, Lange K. "Fast model-based estimation of ancestry in unrelated individuals." *Genome Research*, 2009, 19(9):1655-1664. doi:10.1101/gr.094052.109
  ```bibtex
  @article{alexander2009admixture,
    title={Fast model-based estimation of ancestry in unrelated individuals},
    author={Alexander, David H and Novembre, John and Lange, Kenneth},
    journal={Genome Research},
    volume={19},
    number={9},
    pages={1655--1664},
    year={2009},
    publisher={Cold Spring Harbor Laboratory Press}
  }
  ```

- **FastTree2**: Price MN, Dehal PS, Arkin AP. "FastTree 2—approximately maximum-likelihood trees for large alignments." *PLoS ONE*, 2010, 5(3):e9490. doi:10.1371/journal.pone.0009490
  ```bibtex
  @article{price2010fasttree,
    title={FastTree 2—approximately maximum-likelihood trees for large alignments},
    author={Price, Morgan N and Dehal, Paramjit S and Arkin, Adam P},
    journal={PLoS ONE},
    volume={5},
    number={3},
    pages={e9490},
    year={2010},
    publisher={Public Library of Science}
  }
  ```

- **ANGSD**: Korneliussen TS, Albrechtsen A, Nielsen R. "ANGSD: Analysis of Next Generation Sequencing Data." *BMC Bioinformatics*, 2014, 15:356. doi:10.1186/s12859-014-0356-4
  ```bibtex
  @article{korneliussen2014angsd,
    title={ANGSD: Analysis of Next Generation Sequencing Data},
    author={Korneliussen, Thorfinn S and Albrechtsen, Anders and Nielsen, Rasmus},
    journal={BMC Bioinformatics},
    volume={15},
    pages={356},
    year={2014},
    publisher={BioMed Central}
  }
  ```

- **PopLDdecay**: Zhang C, Dong SS, Xu JY, He WM, Yang TL. "PopLDdecay: a fast and flexible tool for calculating decay of linkage disequilibrium in large population sets." *GigaScience*, 2019, 8(4):giy154. doi:10.1093/gigascience/giy154
  ```bibtex
  @article{zhang2019poplddecay,
    title={PopLDdecay: a fast and flexible tool for calculating decay of linkage disequilibrium in large population sets},
    author={Zhang, Cheng and Dong, Ssu-Shi and Xu, Jia-Yu and He, Wei-Ming and Yang, Tie-Lin},
    journal={GigaScience},
    volume={8},
    number={4},
    pages={giy154},
    year={2019},
    publisher={Oxford University Press}
  }
  ```

### Visualization

- **ggtree**: Yu G, Smith DK, Zhu H, Guan Y, Lam TTY. "ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data." *Methods in Ecology and Evolution*, 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
  ```bibtex
  @article{yu2017ggtree,
    title={ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data},
    author={Yu, Guangchuang and Smith, David K and Zhu, Huan and Guan, Yuting and Lam, Tommy Tsan-Yuk},
    journal={Methods in Ecology and Evolution},
    volume={8},
    number={1},
    pages={28--36},
    year={2017},
    publisher={Wiley}
  }
  ```

### Methods Template Sentence

> "Raw sequencing reads were quality-filtered with fastp (Chen et al. 2018) and aligned to the reference genome using minibwa (Li 2024). Variant calling was performed with freebayes (Garrison & Marth 2012) in joint-calling mode, followed by hard filtering with bcftools (Danecek et al. 2021) for quality (QUAL ≥ 20), genotype depth (DP 4–1000), and minor allele frequency (MAF ≥ 0.05). Population structure was assessed using PCA (plink2; Chang et al. 2015), model-based clustering (ADMIXTURE; Alexander et al. 2009), and neighbor-joining trees (FastTree2; Price et al. 2010). Population differentiation (Fst), nucleotide diversity (π), and Tajima's D were calculated with vcftools (Danecek et al. 2011), and linkage disequilibrium decay was estimated using PopLDdecay (Zhang et al. 2019). Visualizations were generated with ggtree (Yu et al. 2017) and ggplot2."
