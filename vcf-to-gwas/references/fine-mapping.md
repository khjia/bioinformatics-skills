# Fine-mapping (SuSiE and friends)

Companion to `templates/23_finemap_susie.R`. Distilled from the bioSkills
`causal-genomics/fine-mapping` skill (MIT) and adapted to small-N, long-LD
plant/animal germplasm.

## Which method when

| Situation | Method | Notes |
|---|---|---|
| Individual genotypes + phenotype available | `susieR::susie(X, y, L = 10)` | Default in this skill; in-sample LD is exact |
| Summary statistics + matched LD reference | `susie_rss(z, R, n, L = 10)` | Run `estimate_s_rss()` first |
| Locus has a polygenic shoulder | SuSiE-inf / FINEMAP-inf | Separate packages; calibrates non-sparse loci |
| Functional priors (coding, conserved, pathway) | `susie_rss(..., prior_weights = w)` | `w` is a per-SNP probability, not a variance |
| Multi-ancestry cohorts | SuSiEx / MultiSuSiE | Out of scope here |

**Low N (this project, n ≈ 200).** The default 95% credible set is often very
wide. Report `coverage = 0.9` (or 0.8) explicitly in methods; never silently
change coverage. PIPs at low N are optimistic — treat them as a ranking.

## LD diagnostics (summary-statistics route)

```r
s_hat <- susieR::estimate_s_rss(z, R, n)     # global Z/LD inconsistency
# lambda > 0.05 => LD reference does not match the GWAS sample
kriging <- susieR::kriging_rss(z, R, n)      # per-SNP residuals; |z_obs - z_exp| > 3 flags problems
```

Common causes: different cohort for LD, strand/allele flips, wrong genome
build, or an LD matrix from a different ancestry. Fix the reference before
interpreting credible sets. With individual-level `susie()` this diagnostic is
unnecessary (LD is in-sample).

## Reading the output

- A size-1 credible set is a *model choice*, not proof: with long LD and
  moderate N, SuSiE can be overconfident (PIP ≈ 1 while a correlated variant
  is equally plausible). Check `r²` between the top variant and its neighbours.
- Credible sets are per signal; the number of sets ≈ number of independent
  signals the model supports. A signal below the model's detection threshold
  (e.g. p ≈ 3 × 10⁻⁶ at n = 200) can be missed.
- Compare with the multi-locus QTN list (Fast3VmrMLM): QTN are LD proxies;
  SuSiE collapses them to representatives. Report both.

## Pitfalls

1. **`prior_weights` vs `prior_variance`** — per-SNP probability vs per-effect
   variance. Passing the wrong one is silently accepted and meaningless.
2. **Too small `L`** — if the locus has more signals than `L`, extra effects
   are absorbed and PIPs distort. Re-run with `L = 5, 10, 15` and check
   stability (this skill's default check).
3. **LD mismatch** — see diagnostics above.
4. **Long-range LD / tandem gene clusters** — e.g. KCS/FAE1 arrays; check
   whether the credible set spans several paralogs before naming one causal gene.
5. **PIP ≠ posterior probability at low N** — calibrate expectations; report
   the credible set, not just the top SNP.

## Minimal reproducible recipe (individual-level)

```r
R  <- qr.resid(qr(cbind(1, PCs)), cbind(y, X))   # remove PC effects
y  <- scale(R[, 1]); X <- scale(R[, -1])
fit <- susie(X, y, L = 10, estimate_residual_variance = TRUE)
cs  <- susie_get_cs(fit, X = X, coverage = 0.9, min_abs_corr = 0.1)
pip <- susie_get_pip(fit)
```

Record: `susieR` version, `L`, coverage, `min_abs_corr`, region boundaries,
number of SNPs, and whether LD was in-sample.
