## 23_finemap_susie.R ---- per-trait SuSiE fine-mapping of the main peak
## Usage: Rscript 23_finemap_susie.R fa|vee <trait>
## Requires: susieR (>= 0.14), data.table, ggplot2, patchwork
suppressMessages({library(data.table); library(susieR); library(ggplot2); library(patchwork)})
WDIR <- Sys.getenv("GWAS_WDIR", getwd())
CFG  <- Sys.getenv("GWAS_CONFIG", file.path(WDIR, "scripts/00_config.R"))
source(CFG)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; TR <- args[2]
if (prefix == "fa") {
  qtl_f <- qtl_file; CN <- unname(trait_cn[TR])
} else {
  qtl_f <- vee_qtl_file; CN <- unname(vee_cn[TR])
}
L <- as.integer(Sys.getenv("SUSIE_L", "10"))
HALF <- as.numeric(Sys.getenv("SUSIE_HALF_MB", "1")) * 1e6

OUT <- file.path(res_dir, paste0("susie_", TR)); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(fig_dir, paste0("susie_", TR)); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

qtl <- fread(qtl_f)[trait == TR]
qg <- qtl[!is.na(gene_id) & gene_id != "" & gene_id != "intergenic"]
anchor <- if (nrow(qg)) qg[which.max(lod)] else qtl[which.max(lod)]
CHR <- anchor$chr[1]; APOS <- anchor$pos[1]
FROM <- max(0, APOS - HALF); TO <- APOS + HALF
cat(sprintf("%s | anchor %s:%.2f Mb | region %.2f-%.2f Mb\n", TR, CHR, APOS/1e6, FROM/1e6, TO/1e6))

## ---- genotypes ----
pref <- file.path(data_dir, paste0("susie_", TR))
run_plink <- function() {
  for (g in c(0.1, 0.25, 0.5)) {
    cmd <- sprintf("%s --bfile %s --chr %s --from-bp %d --to-bp %d --maf %g --max-alleles 2 --snps-only --geno %g --make-just-bim --export A --out %s",
                   plink2, gt_bed, CHR, FROM, TO, MAF_MIN, g, pref)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (file.exists(paste0(pref, ".raw"))) return(TRUE)
  }
  FALSE
}
if (!run_plink()) stop("no SNPs in region after filters")
raw <- fread(paste0(pref, ".raw")); setnames(raw, 7:ncol(raw), sub("_[ACGT]+$", "", names(raw)[7:ncol(raw)]))
bim <- fread(paste0(pref, ".bim")); setnames(bim, c("chr", "snp", "cm", "pos", "a1", "a2"))
raw[, IID := as.character(IID)]
phe <- fread(trait_file); setnames(phe, 1, "IID"); phe[, IID := as.character(IID)]
pca <- fread(pca_file); setnames(pca, 1, "IID"); pca[, IID := as.character(IID)]
m <- merge(raw, phe, by = "IID"); m <- merge(m, pca, by = "IID")
sn <- bim$snp; X <- as.matrix(m[, ..sn])
cm <- colMeans(X, na.rm = TRUE); cm[is.nan(cm)] <- 0
for (j in seq_len(ncol(X))) { na <- is.na(X[, j]); if (any(na)) X[na, j] <- cm[j] }
y <- m[[TR]]; PCs <- as.matrix(m[, .(pc_1, pc_2, pc_3)])
R <- qr.resid(qr(cbind(1, PCs)), cbind(y, X)); y <- R[, 1]; X <- R[, -1, drop = FALSE]
X <- scale(X); y <- as.numeric(scale(y))

## ---- SuSiE ----
set.seed(1)
fit <- susie(X, y, L = L, estimate_residual_variance = TRUE,
             estimate_prior_variance = TRUE, verbose = FALSE, max_iter = 500)
pip <- susie_get_pip(fit)
cs <- susie_get_cs(fit, X = X, coverage = 0.95, min_abs_corr = 0.1)
sets <- cs$cs
cat(sprintf("SuSiE: %d credible sets | sigma2 %.3f\n", length(sets), fit$sigma2))

genes <- as.data.table(readRDS(file.path(data_dir, "gene_coords.rds")))[chr == CHR]
nearest_gene <- function(pos) {
  d <- pmin(abs(genes$start - pos), abs(genes$end - pos))
  genes$gene_id[which.min(d)]
}
cs_rows <- lapply(seq_along(sets), function(k) {
  idx <- sets[[k]]; top <- idx[which.max(pip[idx])]
  data.table(CS = k, size = length(idx), top_snp = sn[top], top_pos = bim$pos[top],
             top_pip = pip[top], mean_pip = mean(pip[idx]), top_gene = nearest_gene(bim$pos[top]))
})
CS <- if (length(cs_rows)) rbindlist(cs_rows) else data.table(CS = integer(), size = integer(), top_snp = character(), top_pos = integer(), top_pip = numeric(), mean_pip = numeric(), top_gene = character())
fwrite(CS, file.path(OUT, "credible_sets.tsv"), sep = "\t")

pipdt <- data.table(snp = sn, chr = CHR, pos = bim$pos, pip = pip)
pipdt[, cs := NA_integer_]
for (k in seq_along(sets)) pipdt[snp %in% sn[sets[[k]]], cs := k]
pipdt[, gene := vapply(pos, nearest_gene, character(1))]
fwrite(pipdt, file.path(OUT, "pip.tsv"), sep = "\t")

qtl2 <- qtl[chr == CHR & pos >= FROM & pos <= TO]
qtl2[, cs := NA_integer_]
for (k in seq_along(sets)) qtl2[marker %in% sn[sets[[k]]], cs := k]
qtl2[, in_cs := !is.na(cs)]
fwrite(qtl2, file.path(OUT, "qtl_overlap.tsv"), sep = "\t")

sink(file.path(OUT, "summary.txt"))
cat(sprintf("%s SuSiE fine-mapping | %s %.2f-%.2f Mb | SNPs %d | n %d | L %d\n",
            TR, CHR, FROM/1e6, TO/1e6, ncol(X), nrow(X), L))
cat(sprintf("credible sets: %d\n", length(sets)))
print(CS)
cat(sprintf("\nFast3VmrMLM QTN captured by CS: %d / %d\n", sum(qtl2$in_cs), nrow(qtl2)))
sink()

## ---- figure ----
mid_dir <- if (prefix == "fa") file.path(proj, "06.Fast3VmrMLM/results") else vee_mid_dir
mid <- fread(file.path(mid_dir, paste0(TR, "_midresult.csv")))
setnames(mid, c("snp", "chr", "pos", "wald", "p"))
mid <- mid[chr == CHR & pos >= FROM & pos <= TO][, neglog := -log10(pmax(p, 1e-300))]
d <- merge(mid[, .(snp, pos, neglog)], pipdt[, .(snp, pip, cs)], by = "snp")
d[, csg := factor(ifelse(is.na(cs), "非 CS", paste0("CS", cs)),
                  levels = c(paste0("CS", seq_along(sets)), "非 CS"))]
pal <- setNames(rep(c("#2F5597", "#C00000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
                      "#F0E442", "#8C564B", "#7F7F7F", "#17BECF"), length(sets)), paste0("CS", seq_along(sets)))
p1 <- ggplot(d, aes(pos / 1e6, neglog)) +
  geom_point(aes(colour = csg), size = 1.2, alpha = 0.8, stroke = 0) +
  geom_vline(data = qtl2, aes(xintercept = pos / 1e6), colour = "#B2182B", linetype = 3, linewidth = 0.4) +
  scale_colour_manual(values = c(pal, "非 CS" = "grey80"), name = "95% CS") +
  labs(subtitle = "红虚线 = Fast3VmrMLM QTN", x = NULL, y = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", plot.subtitle = element_text(size = 8.4))
p2 <- ggplot(d, aes(pos / 1e6, pip)) +
  geom_point(aes(colour = csg), size = 1.2, alpha = 0.85, stroke = 0) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = c(pal, "非 CS" = "grey85"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(x = sprintf("%s 位置 (Mb)", CHR), y = "PIP") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank())
p <- p1 / p2 + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(title = sprintf("%s SuSiE 精细定位（%s %.1f–%.1f Mb）", CN, CHR, FROM/1e6, TO/1e6),
                  theme = theme(plot.title = element_text(family = cjk_font, hjust = .5,
                                                          face = "bold", size = 12)))
ggsave(file.path(FIG, "susie_pip.png"), p, width = 9.6, height = 6.2, dpi = 220,
       device = grDevices::png, type = "cairo")
file.remove(paste0(pref, c(".raw", ".bim", ".log")))
cat("DONE 23\n")
