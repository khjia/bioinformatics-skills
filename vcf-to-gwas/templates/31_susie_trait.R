## 31_susie_trait.R ---- generic SuSiE fine-mapping for one trait
## Usage: Rscript 31_susie_trait.R fa|vee <trait>
suppressMessages({library(data.table); library(susieR); library(ggplot2); library(patchwork); library(ggrepel)})
WDIR <- Sys.getenv("GWAS_WDIR", getwd())
CFG  <- Sys.getenv("GWAS_CONFIG", file.path(WDIR, "scripts/00_config.R"))
source(CFG)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; TR <- args[2]
if (prefix == "fa") {
  phe_f <- trait_file; pca_f <- pca_file; mid_dir <- mid_dir_fa; qtl_f <- qtl_file
  CN <- unname(trait_cn[TR])
} else {
  phe_f <- vee_trait_file; pca_f <- vee_pca_file; mid_dir <- vee_mid_dir; qtl_f <- vee_qtl_file
  CN <- unname(vee_cn[TR])
}
cat(sprintf("== SuSiE %s %s (%s)\n", prefix, TR, CN))

OUT <- file.path(res_dir, paste0("susie_", TR)); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(fig_dir, paste0("susie_", TR)); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## ---- window from the report summary (anchor +/- 1 Mb) ----
sf <- file.path(res_dir, paste0(prefix, "_", TR, "_report"), "tab_summary.tsv")
if (!file.exists(sf)) { cat("no tab_summary, skip\n"); quit(status = 0) }
st <- fread(sf)
CHR <- as.character(st$anchor_chr[1]); APOS <- as.numeric(st$anchor_pos[1])
if (is.na(CHR) || !nzchar(CHR) || is.na(APOS)) { cat("no anchor, skip\n"); quit(status = 0) }
FROM <- max(1, APOS - 1e6); TO <- APOS + 1e6
cat(sprintf("region %s:%.3f-%.3f Mb\n", CHR, FROM/1e6, TO/1e6))

## ---- genotype ----
pref <- file.path(data_dir, paste0("susie_", TR, "_region"))
get_geno <- function() {
  for (maf in c(MAF_MIN, 0.02)) for (g in c(0.25, 0.5)) {
    cmd <- sprintf("%s --bfile %s --chr %s --from-bp %d --to-bp %d --maf %g --max-alleles 2 --snps-only --geno %g --make-just-bim --export A --out %s",
                   plink2, gt_bed, CHR, FROM, TO, maf, g, pref)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (file.exists(paste0(pref, ".raw")) && file.info(paste0(pref, ".raw"))$size > 1000) return(TRUE)
  }
  FALSE
}
if (!get_geno()) { cat("no genotype, skip\n"); quit(status = 0) }
raw <- fread(paste0(pref, ".raw")); setnames(raw, 7:ncol(raw), sub("_[ACGT]+$", "", names(raw)[7:ncol(raw)]))
bim <- fread(paste0(pref, ".bim")); setnames(bim, c("chr", "snp", "cm", "pos", "a1", "a2"))
raw[, IID := as.character(IID)]

phe <- fread(phe_f); setnames(phe, 1, "IID"); phe[, IID := as.character(IID)]
pca <- fread(pca_f); setnames(pca, 1, "IID"); pca[, IID := as.character(IID)]
m <- merge(raw, phe[, c("IID", TR), with = FALSE], by = "IID")
m <- merge(m, pca, by = "IID")
snps_all <- bim$snp
X <- as.matrix(m[, ..snps_all])
cm <- colMeans(X, na.rm = TRUE); cm[is.nan(cm)] <- 0
for (j in seq_len(ncol(X))) { na <- is.na(X[, j]); if (any(na)) X[na, j] <- cm[j] }
y <- m[[TR]]
pcn <- grep("^pc_", names(m), value = TRUE)
PCs <- as.matrix(m[, ..pcn])

## cap SNP number so SuSiE stays tractable (uniform thinning, keeps order)
MAXS <- 3000
if (ncol(X) > MAXS) {
  keep <- unique(round(seq(1, ncol(X), length.out = MAXS)))
  X <- X[, keep, drop = FALSE]; bim <- bim[keep]; snps_all <- bim$snp
}
cat(sprintf("SNPs %d | samples %d | PCs %d\n", ncol(X), nrow(X), ncol(PCs)))

ok <- !is.na(y)
R <- qr.resid(qr(cbind(1, PCs[ok, , drop = FALSE])), cbind(y[ok], X[ok, , drop = FALSE]))
yy <- as.numeric(scale(R[, 1])); XX <- scale(R[, -1, drop = FALSE])
XX[is.na(XX)] <- 0

## ---- SuSiE ----
L <- if (ncol(XX) < 300) 5 else 10
fit <- tryCatch({
  set.seed(1); t0 <- Sys.time()
  f <- susie(XX, yy, L = L, estimate_residual_variance = TRUE,
             estimate_prior_variance = TRUE, verbose = FALSE, max_iter = 300)
  cat(sprintf("SuSiE done in %.1f min (L=%d)\n", difftime(Sys.time(), t0, units = "mins"), L)); f
}, error = function(e) { cat("SuSiE failed:", conditionMessage(e), "\n"); NULL })
if (is.null(fit)) { cat("skip\n"); quit(status = 0) }

## weak/diffuse signal can drive every prior variance to 0 -> all PIP = 0;
## in that case refit with a fixed prior variance so PIPs stay informative
pip <- susie_get_pip(fit)
if (max(pip, na.rm = TRUE) < 1e-6) {
  cat("all prior variances 0 -> refit with fixed scaled_prior_variance = 0.2\n")
  fit <- tryCatch(susie(XX, yy, L = L, estimate_prior_variance = FALSE,
                        scaled_prior_variance = 0.2, estimate_residual_variance = TRUE,
                        verbose = FALSE, max_iter = 300), error = function(e) fit)
  pip <- susie_get_pip(fit)
}
cs  <- susie_get_cs(fit, X = XX, coverage = 0.95, min_abs_corr = 0.1)
sets <- cs$cs

## if purity filtering removed every set, relax it once and record the size:
## a credible set spanning most of the window means the signal cannot be localised
relaxed <- NULL
if (!length(sets)) {
  rc <- susie_get_cs(fit, X = XX, coverage = 0.90, min_abs_corr = 0)
  if (length(rc$cs)) {
    relaxed <- data.table(CS = seq_along(rc$cs), size = lengths(rc$cs),
                          frac_window = lengths(rc$cs) / ncol(XX))
    fwrite(relaxed, file.path(OUT, "relaxed_cs.tsv"), sep = "\t")
    cat(sprintf("relaxed CS: %d sets, sizes %s\n", nrow(relaxed),
                paste(relaxed$size, collapse = ",")))
  }
}

genes <- as.data.table(readRDS(file.path(data_dir, "gene_coords.rds")))
nearest_gene <- function(pos) {
  d <- pmin(abs(genes$start - pos), abs(genes$end - pos))
  genes$gene_id[which.min(d)]
}

cs_rows <- list()
for (k in seq_along(sets)) {
  idx <- sets[[k]]
  top <- idx[which.max(pip[idx])]
  cs_rows[[k]] <- data.table(
    CS = k, size = length(idx), top_snp = snps_all[top], top_pos = bim$pos[top],
    top_pip = pip[top], mean_pip = mean(pip[idx]), min_pip = min(pip[idx]),
    top_gene = nearest_gene(bim$pos[top]),
    snps = paste(snps_all[idx], collapse = ";"),
    positions = paste(bim$pos[idx], collapse = ";"))
}
CS <- if (length(cs_rows)) rbindlist(cs_rows) else
  data.table(CS = integer(), size = integer(), top_snp = character(), top_pos = numeric(),
             top_pip = numeric(), mean_pip = numeric(), min_pip = numeric(),
             top_gene = character(), snps = character(), positions = character())
fwrite(CS, file.path(OUT, "credible_sets.tsv"), sep = "\t")

pipdt <- data.table(snp = snps_all, chr = CHR, pos = bim$pos, pip = pip)
pipdt[, cs := NA_integer_]
for (k in seq_along(sets)) pipdt[snp %in% snps_all[sets[[k]]], cs := k]
pipdt[, gene := vapply(pos, nearest_gene, character(1))]
fwrite(pipdt, file.path(OUT, "pip.tsv"), sep = "\t")

qtl <- fread(qtl_f)[trait == TR & chr == CHR & pos >= FROM & pos <= TO]
qtl[, cs := NA_integer_]
for (k in seq_along(sets)) qtl[marker %in% snps_all[sets[[k]]], cs := k]
qtl[, in_cs := !is.na(cs)]
fwrite(qtl, file.path(OUT, "qtl_overlap.tsv"), sep = "\t")

sink(file.path(OUT, "summary.txt"))
cat(sprintf("SuSiE: %s %s %s %.2f-%.2f Mb\n", prefix, TR, CHR, FROM/1e6, TO/1e6))
cat(sprintf("SNPs %d | samples %d | L %d | CS %d\n", ncol(XX), nrow(XX), L, length(sets)))
if (!is.null(relaxed)) cat(sprintf("relaxed CS (90%%, no purity): %d sets, sizes %s, frac %s\n", nrow(relaxed), paste(relaxed$size, collapse = ","), paste(sprintf("%.2f", relaxed$frac_window), collapse = ",")))
if (nrow(CS)) print(CS[, .(CS, size, top_snp, top_pos, top_pip = round(top_pip, 3), top_gene)])
cat(sprintf("\nQTN captured: %d / %d\n", sum(qtl$in_cs), nrow(qtl)))
sink()

## ---- figure ----
mid <- fread(file.path(mid_dir, paste0(TR, "_midresult.csv")))
setnames(mid, c("snp", "chr", "pos", "wald", "p"))
mid <- mid[chr == CHR & pos >= FROM & pos <= TO]
mid[, neglog := -log10(pmax(p, 1e-300))]
d <- merge(mid[, .(snp, pos, neglog)], pipdt[, .(snp, pip, cs)], by = "snp")
d[, csg := factor(ifelse(is.na(cs), "非 CS", paste0("CS", cs)),
                  levels = unique(c(paste0("CS", seq_along(sets)), "非 CS")))]
pal <- c("#2F5597", "#C00000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#8C564B")
cols <- c("非 CS" = "grey85")
if (length(sets)) cols <- c(setNames(rep(pal, length.out = length(sets)), paste0("CS", seq_along(sets))),
                            "非 CS" = "grey85")

p1 <- ggplot(d, aes(pos/1e6, neglog)) +
  geom_point(aes(colour = csg), size = 1.2, alpha = 0.8, stroke = 0) +
  geom_vline(data = qtl, aes(xintercept = pos/1e6), colour = "#B2182B", linetype = 3, linewidth = 0.4) +
  scale_colour_manual(values = cols, name = "SuSiE 95% CS") +
  labs(subtitle = "红虚线 = Fast3VmrMLM QTN", x = NULL, y = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", plot.subtitle = element_text(size = 8.4))

p2 <- ggplot(d, aes(pos/1e6, pip)) +
  geom_point(aes(colour = csg), size = 1.2, alpha = 0.85, stroke = 0) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(x = sprintf("%s 位置 (Mb)", CHR), y = "PIP") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank())

p <- p1 / p2 + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(title = sprintf("%s SuSiE 精细定位（%s %.2f–%.2f Mb）", CN, CHR, FROM/1e6, TO/1e6),
                  theme = theme(plot.title = element_text(family = cjk_font, hjust = .5,
                                                          face = "bold", size = 12)))
ggsave(file.path(FIG, "susie_pip.png"), p, width = 9.6, height = 6.2, dpi = 220,
       device = grDevices::png, type = "cairo")
file.remove(paste0(pref, c(".raw", ".bim", ".log")))
cat("DONE 31\n")
