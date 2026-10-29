## 20_trait_report_figs.R ---- per-trait single-trait report figures & tables
## Usage: Rscript 20_trait_report_figs.R fa|vee <trait>
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
WDIR <- Sys.getenv("GWAS_WDIR", getwd())
CFG  <- Sys.getenv("GWAS_CONFIG", file.path(WDIR, "scripts/00_config.R"))
source(CFG)

args   <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
TR     <- args[2]

if (prefix == "fa") {
  pheno_f <- trait_file
  pca_f   <- pca_file
  mid_dir <- mid_dir_fa
  qtl_f   <- qtl_file
  CN      <- unname(trait_cn[TR])
} else {
  pheno_f <- vee_trait_file
  pca_f   <- vee_pca_file
  mid_dir <- vee_mid_dir
  qtl_f   <- vee_qtl_file
  CN      <- unname(vee_cn[TR])
}

OUT <- file.path(res_dir, paste0(prefix, "_", TR, "_report")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(fig_dir, paste0(prefix, "_", TR, "_report")); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

BONF <- 0.05 / if (exists("BONF_N")) BONF_N else 2494579
COL  <- c("奇数染色体" = "#2F5597", "偶数染色体" = "#9DC3E6")

## ---------------- inputs ----------------
qtl  <- fread(qtl_f)[trait == TR]
mid  <- fread(file.path(mid_dir, paste0(TR, "_midresult.csv")))
setnames(mid, c("snp", "chr", "pos", "wald", "p"))
mid  <- mid[!is.na(p) & !is.na(pos)]
phe  <- fread(pheno_f); y <- phe[[TR]]
pca  <- fread(pca_f); setnames(pca, 1, "IID"); pca[, IID := as.character(IID)]
anno <- fread(anno_file, select = c("gene", "Preferred_name", "EC", "Pfam", "A.thaliana",
                                    "eggNOG free text description"))
genes <- as.data.table(readRDS(file.path(data_dir, "gene_coords.rds")))

## anchor QTN: prefer gene-annotated QTN with max LOD
qg <- qtl[!is.na(gene_id) & gene_id != "" & gene_id != "intergenic"]
anchor <- if (nrow(qg)) qg[which.max(lod)] else qtl[which.max(lod)]
if (!nrow(anchor)) {
  top <- mid[which.min(p)]
  anchor <- data.table(trait = TR, marker = top$snp, chr = top$chr, pos = top$pos,
                       lod = NA_real_, significance = NA_character_,
                       gene_id = "intergenic", gene_name = "", dist = NA_integer_)
}
APOS <- anchor$pos[1]; ACHR <- anchor$chr[1]; AMARK <- anchor$marker[1]
WIN1 <- APOS - 1e6; WIN2 <- APOS + 1e6
cat(sprintf("trait %s | anchor %s (%s:%.2f Mb) | window %.2f-%.2f Mb\n",
            TR, AMARK, ACHR, APOS/1e6, WIN1/1e6, WIN2/1e6))

## ---------------- 1. phenotype ----------------
ys <- data.table(y = y)
pst <- data.table(n = length(y), mean = mean(y), sd = sd(y), cv = sd(y) / mean(y),
                  min = min(y), q1 = quantile(y, .25), median = median(y),
                  q3 = quantile(y, .75), max = max(y), shapiro_p = shapiro.test(y)$p.value)
fwrite(pst, file.path(OUT, "tab_pheno.tsv"), sep = "\t")

dy <- density(y)
p1 <- ggplot(ys, aes(y)) +
  geom_histogram(aes(y = after_stat(density), fill = "频数直方图"), bins = 24,
                 colour = "white", linewidth = 0.4) +
  geom_density(aes(colour = "核密度曲线"), linewidth = 0.9) +
  geom_vline(aes(xintercept = mean(y), colour = "均值"), linetype = 2, linewidth = 0.7) +
  annotate("text", x = mean(y), y = 0.95 * max(dy$y), vjust = 0, hjust = -0.05,
           label = sprintf("均值 = %.3f", mean(y)), colour = "#C00000", size = 3.2) +
  scale_fill_manual(values = c("频数直方图" = "#9DC3E6"), name = "分布") +
  scale_colour_manual(values = c("核密度曲线" = "#2F5597", "均值" = "#C00000"), name = "叠加元素") +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  labs(title = sprintf("%s 表型分布（n = %d）", CN, length(y)),
       subtitle = sprintf("SD = %.3f   CV = %.1f%%   范围 = %.2f–%.2f   Shapiro-Wilk p = %.2f",
                          sd(y), 100 * sd(y) / mean(y), min(y), max(y), shapiro.test(y)$p.value),
       x = sprintf("%s 含量", CN), y = "密度") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        legend.position = "top", plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5, size = 9))
ggsave(file.path(FIG, "f1_pheno.png"), p1, width = 7.2, height = 4.4, dpi = 220,
       device = grDevices::png, type = "cairo")

## ---------------- 2. Manhattan ----------------
mid[, neglog := -log10(pmax(p, 1e-300))]
mid[, bin := pos %/% 1000L]
db <- mid[, .(pos = as.numeric(median(pos)), neglog = max(neglog)), keyby = .(chr, bin)]
chrs <- unique(mid$chr)
db[, chr := factor(chr, levels = chrs)]
lens <- db[, .(len = max(pos)), keyby = chr][match(chrs, chr)]
lens[, off := c(0, head(cumsum(len + 5e6), -1))]
db <- merge(db, lens[, .(chr, off)], by = "chr", sort = FALSE)
db[, cum := pos + off]
ctr <- db[, .(c = mean(range(cum))), keyby = chr][match(chrs, chr)]
db[, grp := factor(ifelse(as.integer(chr) %% 2 == 1, "奇数染色体", "偶数染色体"), levels = names(COL))]

qt <- copy(qtl); qt[, chr := factor(chr, levels = chrs)]
qt <- merge(qt, lens[, .(chr, off)], by = "chr", sort = FALSE)
qt[, cum := pos + off]
qt[, neglog := -log10(pmax(mid$p[match(marker, mid$snp)], 1e-300))]
qt[, lab := fifelse(!is.na(gene_id) & gene_id != "intergenic" & gene_id != "", gene_id,
                    sprintf("%s %.2fMb", sub("^chr", "", as.character(chr)), pos/1e6))]
topl <- qt[order(-lod)][, .SD[1], by = lab][order(-lod)][1:4]

p2 <- ggplot(db, aes(cum, neglog)) +
  geom_point(aes(colour = grp), size = 1.35, alpha = 0.7, stroke = 0) +
  geom_point(data = qt, aes(cum, neglog, shape = significance), colour = "#C00000",
             size = 2.7, stroke = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = -log10(BONF), colour = "#FF0000", linetype = "dotted", linewidth = 0.55) +

  geom_text_repel(data = topl, aes(cum, neglog, label = lab), inherit.aes = FALSE,
                  size = 2.7, colour = "#7F0000", min.segment.length = 0,
                  max.overlaps = Inf, segment.size = 0.25, direction = "y",
                  nudge_y = 6, seed = 42) +
  scale_colour_manual(values = COL, name = "染色体") +
  scale_shape_manual(values = c(SIG = 17, SUG = 16), name = "QTN 显著性") +
  scale_x_continuous(breaks = ctr$c, labels = sub("^chr", "", ctr$chr), expand = c(0.006, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = sprintf("%s（%s）全基因组关联 Manhattan 图", CN, TR),
       subtitle = sprintf("红三角 = Fast3VmrMLM 多位点模型检出 QTN（SIG %d 个，SUG %d 个），纵坐标为相应标记的单点 −log₁₀P；红点线 = Bonferroni 2.0×10⁻⁸",
                          sum(qt$significance == "SIG", na.rm = TRUE), sum(qt$significance == "SUG", na.rm = TRUE)),
       x = NULL, y = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = cjk_font), panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.text.x = element_text(size = 6.5), legend.position = "top",
        legend.box = "horizontal", legend.key.width = unit(1.1, "lines"),
        plot.title = element_text(hjust = .5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = .5, size = 8.6),
        plot.margin = margin(6, 10, 4, 6))
ggsave(file.path(FIG, "f2_manhattan.png"), p2, width = 10, height = 4.6, dpi = 220,
       device = grDevices::png, type = "cairo")

## ---------------- 3. QQ ----------------
n  <- nrow(mid)
qq <- data.table(obs = -log10(sort(mid$p, decreasing = FALSE)), exp = -log10(ppoints(n)))
lam <- median(qchisq(1 - mid$p, 1)) / qchisq(0.5, 1)
p3 <- ggplot(qq, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_point(size = 0.6, alpha = 0.45, colour = "#2F5597", stroke = 0) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, max(qq$obs) * 1.02)) +
  labs(title = sprintf("%s 分位-分位（QQ）图", CN),
       subtitle = sprintf("λ_GC = %.3f（基于 %s 个标记）", lam, format(n, big.mark = ",")),
       x = expression(期望~-log[10](italic(p))), y = expression(观测~-log[10](italic(p)))) +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5, size = 9))
ggsave(file.path(FIG, "f3_qq.png"), p3, width = 4.6, height = 4.4, dpi = 220,
       device = grDevices::png, type = "cairo")

## ---------------- 4. QTN per chromosome ----------------
ch <- qt[, .(SIG = sum(significance == "SIG", na.rm = TRUE),
             SUG = sum(significance == "SUG", na.rm = TRUE)), by = chr]
ch <- ch[order(chr)]
chl <- melt(ch, id.vars = "chr", variable.name = "type", value.name = "n")
chl[, type := factor(type, levels = c("SIG", "SUG"), labels = c("SIG（LOD ≥ 7.7）", "SUG（3 < LOD < 7.7）"))]
p4 <- ggplot(chl[n > 0], aes(chr, n, fill = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.6, colour = "white") +
  scale_fill_manual(values = c("SIG（LOD ≥ 7.7）" = "#2F5597", "SUG（3 < LOD < 7.7）" = "#9DC3E6"), name = "QTN 类型") +
  scale_x_discrete(labels = function(x) sub("^chr", "", x)) +
  labs(title = sprintf("%s 检出 QTN 的染色体分布（共 %d 个）", CN, nrow(qt)),
       x = "染色体", y = "QTN 数量") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = cjk_font), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(hjust = .5, face = "bold"))
ggsave(file.path(FIG, "f4_qtn_chr.png"), p4, width = 8, height = 4, dpi = 220,
       device = grDevices::png, type = "cairo")
fwrite(ch, file.path(OUT, "tab_qtn_chr.tsv"), sep = "\t")

## ---------------- 5. locus zoom ----------------
run_plink <- function(chr, s, e, pref) {
  for (g in c(0.1, 0.25, 0.5)) {
    cmd <- sprintf("%s --bfile %s --chr %s --from-bp %d --to-bp %d --maf %g --max-alleles 2 --snps-only --geno %g --make-just-bim --export A --out %s",
                   plink2, gt_bed, chr, max(0, s), e, MAF_MIN, g, pref)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (file.exists(paste0(pref, ".raw"))) return(TRUE)
  }
  FALSE
}
pref <- file.path(data_dir, paste0("tmp_", prefix, "_", TR, "_peak"))
if (!run_plink(ACHR, WIN1, WIN2, pref)) stop("no SNPs in peak window after filters")
raw <- fread(paste0(pref, ".raw"))
setnames(raw, 7:ncol(raw), sub("_[ACGT]+$", "", names(raw)[7:ncol(raw)]))
bim <- fread(paste0(pref, ".bim")); setnames(bim, c("chr", "snp", "cm", "pos", "a1", "a2"))
raw[, IID := as.character(IID)]
snpcols5 <- bim$snp
G <- as.matrix(raw[, ..snpcols5]); cm <- colMeans(G, na.rm = TRUE); cm[is.nan(cm)] <- 0
for (j in seq_len(ncol(G))) { na <- is.na(G[, j]); if (any(na)) G[na, j] <- cm[j] }
lz0 <- merge(data.table(snp = bim$snp, pos = bim$pos), mid[, .(snp, p, neglog)], by = "snp")
LEAD <- lz0[which.max(neglog), snp]
li <- which(bim$snp == LEAD)
r2 <- as.vector(cor(G, G[, li])^2)
lz <- data.table(snp = bim$snp, pos = bim$pos, r2 = r2)
lz <- merge(lz, mid[, .(snp, p, neglog)], by = "snp")
lz[is.na(r2), r2 := 0]
lz[, r2bin := cut(r2, c(-Inf, .2, .4, .6, .8, 1), labels = c("<0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "≥0.8"))]

gl <- genes[chr == ACHR & end > WIN1 & start < WIN2]
gl[, midp := (start + end) / 2]
gl[, y := as.numeric(factor(gene_id))]
gq <- qt[gene_id %in% gl$gene_id][order(-lod)][, .SD[1], by = gene_id][order(-lod)][1:4]
lab_genes <- gq$gene_id
gl[, hl := gene_id %in% lab_genes]

p5a <- ggplot(lz, aes(pos / 1e6, neglog)) +
  geom_point(aes(colour = r2bin), size = 1.25, alpha = 0.9, stroke = 0) +
  geom_vline(xintercept = lz[which.max(neglog), pos] / 1e6, colour = "#C00000", linetype = 3, linewidth = 0.4) +
  scale_colour_manual(values = c("<0.2" = "#DEEBF7", "0.2–0.4" = "#9ECAE1", "0.4–0.6" = "#4292C6",
                                 "0.6–0.8" = "#2171B5", "≥0.8" = "#08306B"),
                      name = expression(与~lead~SNP~的~italic(r)^2), drop = TRUE) +
  scale_x_continuous(limits = c(WIN1/1e6, WIN2/1e6), expand = c(0, 0)) +
  labs(subtitle = sprintf("区域最强单标记：%s（p = %.2g，−log₁₀P = %.1f）", LEAD,
                          lz[which.max(neglog), p], max(lz$neglog)),
       x = NULL, y = expression(-log[10](italic(p)))) +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right", legend.key.size = unit(0.8, "lines"),
        plot.subtitle = element_text(size = 8.4))

p5b <- ggplot(gl) +
  geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = y, yend = y, colour = hl), linewidth = 1.6) +
  geom_text_repel(data = gl[hl == TRUE], aes(midp / 1e6, y, label = gene_id), size = 2.4,
                  colour = "#7F0000", nudge_y = 4, direction = "y",
                  min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#C00000"), guide = "none") +
  scale_x_continuous(limits = c(WIN1/1e6, WIN2/1e6), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = c(4, 12))) +
  labs(x = sprintf("%s 位置 (Mb)", ACHR), y = "基因") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

p5 <- p5a / p5b + plot_layout(heights = c(4, 1.15)) +
  plot_annotation(title = sprintf("%s 主峰区域关联图（%s %.2f–%.2f Mb）", CN, ACHR, WIN1/1e6, WIN2/1e6),
                  theme = theme(plot.title = element_text(family = cjk_font, hjust = .5,
                                                          face = "bold", size = 12)))
ggsave(file.path(FIG, "f5_locus_zoom.png"), p5, width = 9.6, height = 6.4, dpi = 220,
       device = grDevices::png, type = "cairo")

## ---------------- 6. regional QTN signals (50 kb pruning) ----------------
qs <- qtl[chr == ACHR & pos >= WIN1 & pos <= WIN2 & !is.na(lod)][order(-lod)]
sel <- integer(0)
for (i in seq_len(nrow(qs))) {
  if (!length(sel) || all(abs(qs$pos[i] - qs$pos[sel]) >= 50000)) sel <- c(sel, i)
}
sigt <- qs[sel][order(-lod)]
sigt[, significance := factor(significance, levels = c("SIG", "SUG"))]
sigt[, rank := seq_len(.N)]
sigt[, dist_kb := round((pos - APOS) / 1000, 1)]
sigt[, lab := fifelse(!is.na(gene_id) & gene_id != "intergenic" & gene_id != "", gene_id,
                      sprintf("%s %.2fMb", sub("^chr", "", as.character(chr)), pos/1e6))]
fwrite(sigt, file.path(OUT, "tab_signals.tsv"), sep = "\t")

if (nrow(sigt)) {
  p6 <- ggplot(sigt, aes(pos / 1e6, lod)) +
    geom_segment(aes(xend = pos / 1e6, y = 0, yend = lod, colour = significance), linewidth = 0.8) +
    geom_point(aes(colour = significance), size = 2.2) +
    geom_hline(yintercept = 7.7, colour = "grey40", linetype = 2, linewidth = 0.45) +
    geom_text_repel(aes(label = lab), size = 2.4, colour = "#7F0000",
                    nudge_y = 8, direction = "y", min.segment.length = 0,
                    segment.size = 0.2, max.overlaps = Inf, seed = 1) +
    scale_colour_manual(values = c(SIG = "#B2182B", SUG = "#FDAE61"), name = "QTN 显著性", drop = FALSE) +
    scale_x_continuous(limits = c(WIN1/1e6, WIN2/1e6), expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = sprintf("%s 主峰窗口内 QTN 的 LOD 分布", CN),
         subtitle = sprintf("按 50 kb 去冗余后共 %d 个信号；虚线 = LOD 7.7（Bonferroni）", nrow(sigt)),
         x = sprintf("%s 位置 (Mb)", ACHR), y = "LOD") +
    theme_bw(base_size = 11) +
    theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
          legend.position = "top", plot.title = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5, size = 8.6))
  ggsave(file.path(FIG, "f6_qtn_region.png"), p6, width = 9.2, height = 4.2, dpi = 220,
         device = grDevices::png, type = "cairo")
}

## ---------------- 7. haplotype (top genes in peak window) ----------------
HF <- 0.05; K <- 4
gwin <- merge(gl[, .(gene_id, chr, start, end)], qt[, .(n_qtn = .N, max_lod = max(lod)), by = gene_id],
              by = "gene_id", all.x = TRUE)
gwin[is.na(n_qtn), `:=`(n_qtn = 0L, max_lod = NA_real_)]
gwin[, dist := pmin(abs(start - APOS), abs(end - APOS))]
setorder(gwin, -n_qtn, -max_lod, dist)
cand <- gwin[1:min(3, nrow(gwin))]
cat("candidate genes:", paste(cand$gene_id, collapse = ", "), "\n")

phe2 <- copy(phe); setnames(phe2, 1, "IID"); phe2[, IID := as.character(IID)]
assoc <- function(Gm, yv, Cm0) {
  ok <- !is.na(yv); Gm <- Gm[ok, , drop = FALSE]; yv <- yv[ok]
  Cm <- cbind(1, Cm0[ok, , drop = FALSE])
  ry <- qr.resid(qr(Cm), yv); Gc <- qr.resid(qr(Cm), Gm)
  df <- length(yv) - ncol(Cm) - 1
  ssx <- colSums(Gc^2); beta <- as.vector(crossprod(Gc, ry)) / ssx
  s2 <- (sum(ry^2) - beta^2 * ssx) / df
  list(beta = beta, p = 2 * pt(-abs(beta / sqrt(s2 / ssx)), df))
}

hap_res <- list(); hap_freq <- list(); plots <- list()
for (gid in head(gwin$gene_id, 10)) {
  if (length(plots) >= 3) break
  g <- genes[gene_id == gid]
  pref2 <- file.path(data_dir, paste0("tmp_", prefix, "_", TR, "_hap_", gid))
  if (!run_plink(g$chr, g$start - 1000, g$end + 1000, pref2)) next
  hr <- fread(paste0(pref2, ".raw"))
  if (ncol(hr) < 8) { file.remove(paste0(pref2, c(".raw", ".bim", ".log"))); next }
  setnames(hr, 7:ncol(hr), sub("_[ACGT]+$", "", names(hr)[7:ncol(hr)]))
  hb <- fread(paste0(pref2, ".bim")); setnames(hb, c("chr", "snp", "cm", "pos", "a1", "a2"))
  hr[, IID := as.character(IID)]
  mm <- merge(hr, phe2, by = "IID"); mm <- merge(mm, pca, by = "IID")
  hsnps <- hb$snp[hb$snp %in% names(mm)]
  Gh <- as.matrix(mm[, ..hsnps])
  cmh <- round(colMeans(Gh, na.rm = TRUE)); cmh[is.nan(cmh)] <- 0
  for (j in seq_len(ncol(Gh))) { na <- is.na(Gh[, j]); if (any(na)) Gh[na, j] <- cmh[j] }
  Gh <- round(Gh)
  Ch <- as.matrix(mm[, .(pc_1, pc_2, pc_3)]); yh <- mm[[TR]]
  ph <- assoc(Gh, yh, Ch)$p
  cs <- which(ph < 0.05); if (length(cs) < 2) cs <- order(ph)[1:min(2, length(ph))]
  selh <- integer(0)
  for (j in cs[order(ph[cs])]) {
    if (length(selh) >= K) break
    if (length(selh) && any(cor(Gh[, selh, drop = FALSE], Gh[, j])^2 > 0.8)) next
    selh <- c(selh, j)
  }
  if (length(selh) < 2) { file.remove(paste0(pref2, c(".raw", ".bim", ".log"))); next }
  hp <- do.call(paste0, lapply(selh, function(j) as.integer(Gh[, j])))
  fr <- table(hp) / length(hp); keep <- names(fr)[fr >= HF]
  hp2 <- ifelse(hp %in% keep, hp, "其他")
  hap_freq[[length(hap_freq) + 1]] <- data.table(gene_id = gid, haplotype = names(fr),
                                                 freq = as.numeric(fr),
                                                 tag_snps = paste(hsnps[selh], collapse = ";"))
  dt <- data.table(hap = factor(hp2, levels = c(sort(unique(hp2[hp2 != "其他"])), "其他")), y = yh)
  kw <- tryCatch(kruskal.test(y ~ hap, data = dt)$p.value, error = function(e) NA_real_)
  av <- tryCatch(summary(aov(y ~ hap, data = dt))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
  st <- dt[, .(n = .N, mean = mean(y)), by = hap]
  hap_res[[length(hap_res) + 1]] <- data.table(gene_id = gid, n_tags = length(selh),
                                               tag_snps = paste(hsnps[selh], collapse = ";"),
                                               n_hap = length(unique(hp2)), KW_p = kw, ANOVA_p = av,
                                               top_hap = as.character(st[which.max(mean), hap]),
                                               top_mean = max(st$mean))
  ord <- c(as.character(st[hap != "其他"][order(mean), hap]), "其他")
  dt[, hap := factor(as.character(hap), levels = ord)]
  dt[, grp := factor(ifelse(hap == "其他", "低频单倍型（<5%）", "高频单倍型（≥5%）"),
                     levels = c("高频单倍型（≥5%）", "低频单倍型（<5%）"))]
  nm <- dt[, .N, by = hap]; setkey(nm, hap)
  lab <- setNames(sprintf("%s\n(n=%d)", nm$hap, nm$N), nm$hap)
  plots[[gid]] <- ggplot(dt, aes(hap, y, fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.68, linewidth = 0.35) +
    geom_jitter(width = 0.14, size = 0.6, alpha = 0.45, colour = "grey20") +
    scale_x_discrete(labels = lab) +
    scale_fill_manual(values = c("高频单倍型（≥5%）" = "#2F5597", "低频单倍型（<5%）" = "grey78"),
                      name = "单倍型类别") +
    labs(title = gid, subtitle = sprintf("%d 个标签 SNP；KW p = %.2g", length(selh), kw),
         x = "单倍型（标签 SNP 基因型组合）", y = CN) +
    theme_bw(base_size = 9) +
    theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 8), axis.text.x = element_text(size = 7),
          legend.position = "right")
  file.remove(paste0(pref2, c(".raw", ".bim", ".log")))
}
if (length(plots)) {
  p7 <- wrap_plots(plots, nrow = 1, guides = "collect") + plot_layout(guides = "collect") +
    plot_annotation(title = sprintf("%s 主峰候选基因单倍型关联", CN),
                    theme = theme(plot.title = element_text(family = cjk_font, hjust = .5,
                                                            face = "bold", size = 12)))
  ggsave(file.path(FIG, "f7_haplotype.png"), p7, width = 10.5, height = 4.2, dpi = 220,
         device = grDevices::png, type = "cairo")
}
HAP <- rbindlist(hap_res); fwrite(HAP, file.path(OUT, "tab_hap.tsv"), sep = "\t")
fwrite(rbindlist(hap_freq), file.path(OUT, "tab_hap_freq.tsv"), sep = "\t")

## ---------------- tables ----------------
qt2 <- merge(qt, anno, by.x = "gene_id", by.y = "gene", all.x = TRUE)
setorder(qt2, -lod)
fwrite(qt2, file.path(OUT, "tab_qtn_all.tsv"), sep = "\t")
fwrite(head(qt2, 20), file.path(OUT, "tab_qtn_top20.tsv"), sep = "\t")

near <- merge(gl[, .(gene_id, chr, start, end, strand)], anno, by.x = "gene_id", by.y = "gene", all.x = TRUE)
near <- merge(near, qt[, .(n_qtn = .N, max_lod = max(lod)), by = gene_id], by = "gene_id", all.x = TRUE)
near[is.na(n_qtn), `:=`(n_qtn = 0L)]
near[, dist_kb := fifelse(APOS >= start & APOS <= end, 0, pmin(abs(start - APOS), abs(end - APOS)) / 1000)]
near[, at_id := sub("^([A-Za-z0-9]+)\\..*$", "\\1", A.thaliana)]
near[, at_desc := fifelse(grepl('description:"', A.thaliana),
                          sub('.*description:"([^"]*)".*', "\\1", A.thaliana), NA_character_)]
near[, fun := fifelse(!is.na(Preferred_name) & Preferred_name != "-", Preferred_name,
              fifelse(!is.na(`eggNOG free text description`) & `eggNOG free text description` != "-",
                      `eggNOG free text description`, "—"))]
setorder(near, start)
fwrite(near, file.path(OUT, "tab_genes_near_peak.tsv"), sep = "\t")

kcs <- merge(near[gene_id %in% cand$gene_id, .(gene_id, chr, start, end, strand)],
             anno, by.x = "gene_id", by.y = "gene", all.x = TRUE)
kcs <- merge(kcs, qt[, .(n_qtn = .N, max_lod = max(lod)), by = gene_id], by = "gene_id", all.x = TRUE)
fwrite(kcs, file.path(OUT, "tab_candidate_genes.tsv"), sep = "\t")

summ <- data.table(trait = TR, anchor = AMARK, anchor_chr = ACHR, anchor_pos = APOS,
                   window_start = WIN1, window_end = WIN2, n_genes = nrow(near),
                   n_snp = nrow(mid), n_qtn = nrow(qt),
                   n_sig = sum(qt$significance == "SIG", na.rm = TRUE),
                   n_sug = sum(qt$significance == "SUG", na.rm = TRUE),
                   max_lod = max(qt$lod, na.rm = TRUE), min_p = min(mid$p), lambda = lam,
                   n_snp_sig = sum(mid$p < BONF), n_snp_sug = sum(mid$p < 1e-5),
                   n_signals = nrow(sigt), n_haplotypes = nrow(HAP))
fwrite(summ, file.path(OUT, "tab_summary.tsv"), sep = "\t")

file.remove(paste0(pref, c(".raw", ".bim", ".log")))
cat(sprintf("DONE %s %s | genes %d | signals %d | hap %d\n", prefix, TR, nrow(near), nrow(sigt), nrow(HAP)))
