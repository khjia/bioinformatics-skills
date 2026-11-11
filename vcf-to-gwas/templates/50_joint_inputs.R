## 50_joint_inputs.R ---- cross-trait inputs for the joint report (22 traits)
suppressMessages({library(data.table); library(ggplot2); library(reshape2)})
WDIR <- Sys.getenv("GWAS_WDIR", "/media/nfs1/hermes/project/yuanbaofeng/09.postgwas2")
source(file.path(WDIR, "scripts/00_config.R"))

OUT <- file.path(res_dir, "joint"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(fig_dir, "joint");      dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## ---- trait groups (order matters: FA by chain length, then quality, then VE) ----
fa_order <- c("fa_C16_0","fa_C16_1","fa_C18_0","fa_C18_1","fa_C18_2","fa_C18_3n6","fa_C18_3n3",
              "fa_C20_0","fa_C20_1","fa_C20_2","fa_C22_0","fa_C22_1","fa_C24_0","fa_C24_1")
qual     <- c("oil_g100g","protein_g100g","seed_weight_g")
vee_tr   <- c("vee_alpha","vee_beta","vee_gamma","vee_delta","vee_total")
all_tr   <- c(fa_order, qual, vee_tr)
grp <- c(setNames(rep("脂肪酸", length(fa_order)), fa_order),
         setNames(rep("品质", length(qual)), qual),
         setNames(rep("维生素E", length(vee_tr)), vee_tr))
cn_all <- c(trait_cn, vee_cn)
lab_tr <- function(x) unname(cn_all[x])
cat(sprintf("traits: %d (FA %d / quality %d / VE %d)\n", length(all_tr), length(fa_order), length(qual), length(vee_tr)))

## ---- phenotypes ----
p1 <- fread(trait_file);     setnames(p1, 1, "IID")
p2 <- fread(vee_trait_file); setnames(p2, 1, "IID")
phe <- merge(p1, p2, by = "IID")
stopifnot(nrow(phe) == 200, all(all_tr %in% names(phe)))

M <- as.matrix(phe[, ..all_tr])
R <- cor(M, use = "pairwise.complete.obs")
P <- outer(seq_along(all_tr), seq_along(all_tr), Vectorize(function(i, j) {
  if (i == j) return(NA_real_); cor.test(M[, i], M[, j])$p.value }))
cors <- data.table(trait1 = rep(all_tr, each = length(all_tr)),
                   trait2 = rep(all_tr, times = length(all_tr)),
                   r = as.vector(R), p = as.vector(P))
cors <- cors[trait1 != trait2]
fwrite(cors, file.path(OUT, "tab_cor.tsv"), sep = "\t")
fwrite(data.table(trait = all_tr, as.data.table(R)), file.path(OUT, "tab_cor_matrix.tsv"), sep = "\t")

## ---- QTL table (22 traits) ----
q1 <- fread(qtl_file)[, .(trait, marker, chr, pos, lod, significance, gene_id, gene_name)]
q2 <- fread(vee_qtl_file)[, .(trait, marker, chr, pos, lod, significance, gene_id, gene_name)]
q  <- rbind(q1, q2)[trait %in% all_tr]
fwrite(q, file.path(OUT, "tab_qtl_all.tsv"), sep = "\t")
cat(sprintf("QTL rows: %d\n", nrow(q)))

qs <- q[, .(n_qtn = .N, n_sig = sum(significance == "SIG"), n_sug = sum(significance == "SUG"),
            max_lod = max(lod)), by = trait]
qs[, `:=`(group = grp[trait], CN = lab_tr(trait))]
topc <- q[, .SD[which.max(lod)], by = trait][, .(trait, top_chr = chr)]
qs <- merge(qs, topc, by = "trait")
setorder(qs, -n_qtn)
fwrite(qs, file.path(OUT, "tab_qtl_summary.tsv"), sep = "\t")

## ---- shared genes / markers ----
ann <- fread(anno_file, select = c("gene", "Preferred_name", "EC", "Pfam", "A.thaliana",
                                   "eggNOG free text description"))
qg <- q[!is.na(gene_id) & gene_id != "" & gene_id != "intergenic"]
sg <- qg[, .(n_traits = uniqueN(trait), traits = paste(sort(unique(trait)), collapse = ","),
             n_qtn = .N, max_lod = max(lod)), by = .(gene_id)]
sg <- merge(sg, ann, by.x = "gene_id", by.y = "gene", all.x = TRUE)
setorder(sg, -n_traits, -max_lod)
fwrite(sg, file.path(OUT, "tab_shared_genes.tsv"), sep = "\t")

sm <- q[, .(n_traits = uniqueN(trait), traits = paste(sort(unique(trait)), collapse = ","),
            max_lod = max(lod)), by = .(chr, pos)]
setorder(sm, -n_traits, -max_lod)
fwrite(sm, file.path(OUT, "tab_shared_markers.tsv"), sep = "\t")
cat(sprintf("shared genes (>=2 traits): %d ; shared loci (>=2): %d\n",
            sg[n_traits >= 2, .N], sm[n_traits >= 2, .N]))

## gene x trait presence matrix (genes seen in >=3 traits)
gt <- qg[, .(present = 1L), by = .(gene_id, trait)]
topg <- sg[n_traits >= 3, gene_id]
gm <- dcast(gt[gene_id %in% topg], gene_id ~ trait, value.var = "present", fill = 0)
fwrite(gm, file.path(OUT, "tab_gene_trait_matrix.tsv"), sep = "\t")

## ---- hotspots: 1 Mb bins with >= 2 traits ----
q[, bin := floor(pos / 1e6) * 1e6]
hs <- q[, .(n_traits = uniqueN(trait), n_qtn = .N, traits = paste(sort(unique(trait)), collapse = ","),
            max_lod = max(lod)), by = .(chr, bin)][n_traits >= 2]
setorder(hs, -n_traits, -max_lod)
fwrite(hs, file.path(OUT, "tab_hotspots.tsv"), sep = "\t")
cat(sprintf("hotspots (>=2 traits): %d\n", nrow(hs)))

## ---- trait x trait shared-gene counts ----
tt <- qg[, .(present = 1L), by = .(gene_id, trait)]
pair <- tt[tt, on = "gene_id", allow.cartesian = TRUE, nomatch = 0]
pair <- pair[trait != i.trait, .(n_gene = uniqueN(gene_id)), by = .(t1 = trait, t2 = i.trait)]
pair <- pair[t1 %in% all_tr & t2 %in% all_tr]
fwrite(pair, file.path(OUT, "tab_trait_pair_shared.tsv"), sep = "\t")

## ---- KEGG enrichment of shared genes (hypergeometric) ----
kg <- tryCatch(fread(file.path(proj, "07.post_gwas/01.enrichment/term2gene/kegg_term2gene.tsv"),
                     header = FALSE, col.names = c("term", "gene")), error = function(e) NULL)
kn <- tryCatch(fread(file.path(proj, "07.post_gwas/01.enrichment/term2name/kegg_term2name.tsv"),
                     header = FALSE, col.names = c("term", "name")), error = function(e) NULL)
kegg <- NULL
if (!is.null(kg) && nrow(kg)) {
  uni <- unique(kg$gene); setg <- unique(sg[n_traits >= 3, gene_id])
  N <- length(uni); K <- length(intersect(setg, uni)); n <- length(setg)
  kegg <- kg[, .(hit = list(intersect(gene, setg))), by = term]
  kegg[, k := lengths(hit)]
  kegg <- kegg[k >= 2]
  kegg[, p := phyper(k - 1, K, N - K, n, lower.tail = FALSE)]
  if (!is.null(kn)) kegg <- merge(kegg, kn, by = "term", all.x = TRUE)
  kegg[, FDR := p.adjust(p, "BH")]
  setorder(kegg, p)
  kegg[, genes := vapply(hit, function(x) paste(head(x, 6), collapse = ","), character(1))]
  kegg[, hit := NULL]
  fwrite(kegg, file.path(OUT, "tab_kegg_shared.tsv"), sep = "\t")
}

## ================= figures =================
grp_col <- c("脂肪酸" = "#2F5597", "品质" = "#E69F00", "维生素E" = "#009E73")

## f1 correlation heatmap (ordered by group / chain length)
ord <- all_tr
cd <- cors[, .(t1 = factor(trait1, levels = ord), t2 = factor(trait2, levels = ord), r)]
p1f <- ggplot(cd, aes(t2, t1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Pearson r") +
  scale_y_discrete(labels = lab_tr(ord)) +
  scale_x_discrete(labels = lab_tr(ord), position = "top") +
  labs(x = NULL, y = NULL, title = "22 个性状的表型相关矩阵",
       subtitle = "方块颜色为 Pearson r（蓝负红正），* 表示 p < 0.05") +
  theme_bw(base_size = 9) +
  theme(text = element_text(family = cjk_font), axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank(), plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5, size = 8))
ggsave(file.path(FIG, "j1_cor_heatmap.png"), p1f, width = 8.6, height = 7.6, dpi = 220,
       device = grDevices::png, type = "cairo")

## f2 QTL counts per trait
qs[, trait := factor(trait, levels = rev(all_tr))]
p2f <- ggplot(qs, aes(x = n_qtn, y = trait, fill = group)) +
  geom_col(width = .7) +
  geom_text(aes(label = n_qtn), hjust = -0.2, size = 2.8, family = cjk_font) +
  scale_fill_manual(values = grp_col, name = "性状组") +
  scale_y_discrete(labels = lab_tr(levels(qs$trait))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "QTN 数", y = NULL, title = "各性状检出的 QTN 数") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"))
ggsave(file.path(FIG, "j2_qtl_bar.png"), p2f, width = 7.2, height = 5.4, dpi = 220,
       device = grDevices::png, type = "cairo")

## f3 trait x trait shared-gene heatmap
pp <- pair[, .(t1 = factor(t1, levels = rev(ord)), t2 = factor(t2, levels = ord), n_gene)]
p3f <- ggplot(pp, aes(t2, t1, fill = n_gene)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       name = "共享基因数") +
  scale_x_discrete(labels = lab_tr(ord), position = "top") +
  scale_y_discrete(labels = lab_tr(levels(pp$t1))) +
  labs(x = NULL, y = NULL, title = "两两性状共享的 QTL 基因数",
       subtitle = "对角线为空；只在两个性状都检出同一基因时才计数") +
  theme_bw(base_size = 9) +
  theme(text = element_text(family = cjk_font), axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank(), plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5, size = 8))
ggsave(file.path(FIG, "j3_pair_shared.png"), p3f, width = 8.6, height = 7.6, dpi = 220,
       device = grDevices::png, type = "cairo")

## f4 hotspot map
hs[, grp_top := grp[sub(",.*$", "", traits)]]
p4f <- ggplot(hs, aes(x = bin / 1e6, y = chr)) +
  geom_point(aes(size = n_traits, colour = n_traits), alpha = .8) +
  scale_size_continuous(range = c(1.2, 6), name = "涉及性状数") +
  scale_colour_gradientn(colours = c("#FFEDA0", "#FEB24C", "#F03B20", "#BD0026"),
                         name = "涉及性状数") +
  labs(x = "物理位置 (Mb)", y = NULL, title = "多性状热点区间（1 Mb 分箱，≥ 2 个性状）",
       subtitle = sprintf("共 %d 个热点；点越大、颜色越深表示共用的性状越多", nrow(hs))) +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"),
        plot.subtitle = element_text(hjust = .5, size = 8.5))
ggsave(file.path(FIG, "j4_hotspots.png"), p4f, width = 9.4, height = 5.2, dpi = 220,
       device = grDevices::png, type = "cairo")

## f5 top pleiotropic genes
tg <- head(sg[n_traits >= 2], 20)
tg[, lab := fifelse(is.na(Preferred_name) | Preferred_name %in% c("-", "—"), gene_id, Preferred_name)]
tg[, lab := factor(lab, levels = rev(unique(lab)))]
p5f <- ggplot(tg, aes(x = n_traits, y = lab, fill = n_traits)) +
  geom_col(width = .72) +
  geom_text(aes(label = sprintf("%d 个性状 | LOD %.1f", n_traits, max_lod)), hjust = -0.1, size = 2.7,
            family = cjk_font) +
  scale_fill_gradientn(colours = c("#FEB24C", "#F03B20", "#BD0026"), name = "涉及性状数") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(x = "涉及性状数", y = NULL, title = "多效性最强的 20 个基因") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold"))
ggsave(file.path(FIG, "j5_top_genes.png"), p5f, width = 8.4, height = 5.6, dpi = 220,
       device = grDevices::png, type = "cairo")

## f6 VE component correlations
vsub <- cors[trait1 %in% vee_tr & trait2 %in% vee_tr]
p6f <- ggplot(vsub, aes(factor(trait1, levels = vee_tr), factor(trait2, levels = vee_tr), fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.2, family = cjk_font) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(labels = lab_tr(vee_tr), position = "top") +
  scale_y_discrete(labels = lab_tr(vee_tr)) +
  labs(x = NULL, y = NULL, title = "维生素 E 五个组分之间的相关") +
  theme_bw(base_size = 10) +
  theme(text = element_text(family = cjk_font), axis.text.x = element_text(angle = 30, hjust = 0),
        panel.grid = element_blank(), plot.title = element_text(hjust = .5, face = "bold"))
ggsave(file.path(FIG, "j6_vee_cor.png"), p6f, width = 5.4, height = 4.6, dpi = 220,
       device = grDevices::png, type = "cairo")

cat("DONE 50\n")
