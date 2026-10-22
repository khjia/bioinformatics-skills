# 进化树（上）+ 群体结构（下）合二为一图 —— 矩形布局，自算坐标，不依赖 ggtreeExtra
# 用法: Rscript tree_structure_combined.R <P> <PRE> <SUF> <kmin> <kmax>
#   P   分析目录（下面有 tree/<PRE>.nwk, admixture/<PRE>.fam, admixture/<PRE>.<K>.Q）
#   PRE 剪枝前缀（如 wuhuaguo_195_pruned）; SUF 图名后缀（如 _195）; kmin/kmax 画哪些 K
suppressMessages({library(ape); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
.args <- commandArgs(trailingOnly = TRUE)
P    <- ifelse(length(.args) >= 1, .args[1], "/path/to/05.population")
PRE  <- ifelse(length(.args) >= 2, .args[2], "panel_pruned")
SUF  <- ifelse(length(.args) >= 3, .args[3], "")
KS   <- if (length(.args) >= 5) as.integer(.args[4]):as.integer(.args[5]) else 2:6
FIG  <- file.path(P, "figures"); dir.create(FIG, showWarnings = FALSE)
FONT <- "Times New Roman"; POW <- 0.45
PAL  <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")

tree <- ladderize(read.tree(file.path(P, "tree", paste0(PRE, ".nwk"))))
n <- length(tree$tip.label); tips <- tree$tip.label
ids <- read.table(file.path(P, "admixture", paste0(PRE, ".fam")))$V2
stopifnot(setequal(ids, tips))
NMARK <- NA_integer_
.rawf <- file.path(P, "tree", paste0(PRE, ".raw"))
if (file.exists(.rawf)) NMARK <- length(strsplit(readLines(.rawf, n = 1), "\t")[[1]]) - 6
Q <- lapply(KS, function(k) { m <- as.matrix(read.table(sprintf("%s/admixture/%s.%d.Q", P, PRE, k)))
  rownames(m) <- ids; colnames(m) <- paste0("C", seq_len(ncol(m))); m })
names(Q) <- as.character(KS)

## 矩形布局: x = 叶序号, y = 到根距离（根在上）
x_tip <- setNames(seq_len(n), tips)
depth <- node.depth.edgelength(tree); depth <- depth - min(depth)
Y <- (depth / max(depth))^POW
desc_tips <- function(node) { if (node <= n) return(node)
  unlist(lapply(which(tree$edge[, 1] == node), function(i) desc_tips(tree$edge[i, 2]))) }
x_node <- setNames(numeric(n + tree$Nnode), as.character(seq_len(n + tree$Nnode)))
for (i in seq_len(n)) x_node[i] <- x_tip[[tips[i]]]
for (nd in (n + 1):(n + tree$Nnode)) { d <- desc_tips(nd); x_node[nd] <- (min(d) + max(d)) / 2 }
segs <- list(); nn <- 1
for (e in seq_len(nrow(tree$edge))) {
  p <- tree$edge[e, 1]; ch <- tree$edge[e, 2]
  segs[[nn]] <- data.frame(x = c(x_node[p], x_node[ch]), y = c(Y[p], Y[p]), seg = paste0("h", e)); nn <- nn + 1
  segs[[nn]] <- data.frame(x = c(x_node[ch], x_node[ch]), y = c(Y[p], Y[ch]), seg = paste0("v", e)); nn <- nn + 1
}
tree_df <- bind_rows(segs)

bar_df <- bind_rows(lapply(KS, function(k) {
  m <- Q[[as.character(k)]]; ord <- order(-colMeans(m))
  as.data.frame(m[, ord, drop = FALSE]) %>% tibble::rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "cluster", values_to = "prop") %>%
    mutate(K = factor(paste0("K=", k), levels = paste0("K=", KS)),
           cluster = factor(cluster, levels = paste0("C", ord)), xi = x_tip[sample]) }))

XL <- c(0.5, n + 0.5)
p_tree <- ggplot() +
  geom_path(data = tree_df, aes(x, y, group = seg), linewidth = 0.28, color = "grey25") +
  geom_text(data = data.frame(x = seq_len(n), y = 1.015, lab = tips), aes(x, y, label = lab),
            angle = 90, hjust = 0, vjust = 0.5, family = FONT, size = 1.9, color = "grey15") +
  scale_x_continuous(limits = XL, expand = c(0, 0)) +
  scale_y_reverse(limits = c(1.30, 0), expand = c(0, 0)) +
  labs(subtitle = sprintf("Neighbour-joining tree of %d wild accessions (%s LD-pruned SNPs; vertical axis = branch length^%.2f)",
                          n, format(NMARK, big.mark = ","), POW)) +
  theme_void(base_family = FONT) +
  theme(plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey25", margin = margin(b = 3)),
        plot.margin = margin(6, 12, 0, 12))
p_bars <- ggplot(bar_df, aes(xi, prop, fill = cluster)) +
  geom_col(width = 1, color = NA) + facet_grid(K ~ ., switch = "y") +
  scale_x_continuous(limits = XL, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = "Accessions (same order as the tree)", y = "Ancestry proportion") +
  theme_bw(base_size = 6, base_family = FONT) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, size = 6), strip.background = element_rect(fill = "grey93", color = NA),
        plot.margin = margin(0, 12, 6, 12))
fig <- p_tree / p_bars + plot_layout(heights = c(3.1, 2.6))
ggsave(file.path(FIG, paste0("fig5_tree_structure_K", min(KS), "-", max(KS), SUF, ".pdf")), fig,
       width = 15, height = 8.6, device = cairo_pdf)
ggsave(file.path(FIG, paste0("fig5_tree_structure_K", min(KS), "-", max(KS), SUF, ".png")), fig,
       width = 15, height = 8.6, dpi = 600, bg = "white", limitsize = FALSE)
