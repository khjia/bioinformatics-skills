# 环形 NJ 树 + K=2..N ADMIXTURE 外圈（自算极坐标，不依赖 ggtreeExtra）
# 用法: 改下面 CFG 段路径后 Rscript 本脚本；输入 = tree/*.nwk (ape 可读) + admixture/*.fam + admixture/*.<K>.Q
# 规范: 图内文字仅英文; Times New Roman; PDF=cairo_pdf 内嵌字体; PNG 600 dpi
suppressMessages({library(ape); library(ggplot2); library(dplyr); library(tidyr)})

## ---------------- CFG ----------------
P    <- "/path/to/05.population"          # 项目分析目录
PRE  <- "panel_194_pruned"                # tree/PRE.nwk, admixture/PRE.fam, admixture/PRE.<K>.Q
KS   <- 2:10                              # 画哪些 K（内圈→外圈）
TIT  <- "Population structure of 194 wild accessions"
FONT <- "Times New Roman"
PAL  <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
          "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
W <- 0.075; GAP <- 0.012                  # 环宽 / 环间距（半径单位）
POW  <- 0.45                              # 枝长幂变换（NJ 树内部分支极短，必须压缩）
## -------------------------------------

FIG <- file.path(P, "figures"); dir.create(FIG, showWarnings = FALSE)
tree <- ladderize(read.tree(file.path(P, "tree", paste0(PRE, ".nwk"))))
n <- length(tree$tip.label); tips <- tree$tip.label
ids <- read.table(file.path(P, "admixture", paste0(PRE, ".fam")))$V2   # .Q 行序 = fam 序
stopifnot(setequal(ids, tips))
Q <- lapply(KS, function(k) {
  m <- as.matrix(read.table(sprintf("%s/admixture/%s.%d.Q", P, PRE, k)))
  rownames(m) <- ids; colnames(m) <- paste0("C", seq_len(ncol(m))); m })
names(Q) <- as.character(KS)

## ---- 极坐标布局 ----
theta <- setNames(2 * pi * (seq_len(n) - 1) / n, tips)
depth <- node.depth.edgelength(tree); depth <- depth - min(depth)
R <- (depth / max(depth))^POW
desc_tips <- function(node) {
  if (node <= n) return(node)
  unlist(lapply(which(tree$edge[, 1] == node), function(i) desc_tips(tree$edge[i, 2])))
}
node_theta <- setNames(numeric(n + tree$Nnode), as.character(seq_len(n + tree$Nnode)))
for (i in seq_len(n)) node_theta[i] <- theta[[tips[i]]]
for (nd in (n + 1):(n + tree$Nnode)) {            # 内部节点角度 = 后代叶索引中点 → 无环绕问题
  d <- desc_tips(nd); node_theta[nd] <- 2 * pi * ((min(d) - 1) + (max(d) - 1)) / (2 * n)
}
arc <- function(r, t0, t1, npts = 12) {
  tt <- seq(t0, t1, length.out = npts)
  data.frame(x = r * cos(tt), y = r * sin(tt)) }
segs <- list(); nn <- 1
for (e in seq_len(nrow(tree$edge))) {
  p <- tree$edge[e, 1]; ch <- tree$edge[e, 2]
  tp <- node_theta[p]; tc <- node_theta[ch]
  segs[[nn]] <- data.frame(x = c(R[p] * cos(tp), R[ch] * cos(tp)),
                           y = c(R[p] * sin(tp), R[ch] * sin(tp)),
                           seg = paste0("e", e, "r")); nn <- nn + 1
  if (abs(tc - tp) > 1e-9) { a <- arc(R[ch], tp, tc); a$seg <- paste0("e", e, "a")
                             segs[[nn]] <- a; nn <- nn + 1 }
}
tree_df <- bind_rows(segs)          # geom_path 必须 aes(group = seg)，否则各段连成折线

## ---- admixture 外圈 ----
sector <- function(r0, r1, th, half, npts = 6) {
  ts <- seq(th - half, th + half, length.out = npts)
  data.frame(x = c(r0 * cos(ts), rev(r1 * cos(ts))),
             y = c(r0 * sin(ts), rev(r1 * sin(ts)))) }
R0 <- 1.045; HALF <- (2 * pi / n) * 0.45
rp <- list(); rr <- 1
for (j in seq_along(KS)) {
  k <- KS[j]; m <- Q[[as.character(k)]]
  ri <- R0 + (j - 1) * (W + GAP)
  ord <- order(-colMeans(m))                      # 簇按全局比例降序 → 跨 K 颜色语义一致
  for (i in seq_len(n)) {
    s <- tips[i]; pv <- m[s, ord]; pv <- pv / sum(pv); cum <- c(0, cumsum(pv))
    for (c in seq_len(k)) {
      poly <- sector(ri + cum[c] * W, ri + cum[c + 1] * W, theta[[s]], HALF)
      poly$grp <- paste(k, i, c); poly$col <- PAL[c]; rp[[rr]] <- poly; rr <- rr + 1
    }
  }
}
rings <- bind_rows(rp)
RN <- R0 + length(KS) * (W + GAP); lab_r <- RN + 0.035
phi <- theta[tips] * 180 / pi; flip <- phi > 90 & phi < 270
lab_df <- data.frame(x = lab_r * cos(theta[tips]), y = lab_r * sin(theta[tips]), lab = tips,
                     ang = ifelse(flip, phi - 180, phi), hjust = ifelse(flip, 1, 0))
leg <- data.frame(x = -2.55, y = c(-2.02, -2.24),
                  lab = c("cluster 1 (major)", paste0("cluster 2   (K=", KS[1], " ring)")),
                  col = PAL[1:2])

p <- ggplot() +
  geom_path(data = tree_df, aes(x, y, group = seg), linewidth = 0.3, color = "grey30") +
  geom_polygon(data = rings, aes(x, y, group = grp, fill = col), color = NA) +
  geom_point(data = leg, aes(x, y), shape = 15, size = 2.6, color = leg$col) +
  geom_text(data = leg, aes(x + 0.10, y, label = lab), family = FONT, size = 2.6, hjust = 0) +
  annotate("text", x = -2.55, y = -2.46, family = FONT, size = 2.6, hjust = 0,
           label = paste0("Rings, inner to outer: K = ", paste(KS, collapse = ", "))) +
  geom_text(data = lab_df, aes(x, y, label = lab, angle = ang, hjust = hjust),
            family = FONT, size = 2.0, color = "grey15", vjust = 0.5) +
  scale_fill_identity() +
  coord_equal(xlim = c(-2.75, 2.35), ylim = c(-2.75, 2.35), clip = "off") +
  labs(title = TIT,
       subtitle = sprintf("Neighbour-joining tree based on LD-pruned SNPs (radial axis = branch length^%.2f for legibility)\nRing colour = ancestry cluster ranked by mean proportion within each K (cluster 1 = major)", POW)) +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = 17, hjust = .5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 9.5, hjust = .5, color = "grey25", lineheight = 1.15))
ggsave(file.path(FIG, paste0("fig1_circular_tree_admixture_K", min(KS), "-", max(KS), ".pdf")), p,
       width = 15, height = 15, device = cairo_pdf)
ggsave(file.path(FIG, paste0("fig1_circular_tree_admixture_K", min(KS), "-", max(KS), ".png")), p,
       width = 15, height = 15, dpi = 600, bg = "white", limitsize = FALSE)
cat("done\n")
