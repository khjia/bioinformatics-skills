## Level-by-level circles. Only TF-TF and TF-gene edges are drawn
## (all gene-gene co-expression edges removed). Grey lines, no arrows.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

root <- "/media/nfs1/hermes/togcn_official"
outdir <- file.path(root, "figs")

tree <- utils::read.delim(file.path(outdir, "tree_deg93.tsv"), stringsAsFactors = FALSE)
gr <- utils::read.csv(file.path(outdir, "groups_deg93.csv"), stringsAsFactors = FALSE)
gb <- utils::read.csv(file.path(outdir, "gene_best_tf_deg93.csv"), stringsAsFactors = FALSE)
tree <- merge(tree, gr[, c("gene", "group")], by = "gene", all.x = TRUE)

labs <- unique(gr$group[order(gr$group_order)])
gorder <- setNames(seq_along(labs) - 1L, labs)
tfs <- tree %>% filter(is_TF == 1) %>% select(gene, group) %>% distinct() %>% arrange(group, gene)
genes <- gr %>% filter(!gene %in% tfs$gene) %>% select(gene, group)

TWO_PI <- 2 * pi
stats <- data.frame(group = labs) %>%
  left_join(tfs %>% count(group, name = "n_tf"), by = "group") %>%
  left_join(genes %>% count(group, name = "n_g"), by = "group") %>%
  mutate(n = n_tf + n_g)
RMAX <- 30; RMIN <- 7
stats$R <- RMIN + (RMAX - RMIN) * sqrt(stats$n / max(stats$n))
stats$cx <- cumsum(c(0, 2 * head(stats$R, -1) + 6)) + stats$R

tf_pos <- list()
for (i in seq_len(nrow(stats))) {
  g <- stats$group[i]; R <- stats$R[i]; cx <- stats$cx[i]
  d <- tfs %>% filter(group == g)
  rt <- if (nrow(d) > 1) 0.36 * R else 0
  if (nrow(d) == 1) { d$x <- cx; d$y <- 0 } else {
    a <- TWO_PI * (seq_len(nrow(d)) - 1) / nrow(d) + pi / 2
    d$x <- cx + rt * cos(a); d$y <- rt * sin(a)
  }
  d$sz <- min(1.5, max(0.10, 0.85 * TWO_PI * max(rt, 0.8) / nrow(d)))
  tf_pos[[g]] <- d
}
tf_all <- bind_rows(tf_pos)
tfcx <- setNames(tf_all$x, tf_all$gene); tfcy <- setNames(tf_all$y, tf_all$gene)

g_pos <- list()
for (i in seq_len(nrow(stats))) {
  g <- stats$group[i]; R <- stats$R[i]; cx <- stats$cx[i]
  d <- genes %>% filter(group == g)
  d <- merge(d, gb[, c("gene", "tf")], by = "gene", all.x = TRUE)
  ang <- atan2(tfcy[d$tf] - 0, tfcx[d$tf] - cx)
  miss <- is.na(ang)
  if (any(miss)) ang[miss] <- -pi / 2 + TWO_PI * (seq_len(sum(miss)) - 1) / max(sum(miss), 1)
  d <- d[order(ang), ]
  a2 <- TWO_PI * (seq_len(nrow(d)) - 1) / max(nrow(d), 1) - pi / 2
  d$x <- cx + R * cos(a2); d$y <- R * sin(a2)
  d$sz <- min(0.80, max(0.05, 0.85 * TWO_PI * R / max(nrow(d), 1)))
  g_pos[[g]] <- d
}
g_all <- bind_rows(g_pos)

posx <- c(tfcx, setNames(g_all$x, g_all$gene))
posy <- c(tfcy, setNames(g_all$y, g_all$gene))
grp  <- c(setNames(tfs$group, tfs$gene), setNames(g_all$group, g_all$gene))
isf  <- setNames(c(rep(TRUE, nrow(tfs)), rep(FALSE, nrow(g_all))), names(posx))

rel <- utils::read.csv(file.path(root, "rundeg_0.93", "Node_relation.csv"),
                       stringsAsFactors = FALSE, strip.white = TRUE)
colnames(rel) <- c("a", "b", "pcc")
rel$a <- trimws(rel$a); rel$b <- trimws(rel$b)
rel <- rel[rel$a %in% names(posx) & rel$b %in% names(posx), ]
rel <- rel[isf[rel$a] | isf[rel$b], ]          # keep TF-TF and TF-gene only
rel$kind <- ifelse(isf[rel$a] & isf[rel$b], "TF-TF", "TF-gene")
## orientation: TF side first
swap <- rel$kind == "TF-gene" & !isf[rel$a]
if (any(swap)) {
  tmp <- rel$a[swap]; rel$a[swap] <- rel$b[swap]; rel$b[swap] <- tmp
}
ord_a <- gorder[grp[rel$a]]; ord_b <- gorder[grp[rel$b]]
rel$dord <- ord_b - ord_a
keep <- ifelse(rel$kind == "TF-TF",
               rel$dord == 0 | abs(rel$dord) == 1,          # same or adjacent layer
               rel$dord == 0 | rel$dord == 1)               # same layer, or TF(k) -> gene(k+1)
cat("drop  :", sum(!keep), " (TF->previous-layer gene:", sum(!keep & rel$dord == -1), ")\n")
rel <- rel[keep, ]
rel$cross <- grp[rel$a] != grp[rel$b]
cat("edges kept:", nrow(rel), " TF-TF:", sum(rel$kind == "TF-TF"),
    " TF-gene:", sum(rel$kind == "TF-gene"), " cross:", sum(rel$cross), "\n")
cat("genes with >=1 TF edge:", length(unique(c(rel$a[!isf[rel$a]], rel$b[!isf[rel$b]]))), "\n")

e <- data.frame(x = posx[rel$a], y = posy[rel$a], x.p = posx[rel$b], y.p = posy[rel$b],
                tf = rel$kind == "TF-TF", cross = rel$cross)

circ <- do.call(rbind, lapply(seq_len(nrow(stats)), function(i) {
  th <- seq(0, TWO_PI, length.out = 240)
  data.frame(group = stats$group[i], x = stats$cx[i] + stats$R[i] * cos(th), y = stats$R[i] * sin(th))
}))
stats$lab <- sprintf("%s\nTF %d / gene %d", stats$group, stats$n_tf, stats$n_g)

p <- ggplot() +
  geom_path(data = circ, aes(x, y, group = group), colour = "grey80", linewidth = 0.3) +
  geom_segment(data = e[!e$tf, ], aes(x = x, y = y, xend = x.p, yend = y.p),
               colour = "grey85", linewidth = 0.06, alpha = 1) +
  geom_segment(data = e[e$tf, ], aes(x = x, y = y, xend = x.p, yend = y.p),
               colour = "grey55", linewidth = 0.09, alpha = 1) +
  geom_point(data = g_all, aes(x, y, size = sz), colour = "#2C5F8A", alpha = 0.9) +
  geom_point(data = tf_all, aes(x, y, size = sz), colour = "#B2182B") +
  scale_size_identity() +
  geom_text(data = stats, aes(x = cx, y = -RMAX - 2.5, label = lab), size = 2.6, lineheight = 0.95, vjust = 1) +
  coord_fixed(clip = "off", xlim = c(-4, max(stats$cx + stats$R) + 4),
              ylim = c(-RMAX - 12, RMAX + 5)) +
  labs(title = "TO-GCN regulatory network, DEG set (r >= 0.93, software default) TF edges only",
       subtitle = paste0("circle area scales with node number (radius ~ sqrt(n));  red = transcription factor (inner ring), ",
                         "blue = target gene (outer rim);  ", format(nrow(e), big.mark = ","),
                         " edges: TF-TF (dark) + TF-gene (light);  gene-gene edges removed")) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 7.2, colour = "grey30", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5))

W <- max(stats$cx + stats$R) + 8
H <- (2 * RMAX + 30)
ggsave(file.path(outdir, "network_tf_edges_deg93.png"), p, width = W, height = H, units = "mm",
       dpi = 500, type = "cairo", bg = "white", limitsize = FALSE)
cat("DONE\n")
