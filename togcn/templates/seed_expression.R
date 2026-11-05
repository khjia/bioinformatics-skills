## Seed genes of the pea TO-GCN: expression across the 5 time points
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

root <- "/media/nfs1/hermes/togcn_official"
outdir <- file.path(root, "figs")
counts_file <- "/media/nfs1/hermes/project/wandou/RNA-seq/05.featurecounts/raw_counts_matrix.tsv"

seeds <- readLines(file.path(root, "data/seeds.txt"))
seeds <- seeds[nzchar(seeds)]

ca <- utils::read.delim(counts_file, row.names = 1, check.names = FALSE)
sel <- grep("^(0|30|60|90|120)_[123]$", colnames(ca), value = TRUE)
sel <- sel[order(as.numeric(sub("_.*$", "", sel)), sub("^.*_", "", sel))]
cnt <- as.matrix(ca[, sel]); rm(ca)
time <- as.numeric(sub("_.*$", "", sel)); times <- sort(unique(time))
cpm <- sweep(cnt, 2L, colSums(cnt), "/") * 1e6

d <- data.frame(gene = seeds, cpm[match(seeds, rownames(cpm)), , drop = FALSE], check.names = FALSE)
long <- d %>% pivot_longer(-gene, names_to = "sample", values_to = "cpm") %>%
  mutate(time = as.numeric(sub("_.*$", "", sample)))
summ <- long %>% group_by(gene, time) %>%
  summarise(mean = mean(cpm), sd = sd(cpm), .groups = "drop") %>%
  mutate(logcpm = log2(mean + 1), logcpm_sd = sd / ((mean + 1) * log(2)))

## TF family annotation
fam <- utils::read.delim("/tmp/itak_pea/tf_classification.txt", header = FALSE, stringsAsFactors = FALSE)
fam$gene <- sub("\\..*$", "", fam$V1)
summ$family <- fam$family[match(summ$gene, fam$gene)]

## order genes by the time of their expression maximum
pk <- summ %>% group_by(gene) %>% slice_max(mean, n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(time, desc(mean))
summ$gene <- factor(summ$gene, levels = rev(pk$gene))

## gene-level z-score
z <- summ %>% group_by(gene) %>%
  mutate(z = (mean - mean(mean)) / sd(mean)) %>% ungroup()

p1 <- ggplot(z, aes(time, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       name = "z-score\n(across time)") +
  scale_x_continuous(breaks = times, expand = c(0, 0)) +
  labs(x = "Time after treatment (min)", y = NULL,
       title = "Seed genes: relative expression") +
  theme_bw(base_size = 8.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.5),
        plot.title = element_text(face = "bold", size = 9),
        legend.key.height = unit(8, "mm"))

cols <- viridisLite::viridis(length(seeds))
names(cols) <- seeds
p2 <- ggplot(summ, aes(time, logcpm, colour = gene, fill = gene)) +
  geom_ribbon(aes(ymin = logcpm - logcpm_sd, ymax = logcpm + logcpm_sd), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = cols, name = "seed gene") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(breaks = times) +
  labs(x = "Time after treatment (min)", y = "log2(CPM + 1)",
       title = "Seed genes: absolute expression (mean ± SD, n = 3)") +
  theme_bw(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 5.6), legend.key.height = unit(2.6, "mm"),
        legend.title = element_text(size = 7))

utils::write.csv(summ, file.path(outdir, "seed_expression.csv"), row.names = FALSE)

p <- patchwork::wrap_plots(p1, p2, widths = c(1, 1.55))
ggsave(file.path(outdir, "seed_expression.pdf"), p, width = 250, height = 95, units = "mm",
       device = grDevices::cairo_pdf)
ggsave(file.path(outdir, "seed_expression.png"), p, width = 250, height = 95, units = "mm",
       dpi = 600, type = "cairo", bg = "white")
cat("seeds:", length(seeds), "\nDONE\n")
