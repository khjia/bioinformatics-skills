## Final figures with data-driven merged layers (small-layer + trend merging)
suppressPackageStartupMessages({ library(ggplot2) })

root <- "/media/nfs1/hermes/togcn_official"
outdir <- file.path(root, "figs")
counts_file <- "/media/nfs1/hermes/project/wandou/RNA-seq/05.featurecounts/raw_counts_matrix.tsv"

counts_all <- utils::read.delim(counts_file, row.names = 1, check.names = FALSE)
sel <- grep("^(0|30|60|90|120)_[123]$", colnames(counts_all), value = TRUE)
sel <- sel[order(as.numeric(sub("_.*$", "", sel)), sub("^.*_", "", sel))]
counts <- as.matrix(counts_all[, sel, drop = FALSE]); rm(counts_all)
keep <- rowMeans(counts > 0) >= 0.5 & rowMeans(counts) >= 5
counts <- counts[keep, , drop = FALSE]
time <- as.numeric(sub("_.*$", "", sel)); times <- sort(unique(time))
cpm <- sweep(counts, 2L, colSums(counts), "/") * 1e6
logcpm <- log2(cpm + 1)
mu <- vapply(times, function(t) rowMeans(logcpm[, time == t, drop = FALSE]), numeric(nrow(logcpm)))
colnames(mu) <- times
gz <- t(scale(t(mu))); gz[!is.finite(gz)] <- 0

plot_one <- function(tag, title_tag) {
  gr <- utils::read.csv(file.path(outdir, paste0("groups_", tag, ".csv")), stringsAsFactors = FALSE)
  gr <- gr[gr$gene %in% rownames(mu), ]
  labs <- unique(gr$group[order(gr$group_order)])
  res <- do.call(rbind, lapply(labs, function(g) {
    gs <- gr$gene[gr$group == g]
    data.frame(group = g, n = length(gs), time = times,
               z_mean = colMeans(gz[gs, , drop = FALSE]),
               z_median = apply(gz[gs, , drop = FALSE], 2, median),
               log_mean = colMeans(mu[gs, , drop = FALSE]),
               log_median = apply(mu[gs, , drop = FALSE], 2, median),
               stringsAsFactors = FALSE)
  }))
  pk <- times[apply(mu[gr$gene, , drop = FALSE], 1, which.max)]
  comp <- table(factor(gr$group, levels = labs), factor(pk, levels = times))
  summ <- data.frame(group = labs, n = as.integer(table(factor(gr$group, levels = labs))),
                     peak_frac = times[apply(comp, 1, which.max)],
                     pct_peak = round(100 * apply(comp, 1, max) / rowSums(comp), 1))
  summ$peak_0 <- round(100 * comp[, "0"] / rowSums(comp), 1)
  summ$peak_30 <- round(100 * comp[, "30"] / rowSums(comp), 1)
  summ$peak_60 <- round(100 * comp[, "60"] / rowSums(comp), 1)
  summ$peak_90 <- round(100 * comp[, "90"] / rowSums(comp), 1)
  summ$peak_120 <- round(100 * comp[, "120"] / rowSums(comp), 1)
  utils::write.csv(summ, file.path(outdir, paste0("merged_summary_", tag, ".csv")), row.names = FALSE)
  cat("\n===", tag, "=== peak-time composition (%)\n"); print(summ, row.names = FALSE)

  res$gl <- paste0(res$group, "  (n=", res$n, ")")
  res$gl <- factor(res$gl, levels = rev(unique(res$gl)))
  draw <- function(cols, fillname, title, mode, file) {
    a <- res[, c("gl", "time", cols[1])]; colnames(a)[3] <- "value"; a$stat <- "Mean"
    b <- res[, c("gl", "time", cols[2])]; colnames(b)[3] <- "value"; b$stat <- "Median"
    long <- rbind(a, b)
    long$stat <- factor(long$stat, levels = c("Mean", "Median"))
    p <- ggplot(long, aes(time, gl, fill = value)) +
      geom_tile(colour = "white", linewidth = 1.0) +
      facet_wrap(~stat, nrow = 1) +
      scale_x_continuous(breaks = times, expand = c(0, 0)) +
      labs(x = "Time after treatment (min)", y = NULL, title = title, fill = fillname) +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.4),
            strip.text = element_text(face = "bold", size = 10),
            axis.text = element_text(size = 9, colour = "black"),
            plot.title = element_text(face = "bold", size = 11),
            legend.key.height = unit(12, "mm"),
            legend.title = element_text(size = 9.5))
    if (mode == "div") {
      p <- p + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0)
    } else {
      p <- p + scale_fill_gradientn(colours = c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))
    }
    ggsave(file.path(outdir, paste0(file, ".pdf")), p, width = 180, height = 95, units = "mm", device = grDevices::cairo_pdf)
    ggsave(file.path(outdir, paste0(file, ".png")), p, width = 180, height = 95, units = "mm", dpi = 600, type = "cairo", bg = "white")
  }
  draw(c("z_mean", "z_median"), "gene-level z-score\n(log2 CPM, across time)",
       paste0("TO-GCN merged layers (", title_tag, ")"), "div", paste0("final_z_", tag))
  draw(c("log_mean", "log_median"), "log2(CPM+1)",
       paste0("TO-GCN merged layers (", title_tag, "), absolute expression"), "seq", paste0("final_log_", tag))
}

plot_one("all_raw", "r >= 0.95")
plot_one("all_z", "gene z-score input, r >= 0.95")
cat("\nDONE\n")
