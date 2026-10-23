#!/usr/bin/env Rscript
# DESeq2 差异表达 + volcano / MA / PCA / heatmap
# 用法: Rscript run_deseq2.R <counts_matrix_filtered.tsv> <out_dir> [baseline] [padj_cutoff] [lfc_cutoff]
#
#   counts_matrix_filtered.tsv : 第一列 GeneID, 其余列为样本 (制表符分隔)
#   baseline                   : 基准水平（默认 root1）
#
# 样本名约定: <condition>_rep<N>  →  condition = gsub("_rep\\d+$", "", sample)
#
# 注意: lfcShrink(type="ashr") 可能失败 -> tryCatch 回退到未收缩 LFC (padj/pvalue 不变)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript run_deseq2.R <counts.tsv> <out_dir> [baseline] [padj] [lfc]")
counts_file <- args[1]
out_dir <- args[2]
baseline <- if (length(args) >= 3) args[3] else "root1"
padj_cut <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
lfc_cut  <- if (length(args) >= 5) as.numeric(args[5]) else 1
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts <- read.table(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
samples <- colnames(counts)
condition <- gsub("_rep\\d+$", "", samples)
cat("conditions:", paste(unique(condition), collapse = ", "), "\n")

coldata <- data.frame(row.names = samples, condition = factor(condition))
if (!baseline %in% levels(coldata$condition)) {
  stop("baseline '", baseline, "' not among conditions: ",
       paste(unique(condition), collapse = ", "))
}
coldata$condition <- relevel(coldata$condition, ref = baseline)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)

write.csv(data.frame(sample = samples, condition = condition),
          file.path(out_dir, "sample_info.tsv"), sep = "\t", row.names = FALSE)

plot_volcano <- function(res, contrast_label) {
  df <- data.frame(
    gene = rownames(res),
    lfc = res$log2FoldChange,
    padj = res$padj,
    baseMean = res$baseMean
  )
  df$significant <- !is.na(df$padj) & df$padj < padj_cut &
    abs(df$lfc) >= lfc_cut
  df$lfc[is.na(df$lfc)] <- 0
  df$neglog <- -log10(pmax(df$padj, .Machine$double.xmin))
  df$neglog[!is.finite(df$neglog)] <- max(df$neglog[is.finite(df$neglog)]) + 1
  df <- df[order(df$padj), ]

  n_up <- sum(df$significant & df$lfc > 0)
  n_dn <- sum(df$significant & df$lfc < 0)
  p <- ggplot(df, aes(x = lfc, y = neglog, color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey70", "#D73027")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
    labs(title = paste0("Volcano: ", contrast_label),
         x = "log2 fold change", y = "-log10 adjusted p-value",
         subtitle = paste0("Up: ", n_up, "  Down: ", n_dn)) +
    theme_bw()
  ggsave(file.path(out_dir, paste0("volcano_", contrast_label, ".png")), p,
         width = 7, height = 6, dpi = 600)
  ggsave(file.path(out_dir, paste0("volcano_", contrast_label, ".pdf")), p,
         width = 7, height = 6, device = cairo_pdf)
  invisible(df)
}

res_list <- list()
for (lvl in setdiff(levels(coldata$condition), baseline)) {
  ct <- c("condition", lvl, baseline)
  res <- results(dds, contrast = ct)
  resLFC <- tryCatch(
    lfcShrink(dds, contrast = ct, type = "ashr"),
    error = function(e) { message("ashr unavailable, using unshrunken LFC"); res }
  )
  out <- data.frame(gene = rownames(resLFC), resLFC)
  write.csv(out, file.path(out_dir, paste0("DEG_", lvl, "_vs_", baseline, ".csv")),
            row.names = FALSE)
  res_list[[lvl]] <- res
  plot_volcano(resLFC, paste0(lvl, "_vs_", baseline))
}

# summary
sig_summary <- do.call(rbind, lapply(names(res_list), function(lvl) {
  r <- res_list[[lvl]]
  data.frame(comparison = paste0(lvl, "_vs_", baseline),
             total_tested = sum(!is.na(r$padj)),
             significant = sum(!is.na(r$padj) & r$padj < padj_cut),
             up = sum(!is.na(r$padj) & r$padj < padj_cut & r$log2FoldChange > 0),
             down = sum(!is.na(r$padj) & r$padj < padj_cut & r$log2FoldChange < 0))
}))
write.csv(sig_summary, file.path(out_dir, "DEG_summary.csv"), row.names = FALSE)
print(sig_summary)

# MA plot
r0 <- res_list[[1]]
ma <- data.frame(baseMean = r0$baseMean, lfc = r0$log2FoldChange,
                 padj = r0$padj)
ma$significant <- !is.na(ma$padj) & ma$padj < padj_cut
pma <- ggplot(ma, aes(x = baseMean, y = lfc, color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey70", "#D73027")) +
  scale_x_log10() + labs(title = "MA plot", x = "mean normalized count",
                         y = "log2 fold change") + theme_bw()
ggsave(file.path(out_dir, "MA_first_comparison.png"), pma, width = 7, height = 5, dpi = 600)
ggsave(file.path(out_dir, "MA_first_comparison.pdf"), pma, width = 7, height = 5, device = cairo_pdf)

# PCA
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)), percent.usage = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1] * 100, PC2 = pca$x[, 2] * 100,
                     sample = samples, condition = condition)
ppca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of all samples",
       x = sprintf("PC1 (%.1f%%)", sum(pca$sdev[1]^2) / sum(pca$sdev^2) * 100),
       y = sprintf("PC2 (%.1f%%)", sum(pca$sdev[2]^2) / sum(pca$sdev^2) * 100)) +
  theme_bw()
ggsave(file.path(out_dir, "PCA_all_samples.png"), ppca, width = 7, height = 6, dpi = 600)
ggsave(file.path(out_dir, "PCA_all_samples.pdf"), ppca, width = 7, height = 6, device = cairo_pdf)

# heatmap of top DEGs (union across comparisons)
top_genes <- unique(unlist(lapply(res_list, function(r) {
  if (all(is.na(r$padj))) return(character(0))
  head(rownames(r)[order(r$padj)], 50)
})))
mat <- counts[intersect(top_genes, rownames(counts)), , drop = FALSE]
logmat <- log2(mat + 1)
ann <- data.frame(condition = condition)
rownames(ann) <- samples
ph <- pheatmap(logmat, scale_row = "row", annotation_col = ann,
               color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
               fontsize_row = 6, fontsize_col = 7,
               filename = file.path(out_dir, "heatmap_top_DEGs.pdf"))
png(file.path(out_dir, "heatmap_top_DEGs.png"), width = 10, height = 8, res = 600)
pheatmap(logmat, scale_row = "row", annotation_col = ann,
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         fontsize_row = 6, fontsize_col = 7)
dev.off()

cat("done ->", out_dir, "\n")
