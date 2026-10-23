#!/usr/bin/env Rscript
# 多组 clusterProfiler 富集（x = group 的散点图，用户标准样式）
# 用法: Rscript run_multigroup_enrichment.R <deg_dir> <emapper_annotations> <out_dir>
#
#   deg_dir                : 每个比较一个 DEG_*.csv（列名含 gene/log2FoldChange/padj/pvalue）
#   emapper_annotations    : eggNOG-mapper 的 *.emapper.annotations
#   out_dir                : 输出目录（results/ + figures/）
#
# 关键点:
#   - TERM2GENE 列名必须是大写 TERM 和 GENE，否则 enricher() 静默返回 0 条
#   - KEGG 用 KEGG_Pathway (ko#####)，不是 KEGG_ko (K#####)；排除动物通路前缀
#   - 背景集 = 该比较里 pvalue 非 NA 的全部基因
#   - GO 传 TERM2NAME 显示可读名称，只画 BP 子本体
#   - 绘图统一英文 + Arial + 600 dpi PNG + cairo_pdf

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(GO.db)
  library(AnnotationDbi)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript run_multigroup_enrichment.R <deg_dir> <emapper.annotations> <out_dir>")
deg_dir <- args[1]; ann_file <- args[2]; out_dir <- args[3]
res_dir <- file.path(out_dir, "results"); fig_dir <- file.path(out_dir, "figures")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

PADJ <- 0.05

# ---------- 1. 读 eggNOG 注释 ----------
col_names <- c("query","seed_ortholog","evalue","score","eggNOG_OGs",
               "max_annot_lvl","COG_category","Description","Preferred_name",
               "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
               "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy",
               "BiGG_Reaction","PFAMs")
ann <- read.delim(ann_file, comment.char = "#", header = FALSE, col.names = col_names,
                  sep = "\t", stringsAsFactors = FALSE)
# 用 [.] 而不是 \\. 避免 R 脚本转义问题
ann$gene <- sub("[.]mRNA[0-9]+$", "", ann$query)

go_ann <- ann[ann$GOs != "" & ann$GOs != "-", c("gene", "GOs")]
go_list <- strsplit(go_ann$GOs, ",")
term2gene_go <- data.frame(TERM = unlist(go_list),
                           GENE = rep(go_ann$gene, lengths(go_list)),
                           stringsAsFactors = FALSE)

kp_ann <- ann[ann$KEGG_Pathway != "" & ann$KEGG_Pathway != "-", c("gene", "KEGG_Pathway")]
kp_list <- strsplit(kp_ann$KEGG_Pathway, ",")
kp_terms <- unlist(kp_list); kp_genes <- rep(kp_ann$gene, lengths(kp_list))
keep_ko <- grepl("^ko[0-9]{5}$", kp_terms)     # 只保留 ko#####，去掉 map##### 重复
term2gene_ko <- data.frame(TERM = kp_terms[keep_ko], GENE = kp_genes[keep_ko],
                           stringsAsFactors = FALSE)
# 植物过滤：排除动物特异通路
term2gene_ko <- term2gene_ko[!grepl("^(ko05|ko047|ko049|ko046|ko043)",
                                    term2gene_ko$TERM), ]

# GO 名称（TERM2NAME）
go_ids <- unique(term2gene_go$TERM)
go_name <- AnnotationDbi::select(GO.db, keys = go_ids, columns = c("TERM", "ONTOLOGY"),
                                  keytype = "GOID")
term2name_go <- data.frame(TERM = go_ids, NAME = go_name$TERM[match(go_ids, go_name$GOID)],
                           stringsAsFactors = FALSE)

# KEGG 名称（全局 DB）
kegg_db <- "/media/nfs1/hermes/db/kegg_pathway_names.tsv"
term2name_ko <- if (file.exists(kegg_db)) {
  kn <- read.delim(kegg_db, header = FALSE, stringsAsFactors = FALSE)
  data.frame(TERM = kn[[1]], NAME = kn[[2]], stringsAsFactors = FALSE)
} else NULL

# ---------- 2. 逐组富集 ----------
deg_files <- list.files(deg_dir, pattern = "^DEG_.*\\.csv$", full.names = TRUE)
if (!length(deg_files)) stop("no DEG_*.csv found in ", deg_dir)

run_one <- function(f, db) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  gene_col <- if ("gene" %in% names(d)) "gene" else names(d)[1]
  d[[gene_col]] <- sub("[.][0-9]+$", "", d[[gene_col]])
  sig <- d[!is.na(d$padj) & d$padj < PADJ, gene_col]
  bg  <- d[!is.na(d$pvalue), gene_col]
  if (!length(sig)) return(NULL)
  t2g <- if (db == "GO") term2gene_go else term2gene_ko
  t2n <- if (db == "GO") term2name_go else term2name_ko
  sig <- intersect(sig, t2g$GENE); bg <- intersect(bg, t2g$GENE)
  if (!length(sig)) return(NULL)
  eg <- enricher(gene = sig, universe = bg,
                 TERM2GENE = t2g, TERM2NAME = t2n,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 minGSSize = 5, maxGSSize = 500)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
  df <- as.data.frame(eg)
  df$group <- sub("^DEG_", "", sub("\\.csv$", "", basename(f)))
  df$ID <- as.character(df$ID)
  if (db == "GO") {
    onto <- AnnotationDbi::select(GO.db, keys = df$ID, columns = "ONTOLOGY", keytype = "GOID")
    df$ONTOLOGY <- onto$ONTOLOGY[match(df$ID, onto$GOID)]
    df <- df[!is.na(df$ONTOLOGY) & df$ONTOLOGY == "BP", ]
  }
  if (!nrow(df)) return(NULL)
  df
}

all_go  <- do.call(rbind, lapply(deg_files, run_one, db = "GO"))
all_ko  <- do.call(rbind, lapply(deg_files, run_one, db = "KEGG"))
if (!is.null(all_go))  write.csv(all_go,  file.path(res_dir, "GO_enrichment.csv"),  row.names = FALSE)
if (!is.null(all_ko))  write.csv(all_ko,  file.path(res_dir, "KEGG_enrichment.csv"), row.names = FALSE)

# ---------- 3. 绘图（x = group） ----------
plot_grouped <- function(df, title, fname) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  d <- df
  d$Description <- ifelse(is.na(d$Description) | d$Description == "", as.character(d$ID), d$Description)
  d$neglog <- -log10(pmax(d$p.adjust, .Machine$double.xmin))
  d$neglog[!is.finite(d$neglog)] <- max(d$neglog[is.finite(d$neglog)]) + 1
  # 每组取 top 10，合并去重，按各组最小 padj 排序
  top <- do.call(rbind, lapply(split(d, d$group), function(x) head(x[order(x$p.adjust), ], 10)))
  keep <- unique(top$Description)
  d2 <- d[d$Description %in% keep, ]
  ord <- aggregate(p.adjust ~ Description, data = d2, FUN = min)
  ord <- ord[order(ord$p.adjust), ]
  d2$Description <- factor(d2$Description, levels = rev(ord$Description))
  p <- ggplot(d2, aes(x = group, y = Description, size = Count, color = neglog)) +
    geom_point() +
    scale_color_gradient(low = "#4575B4", high = "#D73027") +
    labs(title = title, x = NULL, y = NULL,
         size = "Gene count", color = "-log10 adjusted p") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(fig_dir, paste0(fname, ".pdf")), p, width = 8, height = 7, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(fname, ".png")), p, width = 8, height = 7, dpi = 600)
  invisible(p)
}

plot_grouped(all_go, "GO biological process enrichment", "GO_BP_multigroup")
plot_grouped(all_ko, "KEGG pathway enrichment", "KEGG_multigroup")
cat("done ->", out_dir, "\n")
