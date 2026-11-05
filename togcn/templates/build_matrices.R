## Build TO-GCN input matrices (TPM) for a time-series RNA-seq experiment.
## Filtering rule: keep a gene if its mean TPM in ANY group (time point) is > 1.
## Output format: geneID \t sample1 ... sampleN   (no header, no row names)

counts_file <- "raw_counts_matrix.tsv"
len_file    <- "gene_length.tsv"        # gene_id <tab> length (exon union from GTF)
tf_file     <- "tf_classification.txt"  # iTAK output
outdir      <- "togcn_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- counts ----
counts_all <- utils::read.delim(counts_file, row.names = 1, check.names = FALSE)
sel <- grep("^(0|30|60|90|120)_[123]$", colnames(counts_all), value = TRUE)
sel <- sel[order(as.numeric(sub("_.*$", "", sel)), sub("^.*_", "", sel))]  # keep time x rep order
counts <- as.matrix(counts_all[, sel, drop = FALSE])
rm(counts_all)
time <- as.numeric(sub("_.*$", "", sel))
grp  <- sub("_.*$", "", sel)                     # group label = time point

## ---- gene length ----
len <- utils::read.delim(len_file, header = FALSE, stringsAsFactors = FALSE)
L <- stats::setNames(len$V2, len$V1)
counts <- counts[rownames(counts) %in% names(L), , drop = FALSE]

## ---- TPM (compute BEFORE filtering: the denominator must use all genes) ----
rpk <- counts / (L[rownames(counts)] / 1000)
tpm <- sweep(rpk, 2L, colSums(rpk), "/") * 1e6

## ---- group-level filter: any group with mean TPM > 1 ----
mean_by_grp <- sapply(split(seq_along(sel), grp), function(i)
  rowMeans(tpm[, i, drop = FALSE]))
keep <- apply(mean_by_grp, 1, max) > 1
cat("genes:", nrow(tpm), " kept:", sum(keep), "\n")
tpm <- tpm[keep, , drop = FALSE]

## ---- gene sets ----
tf_tab <- utils::read.delim(tf_file, header = FALSE, stringsAsFactors = FALSE)
tf_tab$gene <- sub("\\..*$", "", tf_tab$V1)
tf_ids <- intersect(unique(tf_tab$gene), rownames(tpm))

## DEG list is optional; when present, use DEGs union TFs as the network gene set
deg_file <- "deg_list.txt"
if (file.exists(deg_file)) {
  deg <- readLines(deg_file)
  genes <- intersect(union(deg, tf_ids), rownames(tpm))
} else {
  genes <- rownames(tpm)
}
cat("network genes:", length(genes), " TFs:", length(tf_ids), "\n")

write_mat <- function(mat, g, file) {
  out <- data.frame(id = g, mat[g, , drop = FALSE], check.names = FALSE)
  utils::write.table(out, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
write_mat(tpm, genes, file.path(outdir, "network_tpm.tsv"))
write_mat(tpm, tf_ids, file.path(outdir, "tf_tpm.tsv"))
writeLines(sel, file.path(outdir, "sample_order.txt"))
cat("DONE\n")
