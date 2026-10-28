## 05_gene_coords.R ---- build data/gene_coords.rds from a GFF3
## Env: PROJECT, GFF, OUT
suppressMessages(library(data.table))

PROJECT <- Sys.getenv("PROJECT"); stopifnot(nzchar(PROJECT))
GFF <- Sys.getenv("GFF", file.path(PROJECT, "final.gene.mRNA.longest.gff3"))
OUT <- Sys.getenv("OUT", file.path(PROJECT, "09.postgwas2/data/gene_coords.rds"))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

g <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", shQuote(GFF)),
           header = FALSE, sep = "\t", quote = "")
setnames(g, c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attrs"))
g[, gene_id := sub(".*ID=([^;]+).*", "\\1", attrs)]
g[, gene_name := fifelse(grepl("Name=", attrs), sub(".*Name=([^;]+).*", "\\1", attrs), gene_id)]
g[, `:=`(start = as.integer(start), end = as.integer(end))]
coords <- g[, .(gene_id, gene_name, chr, start, end, strand)]
setorder(coords, chr, start)
saveRDS(coords, OUT)
cat(sprintf("gene_coords.rds: %d genes -> %s\n", nrow(coords), OUT))
