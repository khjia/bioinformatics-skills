## 32_variant_annotation_trait.R ---- annotate fine-mapped variants and draw gene models
## Usage: Rscript 32_variant_annotation_trait.R fa|vee <trait>
suppressMessages({library(data.table); library(Biostrings); library(ggplot2); library(patchwork)})
WDIR <- Sys.getenv("GWAS_WDIR", getwd())
source(file.path(WDIR, "scripts/00_config.R"))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; TR <- args[2]
cat(sprintf("== variant annotation %s %s\n", prefix, TR))

OUT <- file.path(res_dir, "variant_annotation"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FIGD <- file.path(fig_dir, paste0(prefix, "_", TR, "_report")); dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

sf <- file.path(res_dir, paste0(prefix, "_", TR, "_report"), "tab_summary.tsv")
if (!file.exists(sf)) { cat("no tab_summary\n"); quit(status = 0) }
st <- fread(sf); CHR <- as.character(st$anchor_chr[1]); APOS <- as.numeric(st$anchor_pos[1])
if (is.na(CHR) || !nzchar(CHR)) { cat("no anchor\n"); quit(status = 0) }
FROM <- max(1, APOS - 1e6); TO <- APOS + 1e6

sdir <- file.path(res_dir, paste0("susie_", TR))
if (!file.exists(file.path(sdir, "pip.tsv"))) { cat("no susie output\n"); quit(status = 0) }
pip  <- fread(file.path(sdir, "pip.tsv"))
cset <- fread(file.path(sdir, "credible_sets.tsv"))

vsnp <- unique(c(cset$top_snp, pip[order(-pip)][1:10, snp]))
v <- pip[snp %in% vsnp]
v[, cs_id := fifelse(snp %in% cset$top_snp, paste0("CS", cset$CS[match(snp, cset$top_snp)]), "")]

FA  <- file.path(proj, "genome.a.fasta")
GFF <- file.path(proj, "final.gene.mRNA.longest.gff3")
anno <- fread(anno_file, select = c("gene", "Preferred_name", "Pfam", "A.thaliana",
                                    "eggNOG free text description"))
fa_chr <- file.path(data_dir, paste0(CHR, ".fa"))
if (!file.exists(fa_chr)) system(sprintf("samtools faidx %s %s > %s", shQuote(FA), CHR, shQuote(fa_chr)))
dna <- readDNAStringSet(fa_chr); names(dna) <- CHR

cmd <- sprintf("awk -F'\\t' '$1==\"%s\" && $4>=%d && $5<=%d' %s", CHR, FROM, TO, shQuote(GFF))
g <- fread(cmd = cmd, header = FALSE, sep = "\t", quote = "")
setnames(g, c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attrs"))
g[, id := sub(".*ID=([^;]+).*", "\\1", attrs)]
g[, parent := sub(".*Parent=([^;]+).*", "\\1", attrs)]
genes <- g[type == "gene", .(gene_id = sub("^gene:", "", id), gstart = start, gend = end, strand)]
mrna  <- g[type == "mRNA", .(tx = sub("^transcript:", "", id), gene_id = sub("^gene:", "", parent))]
feat  <- g[type %in% c("exon", "CDS", "five_prime_UTR", "three_prime_UTR"),
           .(tx = sub("^transcript:", "", parent), type, start, end)]
feat  <- merge(feat, mrna, by = "tx", all.x = TRUE)[!is.na(gene_id)]

tx_cds <- function(gid) {
  txid <- mrna[gene_id == gid]$tx[1]
  cd <- feat[gene_id == gid & tx == txid & type == "CDS"][order(start)]
  st <- genes[gene_id == gid]$strand[1]
  if (!nrow(cd)) return(NULL)
  if (st == "-") cd <- cd[order(-start)]
  segs <- lapply(seq_len(nrow(cd)), function(i) {
    s <- subseq(dna[[CHR]], cd$start[i], cd$end[i])
    if (st == "-") reverseComplement(s) else s
  })
  list(cd = cd, seq = DNAString(paste(vapply(segs, as.character, character(1)), collapse = "")),
       lens = cumsum(width(DNAStringSet(segs))), strand = st)
}
codon_change <- function(gid, pos, ref, alt) {
  x <- tx_cds(gid); if (is.null(x)) return(list(codon = NA, aa = NA))
  cd <- x$cd; sp <- x$seq; st <- x$strand
  lens <- x$lens; idx <- NA_integer_
  for (i in seq_len(nrow(cd))) {
    off0 <- if (i == 1) 0 else lens[i - 1]
    if (pos >= cd$start[i] && pos <= cd$end[i]) {
      idx <- off0 + if (st == "+") pos - cd$start[i] else cd$end[i] - pos; break
    }
  }
  if (is.na(idx)) return(list(codon = NA, aa = NA))
  ci <- floor(idx / 3) * 3
  cod <- as.character(subseq(sp, ci + 1, ci + 3))
  if (nchar(cod) < 3) return(list(codon = cod, aa = NA))
  altb <- if (st == "-") as.character(complement(DNAString(alt))) else alt
  codv <- strsplit(cod, "")[[1]]; codv[idx %% 3 + 1] <- altb
  aa1 <- as.character(translate(DNAString(cod), genetic.code = GENETIC_CODE))
  aa2 <- as.character(translate(DNAString(paste(codv, collapse = "")), genetic.code = GENETIC_CODE))
  list(codon = sprintf("%s>%s", cod, paste(codv, collapse = "")),
       aa = sprintf("p.%s%d%s", aa1, ci/3 + 1, aa2))
}
classify <- function(pos, ref, alt) {
  ins <- genes[gstart <= pos & gend >= pos]
  if (nrow(ins)) gid <- ins$gene_id[1] else {
    d <- pmin(abs(genes$gstart - pos), abs(genes$gend - pos)); j <- which.min(d)
    if (d[j] > 50000) return(list(gene = NA, cls = "intergenic", dist = d[j], tss = NA, codon = NA, aa = NA))
    gid <- genes$gene_id[j]
  }
  gg <- genes[gene_id == gid][1]; st <- gg$strand; ft <- feat[gene_id == gid]
  cds <- ft[type == "CDS"]; ex <- ft[type == "exon"]
  u5 <- ft[type == "five_prime_UTR"]; u3 <- ft[type == "three_prime_UTR"]
  tss <- if (st == "+") gg$gstart else gg$gend; dtss <- abs(pos - tss)
  up <- if (st == "+") pos < gg$gstart else pos > gg$gend
  down <- if (st == "+") pos > gg$gend else pos < gg$gstart
  codon <- NA; aa <- NA
  if (nrow(cds[start <= pos & end >= pos])) {
    cls <- "CDS"; cc <- codon_change(gid, pos, ref, alt); codon <- cc$codon; aa <- cc$aa
  } else if (nrow(u5[start <= pos & end >= pos])) cls <- "5'UTR"
  else if (nrow(u3[start <= pos & end >= pos])) cls <- "3'UTR"
  else if (nrow(ex[start <= pos & end >= pos])) cls <- "exon_noncoding"
  else if (pos >= gg$gstart && pos <= gg$gend) cls <- "intron"
  else if (up && dtss <= 2000) cls <- "promoter"
  else if (down && dtss <= 2000) cls <- "downstream"
  else cls <- "intergenic"
  list(gene = gid, cls = cls, dist = if (nrow(ins)) 0L else min(abs(pos - gg$gstart), abs(pos - gg$gend)),
       tss = dtss, codon = codon, aa = aa)
}

parts <- tstrsplit(v$snp, ":", fixed = TRUE)
v[, `:=`(ref = parts[[3]], alt = parts[[4]])]
res <- rbindlist(lapply(seq_len(nrow(v)), function(i) {
  z <- classify(v$pos[i], v$ref[i], v$alt[i])
  data.table(snp = v$snp[i], pos = v$pos[i], cs = v$cs_id[i], pip = round(v$pip[i], 3),
             ref = v$ref[i], alt = v$alt[i], gene = z$gene, region = z$cls,
             dist_gene_bp = z$dist, dist_tss_bp = z$tss, codon = z$codon, aa_change = z$aa)
}))
res <- merge(res, genes[, .(gene = gene_id, strand)], by = "gene", all.x = TRUE)
res <- merge(res, anno, by.x = "gene", by.y = "gene", all.x = TRUE)
setorder(res, -pip)
fwrite(res, file.path(OUT, paste0(TR, "_variant_annotation.tsv")), sep = "\t")
print(res[, .(snp, cs, pip, gene, region, dist_tss_bp, codon, aa_change)])

## ---- gene-model figure around the top 3 variants ----
top3 <- head(res[!is.na(gene)], 3)
if (nrow(top3)) {
  mx <- function(pos, from, to) (pos - from) / (to - from)
  plist <- lapply(seq_len(nrow(top3)), function(i) {
    p <- top3$pos[i]; from <- p - 20000; to <- p + 20000
    gd <- genes[gend >= from & gstart <= to][order(gstart)]
    gd[, y := seq_len(.N)]
    ex <- merge(feat[type == "exon"], gd[, .(gene_id, y)], by = "gene_id")
    ggplot() +
      geom_segment(data = gd, aes(mx(gstart, from, to), y, xend = mx(gend, from, to), yend = y),
                   colour = "#E69F00", linewidth = 2.2) +
      geom_rect(data = ex[!is.na(y)], aes(xmin = mx(start, from, to), xmax = mx(end, from, to),
                                          ymin = y - 0.24, ymax = y + 0.24), fill = "#7B2D8E", colour = NA) +
      geom_vline(xintercept = mx(p, from, to), colour = "#C00000", linetype = 2, linewidth = 0.6) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
                         labels = sprintf("%.2f", c(from, p, to) / 1e6)) +
      scale_y_continuous(expand = expansion(add = c(0.6, 0.8))) +
      labs(x = sprintf("%s 位置 (Mb)", CHR), y = NULL,
           subtitle = sprintf("%s（%s，%s）", sub("^.*:", "", sub(":.*$", "", top3$snp[i])), top3$gene[i], top3$region[i])) +
      theme_bw(base_size = 9) +
      theme(text = element_text(family = cjk_font), panel.grid.minor = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            plot.subtitle = element_text(size = 8.5))
  })
  pv <- wrap_plots(plist, ncol = 1)
  ggsave(file.path(FIGD, "variant_focus.png"), pv, width = 8.2, height = 2.6 * nrow(top3), dpi = 220,
         device = grDevices::png, type = "cairo")
}
cat("DONE 32\n")
