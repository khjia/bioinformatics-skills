## 33_known_loci.R ---- Europe PMC literature check for the leading candidate gene
## Usage: Rscript 33_known_loci.R fa|vee <trait>
suppressMessages({library(data.table); library(jsonlite)})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
WDIR <- Sys.getenv("GWAS_WDIR", getwd())
source(file.path(WDIR, "scripts/00_config.R"))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; TR <- args[2]
cat(sprintf("== known loci %s %s\n", prefix, TR))

OUT <- file.path(res_dir, "known_loci"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cf <- file.path(res_dir, paste0(prefix, "_", TR, "_report"), "tab_candidate_genes.tsv")
if (!file.exists(cf)) { cat("no candidate table\n"); quit(status = 0) }
cand <- fread(cf)
if (!nrow(cand)) { cat("empty candidate table\n"); quit(status = 0) }

KEY <- c(oil_g100g = "seed oil content", protein_g100g = "seed protein content",
         seed_weight_g = "seed weight", fa_C16_0 = "palmitic acid", fa_C16_1 = "palmitoleic acid",
         fa_C18_0 = "stearic acid", fa_C18_1 = "oleic acid", fa_C18_2 = "linoleic acid",
         fa_C18_3n3 = "alpha-linolenic acid", fa_C18_3n6 = "gamma-linolenic acid",
         fa_C20_0 = "arachidic acid", fa_C20_1 = "eicosenoic acid", fa_C20_2 = "eicosadienoic acid",
         fa_C22_0 = "behenic acid", fa_C22_1 = "erucic acid", fa_C24_0 = "lignoceric acid",
         fa_C24_1 = "nervonic acid", vee_alpha = "alpha-tocopherol", vee_beta = "beta-tocopherol",
         vee_gamma = "gamma-tocopherol", vee_delta = "delta-tocopherol",
         vee_total = "vitamin E tocopherol")
kw <- if (TR %in% names(KEY)) KEY[[TR]] else "seed oil fatty acid"
gene <- cand$gene_id[1]

at <- NA_character_
if ("at_id" %in% names(cand)) at <- cand$at_id[1]
if (is.na(at) || at == "-" || at == "") {
  am <- fread(anno_file, select = c("gene", "A.thaliana"))
  hit <- am[gene == cand$gene_id[1], A.thaliana]
  if (length(hit) && nzchar(hit[1]) && hit[1] != "-")
    at <- sub("^([A-Za-z0-9]+)\\..*$", "\\1", hit[1])
}
if (is.na(at) || at == "-" || at == "") at <- gene

run_query <- function(q) {
  url <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%s&format=json&pageSize=20&resultType=core&sort=CITED%%20desc",
                 utils::URLencode(q, reserved = TRUE))
  tmp <- tempfile(fileext = ".json")
  cmd <- sprintf("curl -sS --max-time 90 --retry 3 --retry-delay 3 --retry-all-errors '%s' -o '%s'", url, tmp)
  rc <- system(cmd)
  if (rc != 0 || !file.exists(tmp) || file.info(tmp)$size < 50) return(NULL)
  js <- tryCatch(fromJSON(tmp, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(js)) return(NULL)
  hits <- js$resultList$result
  rows <- rbindlist(lapply(hits, function(h) data.table(
    pmid = h$pmid %||% "", doi = h$doi %||% "", year = h$pubYear %||% "",
    title = gsub("\t", " ", h$title %||% ""),
    journal = (h$journalInfo$journal$title %||% ""),
    first_author = strsplit(h$authorString %||% "", ",")[[1]][1],
    cited_by = h$citedByCount %||% 0)))
  if (!nrow(rows)) rows <- data.table(pmid = character(), doi = character(), year = character(),
                                      title = character(), journal = character(),
                                      first_author = character(), cited_by = numeric())
  list(rows = rows, n = js$hitCount %||% nrow(rows))
}

q1 <- sprintf("%s AND (%s)", at, kw)
cat(sprintf("query 1: %s\n", q1))
r <- run_query(q1)
used <- q1
if (is.null(r) || !nrow(r$rows) || r$n == 0) {
  q2 <- at
  cat(sprintf("query 2 (gene only): %s\n", q2))
  r2 <- run_query(q2)
  if (!is.null(r2)) { r <- r2; used <- q2 }
}
if (is.null(r)) { cat("query failed\n"); quit(status = 0) }

fwrite(r$rows, file.path(OUT, paste0(TR, "_known.tsv")), sep = "\t")
meta <- data.table(gene = at, planted_gene = gene, trait = TR, keywords = kw,
                   query_used = used, n_hits = r$n, returned = nrow(r$rows), source = "EuropePMC")
fwrite(meta, file.path(OUT, paste0(TR, "_meta.tsv")), sep = "\t")
cat(sprintf("%d papers (hitCount %s)\n", nrow(r$rows), r$n))
cat("DONE 33\n")
