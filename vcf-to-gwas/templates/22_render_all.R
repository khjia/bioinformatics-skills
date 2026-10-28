## 22_render_all.R ---- generate figures + render one report per trait
## Env vars: GWAS_WDIR (analysis dir), GWAS_CONFIG (path to 00_config.R)
## Config must define: report_traits_fa, report_traits_vee (optional), vee_cn (optional)
Sys.setenv(RSTUDIO_PANDOC = Sys.getenv("RSTUDIO_PANDOC", "/usr/bin"))
suppressMessages(library(rmarkdown))
W <- Sys.getenv("GWAS_WDIR", getwd())
CFG <- Sys.getenv("GWAS_CONFIG", file.path(W, "scripts/00_config.R"))
source(CFG)

traits_fa  <- if (exists("report_traits_fa"))  report_traits_fa  else traits
traits_vee <- if (exists("report_traits_vee")) report_traits_vee else character(0)
cn_vee     <- if (exists("vee_cn")) vee_cn else setNames(traits_vee, traits_vee)
safe <- function(x) gsub("[ /\\\\*?\"<>|]", "", gsub(":", "-", x))

todo <- rbind(
  data.table::data.table(prefix = "fa",  trait = traits_fa,  cn = unname(trait_cn[traits_fa])),
  data.table::data.table(prefix = "vee", trait = traits_vee, cn = unname(cn_vee[traits_vee])))

for (i in seq_len(nrow(todo))) {
  p <- todo$prefix[i]; tr <- todo$trait[i]; cn <- todo$cn[i]
  cat(sprintf("[%d/%d] %s %s (%s)\n", i, nrow(todo), p, tr, cn)); flush.console()
  RSCRIPT <- Sys.getenv("RSCRIPT", "Rscript")
  r <- system2(RSCRIPT, c(file.path(W, "scripts/20_trait_report_figs.R"), p, tr),
               stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(r, "status")) && attr(r, "status") != 0) {
    cat("FIG FAIL:", tr, "\n", tail(r, 5), "\n"); next
  }
  out <- sprintf("元宝枫_%s（%s）_GWAS单性状报告_%s.html", safe(cn), tr, Sys.Date())
  ok <- tryCatch({
    render(file.path(W, "report/21_trait_report.Rmd"), output_file = out,
           output_dir = file.path(W, "report"),
           params = list(prefix = p, trait = tr, cn = cn),
           envir = new.env(parent = globalenv()), quiet = TRUE)
    TRUE
  }, error = function(e) { cat("RENDER FAIL:", tr, conditionMessage(e), "\n"); FALSE })
  cat(if (ok) "OK\n" else "FAIL\n"); flush.console()
}
cat("ALL DONE\n")
