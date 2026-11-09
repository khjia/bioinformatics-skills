## 41_render_all.R ---- render one HTML report per trait from the step-based template
Sys.setenv(RSTUDIO_PANDOC = "/media/nfs1/hermes/miniforge3/bin")
suppressMessages(library(rmarkdown))
W <- "/media/nfs1/hermes/project/yuanbaofeng/09.postgwas2"
Sys.setenv(GWAS_WDIR = W)
source(file.path(W, "scripts/00_config.R"))

traits_fa <- c("oil_g100g", "protein_g100g", "seed_weight_g",
               "fa_C16_0", "fa_C16_1", "fa_C18_0", "fa_C18_1", "fa_C18_2",
               "fa_C18_3n6", "fa_C18_3n3", "fa_C20_0", "fa_C20_1", "fa_C20_2",
               "fa_C22_0", "fa_C22_1", "fa_C24_0", "fa_C24_1")
traits_vee <- c("vee_alpha", "vee_beta", "vee_gamma", "vee_delta", "vee_total")
cn_vee <- c(vee_alpha = "α-生育酚", vee_beta = "β-生育酚", vee_gamma = "γ-生育酚",
            vee_delta = "δ-生育酚", vee_total = "VE 总量")
safe <- function(x) gsub("[ /\\\\*?\"<>|]", "", gsub(":", "-", x))

todo <- rbind(
  data.table::data.table(prefix = rep("fa", length(traits_fa)), trait = traits_fa,
                         cn = unname(trait_cn[traits_fa])),
  data.table::data.table(prefix = rep("vee", length(traits_vee)), trait = traits_vee,
                         cn = unname(cn_vee[traits_vee])))

for (i in seq_len(nrow(todo))) {
  p <- todo$prefix[i]; tr <- todo$trait[i]; cn <- todo$cn[i]
  cat(sprintf("[%d/%d] %s %s (%s)\n", i, nrow(todo), p, tr, cn)); flush.console()
  out <- sprintf("元宝枫_%s_单性状GWAS分析报告_%s.html", safe(cn), format(Sys.Date(), "%Y-%m-%d"))
  ok <- tryCatch({
    render(file.path(W, "report/21_trait_report_v2.Rmd"), output_file = out,
           output_dir = file.path(W, "report"),
           params = list(prefix = p, trait = tr, cn = cn),
           envir = new.env(parent = globalenv()), quiet = TRUE)
    TRUE
  }, error = function(e) { cat("RENDER FAIL:", tr, conditionMessage(e), "\n"); FALSE })
  cat(if (ok) "OK\n" else "FAIL\n"); flush.console()
}
cat("ALL RENDER DONE\n")
