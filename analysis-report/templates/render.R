#!/usr/bin/env Rscript
## Render an analysis report to a self-contained HTML.
## Usage: Rscript render.R report.Rmd "项目_分析内容_2026-11-05.html"
##
## Pitfall: rmarkdown needs pandoc on PATH. If pandoc lives in a conda env,
## prepend that env's bin BEFORE calling Rscript, e.g.
##   export PATH=/path/to/envs/xxx/bin:$PATH

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript render.R <input.Rmd> [output.html]")

input <- args[1]
output <- if (length(args) >= 2) args[2] else sub("\\.Rmd$", ".html", basename(input))

if (!nzchar(Sys.which("pandoc"))) {
  cand <- c("/usr/local/bin", "/usr/bin",
            file.path(Sys.getenv("CONDA_PREFIX"), "bin"))
  add <- cand[file.exists(file.path(cand, "pandoc"))]
  if (length(add)) Sys.setenv(PATH = paste(c(add, Sys.getenv("PATH")), collapse = .Platform$path.sep))
}

if (!rmarkdown::pandoc_available()) {
  stop("pandoc not found. Activate the conda env that ships pandoc before running this script.")
}

rmarkdown::render(
  input,
  output_file = output,
  output_format = rmarkdown::html_document(
    toc = TRUE, toc_float = list(collapsed = FALSE), toc_depth = 3,
    theme = "cosmo", highlight = "tango", self_contained = TRUE
  ),
  encoding = "UTF-8"
)

cat("\nOutput:", output, "\n")
