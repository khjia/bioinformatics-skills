## 03_run_fast3vmrmlm.R ---- multi-locus GWAS with Fast3VmrMLM v2.0
## Env: PROJECT, TRAITS (comma-separated, or "all"), NTHREADS,
##      SVRAD (20), SVPAL (0.01), SVMLOD (3)
## Note: the pre-compiled .so may need an old TBB; if loading fails, run with
##       LD_PRELOAD=/path/to/libtbb.so.2 (see project notes).
suppressMessages(library(Fast3VmrMLM))

PROJECT  <- Sys.getenv("PROJECT"); stopifnot(nzchar(PROJECT))
INPUT    <- file.path(PROJECT, "06.Fast3VmrMLM/inputs")
OUTPUT   <- file.path(PROJECT, "06.Fast3VmrMLM/results")
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)

NTHREADS <- as.integer(Sys.getenv("NTHREADS", "40"))
SVRAD    <- as.numeric(Sys.getenv("SVRAD", "20"))
SVPAL    <- as.numeric(Sys.getenv("SVPAL", "0.01"))
SVMLOD   <- as.numeric(Sys.getenv("SVMLOD", "3"))

phe <- read.csv(file.path(INPUT, "trait.csv"), check.names = FALSE)
traits_arg <- Sys.getenv("TRAITS", "all")
idx <- if (traits_arg == "all") seq_len(ncol(phe) - 1) else {
  tr <- trimws(strsplit(traits_arg, ",")[[1]])
  m <- match(tr, names(phe)[-1])
  if (anyNA(m)) stop("unknown traits: ", paste(tr[is.na(m)], collapse = ", "))
  m
}

cat(sprintf("Fast3VmrMLM v%s | %d traits | threads %d\n",
            as.character(packageVersion("Fast3VmrMLM")), length(idx), NTHREADS))
t0 <- Sys.time()
Fast3VmrMLM(
  fileGen    = file.path(INPUT, "genotype"),
  filePhe    = file.path(INPUT, "trait.csv"),
  filePS     = file.path(INPUT, "PCA.csv"),
  PopStrType = "PC",
  fileOut    = file.path(OUTPUT, ""),
  genoType   = "SNP",
  trait      = idx,
  svrad      = SVRAD,
  svpal      = SVPAL,
  svmlod     = SVMLOD,
  nThreads   = NTHREADS,
  DrawPlot   = TRUE,
  Plotformat = "*.tiff"
)
cat(sprintf("DONE 03 in %.1f min\n", difftime(Sys.time(), t0, units = "mins")))
