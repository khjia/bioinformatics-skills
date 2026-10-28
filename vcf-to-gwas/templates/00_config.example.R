## 00_config.R ---- example configuration for the vcf-to-gwas skill (report phase)
## Copy to <analysis_dir>/scripts/00_config.R and edit the paths/parameters.

## ---- directories ----
proj     <- "/path/to/project"
wdir     <- file.path(proj, "09.postgwas2")
data_dir <- file.path(wdir, "data")     # gene_coords.rds + temporary plink files
res_dir  <- file.path(wdir, "results")  # per-trait tables
fig_dir  <- file.path(wdir, "figures")  # per-trait figures
rep_dir  <- file.path(wdir, "report")   # Rmd + rendered HTML
for (d in c(data_dir, res_dir, fig_dir, rep_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- inputs: FA / quality traits ----
trait_file <- file.path(proj, "06.Fast3VmrMLM/inputs/trait.csv")       # first col = ID, one col per trait
pca_file   <- file.path(proj, "06.Fast3VmrMLM/inputs/PCA.csv")         # ID + pc_1..pc_3
mid_dir_fa <- file.path(proj, "06.Fast3VmrMLM/results")                # <trait>_midresult.csv (MarkerID,CHR,POS,Waldst,pval)
qtl_file   <- file.path(proj, "06.Fast3VmrMLM/report/qtl_annotated.tsv") # trait,marker,chr,pos,lod,significance,gene_id,...

## ---- inputs: VE traits (optional; set report_traits_vee = character(0) to skip) ----
vee_trait_file <- file.path(proj, "08.vee/inputs/trait_vee.csv")
vee_pca_file   <- file.path(proj, "08.vee/inputs/PCA.csv")
vee_mid_dir    <- file.path(proj, "08.vee/results")
vee_qtl_file   <- file.path(proj, "08.vee/report/qtl_annotated.tsv")

## ---- annotation / genotype ----
gt_bed    <- file.path(proj, "02.ld_filter/beds/full_snps_renamed")    # PLINK1 binary
anno_file <- file.path(proj, "ref/3.annotation/function/gene_annotation_combined.xls")
## data_dir/gene_coords.rds must contain: gene_id, chr, start, end, strand
plink2    <- "/usr/local/bin/plink2"
cjk_font  <- "Noto Sans CJK SC"

## ---- trait lists & Chinese names ----
traits <- c("oil_g100g", "protein_g100g", "seed_weight_g",
            "fa_C16_0", "fa_C16_1", "fa_C18_0", "fa_C18_1", "fa_C18_2",
            "fa_C18_3n6", "fa_C18_3n3", "fa_C20_0", "fa_C20_1", "fa_C20_2",
            "fa_C22_0", "fa_C22_1", "fa_C24_0", "fa_C24_1")
trait_cn <- c(oil_g100g = "含油量", protein_g100g = "蛋白质", seed_weight_g = "千粒重",
              fa_C16_0 = "棕榈酸C16:0", fa_C16_1 = "棕榈油酸C16:1",
              fa_C18_0 = "硬脂酸C18:0", fa_C18_1 = "油酸C18:1", fa_C18_2 = "亚油酸C18:2",
              fa_C18_3n6 = "γ-亚麻酸C18:3n6", fa_C18_3n3 = "α-亚麻酸C18:3n3",
              fa_C20_0 = "花生酸C20:0", fa_C20_1 = "二十碳烯酸C20:1", fa_C20_2 = "二十碳二烯酸C20:2",
              fa_C22_0 = "山嵛酸C22:0", fa_C22_1 = "芥酸C22:1",
              fa_C24_0 = "木蜡酸C24:0", fa_C24_1 = "神经酸C24:1")

vee_cn <- c(vee_alpha = "α-生育酚", vee_beta = "β-生育酚", vee_gamma = "γ-生育酚",
            vee_delta = "δ-生育酚", vee_total = "VE总量")

## traits rendered by 22_render_all.R
report_traits_fa  <- traits
report_traits_vee <- names(vee_cn)   # use character(0) to skip VE

## ---- analysis parameters ----
MAF_MIN <- 0.05     # MAF filter when extracting regional genotypes
BONF_N  <- 2494579  # number of SNPs used for the Bonferroni threshold
