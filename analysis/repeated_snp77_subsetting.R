# Goal: Generate 20 subsets of genomic data to average out the influence of 
# sampling on all prediction scenarios including 'snp77'. 

if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("BGLR","data.table", "parallel", "magrittr", "dplyr", "tidyr",
               "purrr", "testthat", "methods", "tibble", "ggplot2", "stringr",
               "stringi", "lubridate", "readr")
devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")



# -- GENOTYPE INFORMATION -----------------------------------------------
# Information on the set of genotypes.
genos <- readRDS("./data/processed/common_genotypes.RDS")

# -- LOAD PREDICTOR DATA ------------------------------------------------
# pedigree
cat_lst <- readRDS("./data/derived/pred_sub_list.RDS")
hybrids <- cat_lst %>%
  purrr::transpose() %>%
  .[names(.) == "ped100_snp77"] %>%
  at_depth(.depth = 2, .f = "geno") %>%
  .[[1]]

ped100 <- cat_lst %>% 
  purrr::transpose() %>%
  .[names(.) == "ped100_snp77"] %>%
  at_depth(.depth = 2, .f = 1) %>%
  .[[1]]
names(ped100) <- c("Dent", "Flint")

# snp
snp100 <- readRDS("./data/processed/snp_mat.RDS")
snp_nms <- genos %>%
  filter(Data_Type == "snp") %>%
  .$G
snp100 <- snp100 %>%
  .[rownames(.) %in% snp_nms, ]

# mrna
mrna42 <- readRDS("./data/derived/predictor_subsets/transformed-mrna42.RDS")



# -- RESAMPLE PREDICTOR DATA -------------------------------------------- 
# Remove some genomic records so that we can impute them via pedigree
# information.
mrna_genos <- genos %>%
  filter(Data_Type == "mrna") %>%
  split(.$Pool) %>%
  map("G")

set.seed(314)
genos %>%
  filter(Pool %in% c("Dent", "Flint"),
         Data_Type == "snp",
         !G %in% flatten_chr(mrna_genos)) %>%
  split(.$Pool) %>%
  at_depth(1, "G") %>%
  map(~rerun(20, sample(., size = length(.) * 0.6))) %>%
  map(~do.call(cbind, .)) %>%
  map(as_data_frame) %>%
  bind_rows(.id = "Group") %>%
  gather(key = Run, value = Genotype, -Group) -> sample_df

sample_df %>%
  split(.$Run) %>%
  map(~split(., .$Group)) %>%
  at_depth(.depth = 2, "Genotype") -> geno_lst

geno_lst %>%
  transpose() %>%
  names() -> grp_nms

storage_dir <- "./data/derived/predictor_subsets/snp77_repetition_"
lapply(seq_along(geno_lst), FUN = function(i) {
  nm_lst <- geno_lst[[i]]
  snp77 <- lapply(grp_nms, FUN = function(j) {
    snp_geno <- nm_lst[[j]]
    mrna_geno <- mrna_genos[[j]]
    snp35 <- snp100[snp_geno, ]
    snp_mat <- snp100 %>%
      .[rownames(.) %in% mrna_geno, ] %>%
      rbind(., snp35) %>%
      unique(., MARGIN = 2)
    snp_mat
  })
  names(snp77) <- grp_nms
  saveRDS(snp77, file = paste0(storage_dir, i, ".RDS"))
})


