# Goal: Sample core sets of genotypes for the Yan-lab data and scenario two of
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
options(java.parameters = "-Xmx25G")
pacman::p_load("tidyverse", "data.table", "dtplyr", "corehunter")
pacman::p_load_gh("mwesthues/sspredr")

# For the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- common_genotypes %>%
  reduce(intersect)
snp <- readRDS("./data/processed/maizego/imputed_snp_mat.RDS")
snp <- snp[rownames(snp) %in% common_genotypes, ]

## --- Scenario 2 SNP preparation.
# Explore the kinship matrix.
geno <- snp %>%
  genotypes(format = "biparental") %>%
  coreHunterData(genotypes = .)
saveRDS(geno, "./data/derived/maizego/scenario2_biparental_genotypes.RDS")
geno <- readRDS("./data/derived/maizego/scenario2_biparental_genotypes.RDS")

core_seq <- seq(from = 0.1, to = 0.9, by = 0.1)
core_lst <- core_seq %>%
  map(., ~sampleCore(data = geno,
           obj = objective(type = "EN", measure = "MR", weight = 1),
           size = .)
  )
names(core_lst) <- paste0("selected_fraction_", core_seq)

# Add the full set of genotypes that are covered by mRNA data to the list for 
# the equivalent of a "1.0" scenario to exist.
# This is of interest because it simplifies the prediction script.
core_lst[["selected_fraction_1.0"]][["sel"]] <- common_genotypes
core_lst[["selected_fraction_1.0"]][["EN"]][["MR"]] <- NA_real_
saveRDS(core_lst, "./data/derived/maizego/scenario2_snp_core_list.RDS")
