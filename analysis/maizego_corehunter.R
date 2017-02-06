if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
options(java.parameters = "-Xmx20G")
pacman::p_load("tidyverse", "data.table", "dtplyr", "viridis", "stringr",
               "corehunter")
pacman::p_load_gh("mwesthues/sspredr")

snp <- readRDS("./data/derived/maizego/numeric_snp_matrix.RDS")
snp[snp == 1] <- 2
storage.mode(snp) <- "integer"
geno <- genotypes(snp, format = "biparental")

