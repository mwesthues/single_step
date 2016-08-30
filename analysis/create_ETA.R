if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("data.table", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
Pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
Pheno <- Pheno[G %in% genos$Hybrid, , ]

# Endophenotypes
snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna <- mrna[rownames(mrna) %in% rownames(snp), ]


imputation <- TRUE

# Build a kernel from the genomic data
hetgrps <- c("Dent", "Flint")
if (isTRUE(imputation)) {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    grp_geno <- genos[[hetgrp]][["snp"]]
    x <- snp[rownames(snp) %in% grp_geno, ]
    y <- mrna[rownames(mrna) %in% grp_geno, ]
    hyb_id <- ifelse(hetgrp == "Dent", yes = 2, no = 3)
    geno <- vapply(strsplit(genos$Hybrid, split = "_"), FUN = "[[", hyb_id, 
                   FUN.VALUE = character(1))
    impute_eta(x = x,
               y = y,
               geno = geno,
               bglr_model = "BRR")
  })
} else {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    grp_geno <- genos[[hetgrp]][["snp"]]
    x <- snp[rownames(snp) %in% grp_geno, ]
    hyb_id <- ifelse(hetgrp == "Dent", yes = 2, no = 3)
    geno <- vapply(strsplit(genos$Hybrid, split = "_"), FUN = "[[", hyb_id, 
                   FUN.VALUE = character(1))
    complete_eta(x = x, 
                 geno = geno, 
                 bglr_model = "BRR")
  })
}
ETA <- unlist(ETA, recursive = FALSE)



