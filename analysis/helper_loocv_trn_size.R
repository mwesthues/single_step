if (!require("pacman")) install.packages("pacman")
pacman::p_load("Matrix", "data.table", "tidyr", "matrixStats", "dtplyr", 
               "dplyr", "tibble", "ggplot2")

pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
# Common genotypes
geno <- readRDS("./data/processed/common_genotypes.RDS")
pheno <- pheno %>%
  filter(G %in% geno$Hybrid) %>%
  select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  select(-ADL) %>%
  as.data.frame %>%
  column_to_rownames(var = "G") %>%
  as.matrix

geno <- rownames(pheno)
# Position of Dent in rownames(pheno): 2
# Position of Flint in rownames(pheno): 3
mother_idx <- 2
father_idx <- 3

runs <- seq_len(nrow(pheno))
# Paternal IDs
father <- vapply(strsplit(geno, split = "_"), FUN = "[[",
                 father_idx,
                 FUN.VALUE = character(1))
# Maternal IDs
mother <- vapply(strsplit(geno, split = "_"), FUN = "[[",
                 mother_idx,
                 FUN.VALUE = character(1))

trn_fraction <- vapply(seq_len(nrow(pheno)), FUN = function(run) {
  # Use only T0 hybrids for the training set.
  tst <- intersect(grep(mother[run], x = geno, invert = TRUE),
                   grep(father[run], x = geno, invert = TRUE))
  length(tst) / nrow(pheno)
}, FUN.VALUE = double(1))
