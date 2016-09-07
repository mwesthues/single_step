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



# CV SIMILARITY -----------------------------------------------------------
cv_scheme <- readRDS(
  "./data/processed/cv1000_ps8081_trn=500_min_size=50_m=114_f=83.RDS"
  )
sets <- c("T0", "T1", "T2", "TRN")
cv_mat <- cv_scheme %>%
  as_data_frame() %>%
  spread(key = Run, value = Set) %>%
  remove_rownames() %>%
  as.data.frame() %>%
  column_to_rownames(var = "Sample_ID") %>%
  as.matrix()
# Function for the computation of the absolute correlation between each 
# cross-validation run, separately for each set (T0, T1, T2, TRN).
get_set_sim <- function(x, set) {
  x[x == set] <- "1"
  x[x != "1"] <- "0"
  storage.mode(x) <- "numeric"
  x %>%
    cor() %>%
    abs() %>%
    .[upper.tri(.)] %>%
    c()
}
set_sim_lst <- lapply(sets, FUN = function(set) {
  get_set_sim(x = cv_mat, set = set)
})
names(set_sim_lst) <- sets
sim_df <- do.call(cbind, set_sim_lst)
# Histograms of the similarity between CV-runs, separately for each set.
sim_df %>%
  as_data_frame() %>%
  gather(key = Set, value = Similarity) %>%
  ggplot(aes(x = Similarity)) +
  geom_histogram() +
  facet_wrap(~ Set)





















