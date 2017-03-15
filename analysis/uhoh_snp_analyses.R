# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis")


snp <- readRDS("./data/processed/maizego/imputed_snp_mat.RDS")
unique_genotypes <- "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS" %>%
  readRDS() %>%
  keep(names(.) == "snp") %>%
  flatten_chr()


# PCA ---------------------------------------------------------------------
write.lfmm(snp, "./data/derived/maizego/snp.lfmm")

# Determine the structure of the data using genotypic data.
snp_pc <- LEA::pca(input.file = "./data/derived/maizego/snp.lfmm",
                   scale = TRUE)
# Perform Tracy-Widom tests on all eigenvalues to determine the optimal number
# of principal components.
snp_tw <- tracy.widom(snp_pc)
tw_K <- snp_tw %>%
  .["pvalues"] %>%
  flatten_dbl() %>%
  keep(~ .x <= 0.0001) %>%
  length()
tw_K <- min(c(tw_K, 5))

snp_tw %>%
  .["percentage"] %>%
  flatten_dbl() %>%
  plot()

# Qualitatively assign the genotypes to the previously determined 
# subpopulations.
lfmm2geno("./data/derived/maizego/snp.lfmm", 
          output.file = "./data/derived/maizego/snp.geno")
project <- snmf(input.file = "./data/derived/maizego/snp.geno",
                K = seq_len(tw_K),
                project = "new",
                repetitions = 25,
                CPU = 3,
                entropy = TRUE)

# Based on the results of a PCA, decide whether it is even necessary to treat
# subpopulations as such or whether all inbred lines can be treated as a single
# group.
project <- load.snmfProject("./data/derived/maizego/snp.snmfProject")


# Get labels for the different subpopulations of the genotypes.
pop_struc <- fread("./data/input/maizego/genotype_annotation.txt")

pc_mat <- snp_pc$projections
rownames(pc_mat) <- rownames(snp)
pc_with_assignment <- pc_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "Lines") %>%
  gather(key = PC, value = Score, -Lines) %>%
  left_join(y = pop_struc, by = "Lines") %>%
  mutate(PC = str_replace_all(PC, pattern = "V", replacement = "PC"))

pc_with_assignment %>%
  filter(PC %in% c("PC1", "PC2")) %>%
  spread(key = PC, value = Score) %>%
  mutate(Lines = as.factor(Lines)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = Origin)) +
  geom_point() 

pc_with_assignment %>%
  filter(PC %in% c("PC1", "PC3")) %>%
  spread(key = PC, value = Score) %>%
  mutate(Lines = as.factor(Lines)) %>%
  ggplot(., aes(x = PC1, y = PC3, color = Origin)) +
  geom_point() 



# ADMIXTURE ---------------------------------------------------------------
K_lst <- seq(from = 2, to = tw_K, by = 1) %>%
  as.list()

# Function to compute the cross entropy.
ce_fun <- function(x, k) {
  x %>%
    cross.entropy(K = k) %>%
    which.min()
}

# For each K, get the run with the smallest cross entropy criterion.
ce_lst <- K_lst %>%
  map(., .f = ~ce_fun(x = project, k = .))

# Return the admixture coefficients for each chosen run with K ancestral 
# populations.
structure_df <- K_lst %>%
  map2(.y = ce_lst, .f = ~Q(object = project, K = .x, run = .y)) %>%
  map(as_data_frame) %>%
  map(function(x) {
    x$G <- rownames(pc_mat)
    x
  }) %>%
  map(., ~gather(., key = Component, value = Value, -G)) %>%
  bind_rows(.id = "K")

# For each genotype and for each run with K ancestral populations, plot the 
# admixture coefficients by genotype.
structure_df %>%
  mutate(
    K = as.numeric(K),
    K = K + 1,
    K = paste0("K=", as.character(K))
  ) %>%
  mutate_at(vars(K, G, Component), as.factor) %>%
  ggplot(aes(x = G, y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~K, ncol = 1, strip.position = "right") +
  scale_fill_viridis(discrete = TRUE) +
  ggthemes::theme_pander() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")
