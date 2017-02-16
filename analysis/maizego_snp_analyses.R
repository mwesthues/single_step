# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis")
snp <- readRDS("./data/processed/maizego/imputed_snp_mat.RDS")
unique_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
  ) %>%
  keep(names(.) == "snp") %>%
  flatten_chr()


# PCA ---------------------------------------------------------------------
write.lfmm(snp, "./data/derived/maizego/snp.lfmm")

# Determine the structure of the data using genotypic and gene expression data,
# respectively.
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
tw_K <- min(c(tw_K, 10))

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
                repetitions = 10,
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


# Evaluate the population structure when considering the first K principal 
# components.
K <- 10
ce <- cross.entropy(project, K = K)
best <- which.min(ce)
barplot(t(Q(project, K = K, run = best)), col = viridis(n = K))


