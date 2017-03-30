# Goal: Run a pricipal component analysis on the SNP marker data of the 
# parental inbred lines of the UHOH maize hybrids.
# Use the plot of the first two principal components to assess how well the 
# total genetic space of the material is covered by genotypes that have SNP and
# mRNA data.

# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis", "devtools", "ggthemes")
devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")


# ************************ PARENTAL INBRED LINES ***************************
# Group-wise SNP quality checks.
# Load SNP marker data, which were imputed using Beagle.
snp <- "./data/processed/uhoh/imputed_snp_mat.RDS" %>%
  readRDS()
rownames(snp) <- gsub("X", replacement = "", x = rownames(snp))

# Data frame with information group membership (Dent or Flint) of genotypes and
# with information on predictor coverage of genotypes.
common_genotypes <- "./data/processed/common_genotypes.RDS" %>%
  readRDS()
  
# For each heterotic group, apply quality checks to the SNP data and keep only
# markers that pass those criteria.
quality_snps <- common_genotypes %>%
  filter(Data_Type == "snp") %>%
  split(.$Pool) %>%
  map("G") %>%
  map(function(x) {
    snp[rownames(snp) %in% x, ]
  }) %>%
  map(., ~ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  ))


# Save the genotypic data to files so that the PCA can be run with the pca 
# function from the LEA package.
lfmm_path <- "./data/derived/uhoh/"
quality_snps %>%
  seq_along() %>%
  map(function(i) {
    x <- quality_snps[[i]]
    x_nm <- names(quality_snps)[i]
    write.lfmm(x, paste0(lfmm_path, "snp_", x_nm, ".lfmm"))
  })

# Determine the structure of the data using genotypic data.
pca_lst <- lfmm_path %>%
  list.files(pattern = "lfmm") %>%
  discard(. == "snp_Hybrid.lfmm") %>%
  map(., ~paste0(lfmm_path, .)) %>%
  map(., ~LEA::pca(input.file = ., scale = TRUE)) %>%
  map(function(x) x$projections) %>%
  map2(.y = quality_snps, .f = function(x, y) {
    rownames(x) <- rownames(y)
    x
  })

# Get the names of genotypes that are covered only by mRNA data and the ones
# that are covered by both, SNP and mRNA data.
# These information will be required for assigning genotypes to the complete 
# (i.e., SNP) and the incomplete (i.e., mRNA) predictor in the PCA.
predictor_coverage <- common_genotypes %>%
  filter(Data_Type %in% c("snp", "mrna")) %>%
  split(.$Data_Type) %>%
  map("G")

# Convert the PCA list for Dent and Flint material to a data frame and add 
# information on predictor coverage of the genotypes.
hybrids <- pca_lst %>%
  map(as.data.frame) %>%
  map(., ~rownames_to_column(., var = "G")) %>%
  map(., ~gather(., key = PC, value = Score, -G)) %>%
  set_names(names(quality_snps)) %>%
  bind_rows(.id = "Group") %>%
  as_data_frame() %>%
  mutate(PC = str_replace_all(PC, pattern = "V", replacement = "PC")) %>%
  filter(PC %in% paste0("PC", seq_len(2))) %>%
  mutate(
    Cluster = if_else(
      G %in% predictor_coverage$mrna,
      true = "Core",
      false = "Full"
    )
  )


# PCA plot
g1 <- pca_df %>%
  spread(key = PC, value = Score) %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_tableau() +
  facet_grid(. ~ Group) +
  theme_pander(base_size = 10) +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    axis.line = element_line()
    )


# DIVERSITY INBRED LINE CORE SET ------------------------------------------
all_geno_pc <- "./data/derived/maizego/maizego_pca.RDS" %>% 
  readRDS()
common_inbreds <- "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS" %>% 
  readRDS() %>% 
  reduce(intersect)
inbreds <- all_geno_pc %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  gather(key = PC, value = Score, -G) %>%
  as_data_frame() %>%
  filter(PC %in% paste0("PC_", seq_len(2))) %>%
  mutate(Cluster = if_else(
    G %in% common_inbreds, true = "Core", false = "Full"),
    PC = gsub("PC_", replacement = "PC", x = PC),
    Group = "Diversity"
  )

g1 <- rbind(inbreds, hybrids) %>%
  spread(key = PC, value = Score) %>%
  mutate(
    Group = factor(Group, levels = c("Diversity", "Dent", "Flint"))
  ) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
  geom_point() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_tableau() +
  facet_wrap(~ Group, scales = "free") +
  theme_pander(base_size = 10) +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    axis.line = element_line()
    )

ggsave(plot = g1, 
       filename = "./paper/tables_figures/diversity_dent_flint_pca.pdf",
       height = 3, 
       width = 6,
       units = "in")
