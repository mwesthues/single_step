# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis", "cowplot",
               "ggthemes", "pophelper", "forcats")
#devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")

# For the second scenario, keep only names of genotypes, which are covered by
# all data types (i.e. phenotypic, genotypic and transcriptomic).
genotype_sets <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- genotype_sets %>%
  reduce(intersect)

# keep the largest set of genotypes for scenario 'a'
max_geno_number <- genotype_sets %>%
  map_int(length) %>%
  max()
geno_a <- genotype_sets %>%
  keep(~ length(.x) == max_geno_number)
if (length(geno_a) != 1) {
  geno_a <- reduce(geno_a, intersect)
}


snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )


# PCA ---------------------------------------------------------------------
write.lfmm(snp, "./data/derived/maizego/snp.lfmm")

# Determine the structure of the data using genotypic data.
snp_pc <- LEA::pca(input.file = "./data/derived/maizego/snp.lfmm",
                   scale = TRUE)

pc_mat <- snp_pc$projections
rownames(pc_mat) <- rownames(snp)
colnames(pc_mat) <- paste0("PC_", seq_len(ncol(pc_mat)))


# Explained variance of the first five principal components
# Assumed number of ancestral populations: 4
k <- 4
snp_pc %>%
  tracy.widom() %>%
  pull(percentage) %>%
  .[seq_len(k)]


pc_df <- pc_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  gather(key = PC, value = Score, -G) %>%
  as_data_frame() %>%
  filter(PC %in% paste0("PC_", seq_len(2))) %>%
  mutate(
    PC = gsub("PC_", replacement = "PC", x = PC)
  ) %>%
  spread(key = PC, value = Score)



# PCA ---------------------------------------------------------------------
genos <- read_tsv(
  "./data/processed/maizego/popstruc.names",
  col_names = FALSE
  ) %>%
  select(2) %>%
  flatten_chr()




# STRUCTURE FILES ---------------------------------------------------------
slist <- list.files(
  path = "./data/processed/maizego",
  pattern = "popstruc[0-9]",
  full.names = TRUE
  ) %>%
  readQ(
    files = .,
    indlabfromfile = TRUE,
    filetype = "structure"
  )
all_geno_pc <- all_geno_pc$projections
rownames(all_geno_pc) <- rownames(all_geno_snp)
colnames(all_geno_pc) <- paste0("PC_", seq <- len(ncol(all_geno_pc)))
all_geno_pc %>%
  saveRDS(., file = "./data/derived/maizego/maizego_pca.RDS")

# Select K=4
# Assign a color to each cluster.
# Plot the PCA results with colors corresponding to the cluster with the largest
# coefficient for that individual.
get_max_cluster <- function(x, genos) {
  #---
  # for each individual, return the cluster with the largest coefficient.
  #
  # x: list with STRUCTURE objects
  # genos: vector with genotype names pertaining to the STRUCTURE objects
  #---
  x %>%
    mutate(G = genos) %>%
    gather(key = "Cluster", value = "AncCoef", -G) %>%
    group_by(G) %>%
    summarize(main_cluster = which.max(AncCoef)) %>%
    mutate(main_cluster = main_cluster %>% as.character() %>% as.factor()) %>%
    ungroup()
}

cluster_df <- slist %>%
  map(get_max_cluster, genos = genos) %>%
  map(.f = ~right_join(., y = pc_df, by = "G")) %>%
  bind_rows(.id = "K")



change_to_k <- function(x) {
  x %>%
    gsub(x = ., replacement = "", pattern = "_f") %>%
    gsub(x = ., replacement = "", pattern = "popstruc")
}


g1 <- slist %>%
  set_names(., change_to_k(names(.))) %>%
  map(.f = ~mutate(., G = genos)) %>%
  keep(names(.) == k) %>%
  map(.f = ~mutate(., Central_Cluster = Cluster1 + Cluster4)) %>%
  map(.f = ~gather(
    ., key = "Cluster", value = "AncCoef", -G, -Central_Cluster)
    ) %>%
  bind_rows(.id = "K") %>%
  mutate(G = as.factor(G)) %>%
  ggplot(aes(
    x = fct_reorder(G, x = Central_Cluster), y = AncCoef, fill = Cluster)
    ) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )

g2 <- cluster_df %>%
  mutate(K = change_to_k(K)) %>%
  filter(K == k) %>%
  ggplot(aes(x = PC1, y = PC2, color = main_cluster)) +
  geom_point(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_brewer(palette = "Set1") +
  theme_pander(base_size = 10) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.placement = "outside"
  )

plot_grid(
  g1, g2, labels = c("A", "B"), ncol = 2
)


# Select genotypes from the first and the fourth cluster, respectively, given
# that they cluster together and look homogeneous.
selected_geno <- slist %>%
  set_names(., change_to_k(names(.))) %>%
  map(.f = ~mutate(., G = genos)) %>%
  keep(names(.) == k) %>%
  map(get_max_cluster, genos = genos) %>%
  .[[1]] %>%
  filter(main_cluster %in% c("1", "4")) %>%
  pull(G) %>%
  unique()
saveRDS(selected_geno, "./data/derived/maizego/cluster_14_genotypes.RDS")



sel_snps <- snp %>%
  .[rownames(snp) %in% selected_geno, ] %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )
write.lfmm(sel_snps, "./data/derived/maizego/sel_snps.lfmm")
# Determine the structure of the data using genotypic data.
sel_snp_pc <- LEA::pca(
  input.file = "./data/derived/maizego/sel_snps.lfmm",
  scale = TRUE
  )

sel_pc_mat <- sel_snp_pc$projections
rownames(sel_pc_mat) <- rownames(sel_snps)
colnames(sel_pc_mat) <- paste0("PC_", seq_len(ncol(sel_pc_mat)))

sel_pc_df <- sel_pc_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  gather(key = PC, value = Score, -G) %>%
  as_data_frame() %>%
  filter(PC %in% paste0("PC_", seq_len(2))) %>%
  mutate(
    PC = gsub("PC_", replacement = "PC", x = PC)
  ) %>%
  spread(key = PC, value = Score)

sel_snp_pc %>%
  tracy.widom() %>%
  pull(percentage) %>%
  .[seq_len(k)]


