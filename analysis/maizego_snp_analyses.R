# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis", "cowplot",
               "ggthemes", "pophelper", "forcats", "svglite")
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
pc_var_explained <- snp_pc %>%
  tracy.widom() %>%
  pull(percentage)


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

cluster_df %>%
  dplyr::filter(K == "popstruc4_f") %>%
  saveRDS(
    object = .,
    file = "./data/derived/maizego/cluster_df_4pcs.RDS"
  )


change_to_k <- function(x) {
  x %>%
    gsub(x = ., replacement = "", pattern = "_f") %>%
    gsub(x = ., replacement = "", pattern = "popstruc")
}

# Label PC plot components.
label_pca <- function(explained_var, component_number) {
  rounded_var <- round(explained_var, digits = 2)
  paste0(
    "PC", component_number, " (", rounded_var * 100, "%", ")"
    )
}

g1 <- cluster_df %>%
  mutate(K = change_to_k(K)) %>%
  filter(K == k) %>%
  ggplot(aes(x = PC1, y = PC2, color = main_cluster)) +
  geom_point(size = 1) +
  guides(color = guide_legend(
    override.aes = list(size = 3),
    title = "Cluster",
    title.theme = element_text(size = 12, angle = 0)
    )
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_pander(base_size = 12) +
  theme(
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  xlab(label_pca(pc_var_explained[1], 1)) +
  ylab(label_pca(pc_var_explained[2], 2))

ggsave(
  filename = "./paper/tables_figures/inbred_pca.svg",
  plot = g1,
  device = "svg",
  width = 4,
  height = 4,
  units = "in"
  )

