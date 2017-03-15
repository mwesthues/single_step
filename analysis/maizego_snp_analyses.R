# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis", "cowplot", 
               "ggthemes")


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


# SNP VS MRNA PCA ---------------------------------------------------------
# For the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- common_genotypes %>%
  reduce(intersect)

# Explore the kinship matrix.
all_geno_snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )

write.lfmm(all_geno_snp, "./data/derived/maizego/all_tst_snp.lfmm")
all_geno_pc <- LEA::pca("./data/derived/maizego/all_tst_snp.lfmm", scale = TRUE)

all_geno_pc <- all_geno_pc$projections
rownames(all_geno_pc) <- rownames(all_geno_snp)
colnames(all_geno_pc) <- paste0("PC_", seq_len(ncol(all_geno_pc)))

g1 <- all_geno_pc %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  gather(key = PC, value = Score, -G) %>%
  as_data_frame() %>%
  filter(PC %in% paste0("PC_", seq_len(3))) %>%
  mutate(Group = if_else(
    G %in% common_genotypes, true = "mRNA", false = "All"),
    PC = gsub("PC_", replacement = "PC", x = PC)
  ) %>%
  spread(key = PC, value = Score) %>%
  ggplot(aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point() +
  scale_color_tableau() +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
ggsave(plot = g1,
       filename = "./paper/tables_figures/maizego_pca.pdf",
       width = 6,
       height = 4, 
       units = "in")




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
    x$G <- rownames(snp)
    x
  }) %>%
  map(., ~gather(., key = Component, value = Value, -G)) %>%
  bind_rows(.id = "K")

# For each genotype and for each run with K ancestral populations, plot the 
# admixture coefficients by genotype.
g2 <- structure_df %>%
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
  ggthemes::theme_pander(base_size = 10) +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  xlab("Genotype") +
  ylab("Ancestry Coefficient")

# Core sampling PCA results.
pc_df <- "./data/derived/maizego/core_sampling_pca.RDS" %>%
  readRDS()

g3 <- pc_df %>%
  filter(Fraction != "1.0") %>%
  rename(`Core Group` = Core_Group) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Core Group`, shape = `Core Group`)) +
  geom_point(size = 0.5) +
  ggthemes::scale_color_tableau() +
  facet_wrap(~Fraction, nrow = 3, ncol = 3) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")


g23 <- plot_grid(g2, g3, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave(plot = g23,
       filename = "./paper/tables_figures/maizego_admixture_pca.pdf",
       width = 7,
       height = 4,
       units = "in")
