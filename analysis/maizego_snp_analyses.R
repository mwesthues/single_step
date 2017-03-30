# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "LEA", "stringr", "viridis", "cowplot", 
               "ggthemes")

# For the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- common_genotypes %>%
  reduce(intersect)

snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  .[rownames(.) %in% common_genotypes, ] %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05, 
    any_missing = FALSE, remove_duplicated = TRUE
  )


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
# This code snippet is commented out because it takes a long time to run and 
# can easily and inadvertantly be overridden.
#project <- snmf(input.file = "./data/derived/maizego/snp.geno",
#                K = seq_len(tw_K),
#                project = "new",
#                repetitions = 25,
#                CPU = 3,
#                entropy = TRUE)

# Based on the results of a PCA, decide whether it is even necessary to treat
# subpopulations as such or whether all inbred lines can be treated as a single
# group.
project <- load.snmfProject("./data/derived/maizego/snp.snmfProject")


# SNP VS MRNA PCA ---------------------------------------------------------

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
all_geno_pc %>% 
  saveRDS(., file = "./data/derived/maizego/maizego_pca.RDS")



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
  bind_rows(.id = "K") %>%
  mutate(
    K = as.numeric(K),
    K = K + 1,
    K = paste0("K=", as.character(K))
  ) %>%
  mutate_at(vars(K, G, Component), as.factor)

# For each genotype and for each run with K ancestral populations, plot the 
# admixture coefficients by genotype.
g1 <- structure_df %>%
  ggplot(aes(x = G, y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~K, ncol = 1, strip.position = "right") +
  scale_fill_viridis(discrete = TRUE) +
  theme_pander(base_size = 10) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  xlab("Genotype") +
  ylab("Ancestry Coefficient")

# Core sampling PCA results.
pc_df <- "./data/derived/maizego/core_sampling_pca.RDS" %>%
  readRDS()

g2 <- pc_df %>%
  filter(Fraction != "1.0") %>%
  rename(`Core Group` = Core_Group) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Core Group`, shape = `Core Group`)) +
  geom_point(size = 1) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_tableau() +
  facet_wrap(~Fraction, nrow = 3, ncol = 3) +
  theme_pander(base_size = 10) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.placement = "outside"
  )

ggsave(plot = g2,
       filename = "./paper/tables_figures/maizego_core_sampling_pca.pdf",
       width = 7,
       height = 7,
       units = "in"
       )


g12 <- plot_grid(g1, g2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave(plot = g12,
       filename = "./paper/tables_figures/maizego_admixture_pca.pdf",
       width = 7,
       height = 4,
       units = "in")




# ALTERNATIVE CORE SAMPLING -----------------------------------------------
# Select genotypes whose ancestry coefficients indicate a share of the second 
# ancestral population is equal to or greater than 0.5.
# greater than 
ancestry_df <- structure_df %>%
  filter(K == "K=3") %>%
  droplevels() %>%
  split(.$Component) %>%
  map(droplevels) %>%
  map(., ~filter(., Value >= 0.5)) %>%
  map(., ~select(., G)) %>%
  bind_rows(.id = "Ancestor") %>%
  mutate(
    G = as.character(G),
    Ancestor = gsub(pattern = "V", replacement = "A", x = Ancestor)
  )
ancestry_df %>%
  saveRDS(file = "./data/derived/maizego/ancestry_k3.RDS")
 

g3 <- structure_df %>%
  filter(K == "K=3") %>%
  mutate_at(vars(K, G, Component), .funs = funs(as.character)) %>%
  mutate(Component = gsub(pattern = "V", replacement = "A", x = Component)) %>% 
  ggplot(aes(x = G, y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~K, ncol = 1) +
  scale_fill_manual(values = c("#258039", "#F5BE41", "#31A9B8")) +
  theme_pander(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Genotype") +
  ylab("Ancestry Coefficient")


g4 <- pc_df %>% 
  filter(Fraction == "1.0") %>% 
  full_join(y = ancestry_df, by = "G") %>%
  mutate(Ancestor = if_else(
    is.na(Ancestor), true = "mixed", false = Ancestor
  )) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Ancestor, shape = Ancestor)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#258039", "#F5BE41", "#31A9B8", "red")) +
  scale_shape_manual(values = c(15, 17, 19, 3)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_pander(base_size = 10)


g34_wo_legend <- plot_grid(
  g4, 
  g5 + theme(legend.position = "none"),
  labels = c("A", "B"),
  ncol = 2,
  nrow = 1
)

g34_legend <- get_legend(g4 + theme(legend.position = "top"))
g34 <- plot_grid(g34_legend, g34_wo_legend, ncol = 1, rel_heights = c(0.15, 1))
g34

ggsave(plot = g34,
       filename = "./paper/tables_figures/maizego_admixture_k3_pca.pdf",
       width = 7,
       height = 4,
       units = "in"
       )

