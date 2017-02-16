
# Data and packages -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes", "dplyr", "forcats", 
               "viridis", "tidyr", "knitr", "tibble", "purrr")

pheno <- readRDS("./data/derived/maizego/tst_pheno_tibble.RDS")
# Load the names of unique genotypes for each data type (phenotypic, agronomic,
# transcriptomic).
geno_lst <- readRDS("./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS")

# Specify which data are available for each inbred line.
pheno <- pheno %>%
  mutate(Trait = as.factor(Trait),
         Reduced_Data = "no",
         Reduced_Data = ifelse(Genotype %in% geno_lst$mrna, yes = "yes",
                               no = "no"))

reduced_pheno <- pheno %>%
  filter(Reduced_Data == "yes")



# Data distribution -------------------------------------------------------
# Histogram of the distribution of seven agronomic traits.
# Black bins show the distribution of the reduced set of hybrids whereas light 
# grey bins show the distribution of the full set of hybrids.
ggplot(pheno, aes(x = Value)) +
  geom_histogram(fill = "grey", alpha = 0.5) +
  geom_histogram(data = reduced_pheno, color = "black") +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw()


# Skewness ----------------------------------------------------------------
# Generate 'iter' vectors of length 'number of genotypes with agronomic data', 
# sampled at random from a normal distribution with mean 0 and variance 1.
# Then compare observed skewness for traits of interest to this simulated 
# distribution to decide whether the observed skewness is significantly 
# different from the expected skewness.
iter <- 1e+05
set.seed(3409213)
random_skewness <- iter %>%
  replicate(., rnorm(geno_lst %>% .$pheno %>% length())) %>%
  as_data_frame() %>%
  gather(key = Run, value = Value) %>%
  group_by(Run) %>%
  summarize(Skewness = e1071::skewness(Value),
            Skewness = round(Skewness, digits = 2)) %>%
  select(Skewness) %>%
  flatten_dbl()
hist(random_skewness)
 
compare_skewness <- function(observed, simulated) {
  mean(abs(observed) <= abs(simulated))
}

all_skewness <- pheno %>%
  group_by(Trait) %>%
  summarize(Skewness = e1071::skewness(Value),
            Skewness = round(Skewness, digits = 2),
            p_value = compare_skewness(observed = Skewness,
                                       simulated = random_skewness))

# Correlations ------------------------------------------------------------
# Pairwise correlations among seven agronomic traits for the full set of hybrids 
# (upper diagonal) and the reduced set of hybrids (lower diagonal), 
# respectively.
cor_coeffs <- pheno %>%
  spread(key = Trait, value = Value) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "Genotype") %>%
  select(-Reduced_Data) %>%
  cor() %>%
  .[upper.tri(., diag = FALSE)]
traits <- levels(pheno$Trait)
cor_mat <- diag(nrow = length(traits), ncol = length(traits))
cor_mat[upper.tri(cor_mat, diag = FALSE)] <- cor_coeffs
cor_mat[lower.tri(cor_mat, diag = FALSE)] <- cor_coeffs
dimnames(cor_mat) <- list(traits, traits)
cor_mat %>%
  corrplot::corrplot(col = viridis(256),
                     method = "color",
                     diag = FALSE,
                     tl.srt = 60,
                     cl.cex = 0.8,
                     number.cex = 0.8,
                     order = "hclust",
                     addCoef.col = "black",
                     mar = c(2, 0, 0, 0))
