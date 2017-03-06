if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "caret")

mrna <- readRDS("./data/derived/maizego/tst_raw_mrna_datatable.RDS")
# Load the names of unique genotypes for each data type (phenotypic, agronomic,
# transcriptomic).
geno_lst <- readRDS("./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS")


iter <- 1e+05
set.seed(3409213)
random_skewness <- iter %>%
  replicate(., rnorm(geno_lst %>% .$mrna %>% length())) %>%
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

all_skewness <- mrna %>%
  as_data_frame() %>%
  mutate(Gene_ID = as.factor(Gene_ID)) %>%
  group_by(Gene_ID) %>%
  summarize(Skewness = e1071::skewness(Expression),
            Skewness = round(Skewness, digits = 2),
            p_value = compare_skewness(observed = Skewness,
                                       simulated = random_skewness))
all_skewness %>%
  filter(p_value <= 0.05) %>%
  count()


## Data pre-processing
# Pre-process the data in the following order:
# 
# 1.    Near-zero-variance filtering.
# 2.    Box-Cox tranformation
# 3.    Centering
# 4.    Scaling
mrna %>% 
  as_data_frame() %>%
  mutate(Gene_ID = as.factor(Gene_ID)) %>%
  spread(key = Gene_ID, value = Expression) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Genotype") %>%
  as.matrix() %>%
  scale(., center = TRUE, scale = TRUE) %>%
  saveRDS(., file = "./data/derived/maizego/mrna.RDS")
