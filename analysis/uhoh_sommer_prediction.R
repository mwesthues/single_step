if (!require("pacman")) install.packages("pacman")
pacman::p_load("sommer", "tidyverse", "data.table", "dtplyr", "BGLR", 
               "parallel", "ggthemes", "stringr", "methods")
devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")


source("./software/prediction_helper_functions.R")

pred_sets <- "mrna"
traits <- "GTM"
use_cores <- 3
runs <- "1-10"




if (nchar(runs) != 0) {
  runs_split_char <- runs %>%
    str_extract(., pattern = "[^[:alnum:]]") 
  
  run_length <- runs %>%
    strsplit(runs, split = runs_split_char) %>%
    .[[1]] %>%
    as.numeric() %>%
    (function(x) seq(from = x[1], to = x[2], by = 1)) 
}


traits_split_char <- traits %>%
  str_extract(., pattern = "[^[:alnum:]]") 

traits <- traits %>%
  strsplit(., split = traits_split_char) %>%
  .[[1]] %>%
  as.character()


# DATA ------------------------------------------------------------------
snp_path <- "./data/derived/uhoh/snp_matrix.RDS"
agro_path <- "./data/derived/uhoh/agro_tibble.RDS"
mrna_path <- "./data/derived/uhoh/mrna.RDS"
ped_path <- "./data/derived/uhoh/pedigree_matrix.RDS"
named_list <- create_named_list(snp_path, agro_path, mrna_path, ped_path)


named_df <- data.frame(Predictor = names(named_list),
                       Path = unlist(named_list),
                       stringsAsFactor = FALSE)
named_df$Predictor <- gsub(named_df$Predictor,
                           pattern = "_path",
                           replacement = "")
named_df <- named_df %>%
  as_data_frame() %>%
  mutate(Path = as.character(Path))

pred_lst <- named_df %>%
  as_data_frame() %>%
  filter(Predictor %in% pred_sets) %>%
  select(Path) %>%
  flatten_chr() %>%
  map(readRDS)

# The order in the path data frame may not coincide with the expected predictor 
# order. Therefore, we'll extract the predictor names - in the correct order - 
# from the path data frame and assign them to the loaded predictor data.
pred_lst_names <- named_df %>%
  filter(Predictor %in% pred_sets) %>%
  select(Predictor) %>%
  flatten_chr()
names(pred_lst) <- pred_lst_names


## -- PREDICTOR TRANSFORMATION FOR HYBRID DATA ---------------------------
pheno <- named_df %>%
  filter(Predictor == "agro") %>%
  select(Path) %>%
  flatten_chr() %>%
  readRDS()

high_coverage_geno_name_id <- pred_lst %>%
  map(rownames) %>%
  map(length) %>%
  flatten_int() %>%
  which.max()
pred_genotypes <- pred_lst %>%
  .[[get("high_coverage_geno_name_id")]]  %>%
  rownames() 


# Get the names of genotypes for which there are (parental) inbred lines with 
# information on at least one predictor.
# If we do not do this, the BGLR-input matrices will be augmented to genotypes
# for which we cannot make predictions and the algorithm will fail.
covered_genotypes <- pheno %>%
  separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
  filter(Dent %in% pred_genotypes & Flint %in% pred_genotypes) %>%
  unite(Hybrid, Dent, Flint, sep = "_") %>%
  mutate(Hybrid = paste0("DF_", Hybrid)) %>%
  select(Hybrid) %>%
  unique() %>%
  flatten_chr()

pheno <- pheno %>%
  filter(Genotype %in% covered_genotypes)

hybrid_parents <- pheno %>%
  separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
  select(Dent, Flint) %>%
  unique() %>%
  gather(key = Pool, value = Genotype) %>%
  split(.$Pool) %>%
  map(., ~select(., Genotype)) %>%
  map(flatten_chr)

hybrid_names <- pheno %>%
  select(Genotype) %>%
  unique() %>%
  flatten_chr()


pheno <- pheno %>%
  filter(Trait %in% traits,
         Genotype %in% hybrid_names) %>%
  spread(key = Trait, value = Value) %>%
  rename(Hybrid = Genotype) %>%
  mutate(Hybrid = gsub("DF_", replacement = "", x = Hybrid)) %>%
  separate(col = Hybrid, into = c("Dent", "Flint"),
           sep = "_", remove = FALSE) %>%
  mutate(Dent = as.factor(Dent),
         Flint = as.factor(Flint))



# In the case of hybrid data, we still need to split all predictor matrices
# into Flint and Dent components first.
# 1. map_if
# In the case of pedigree data, we need to split the data by genotypes in the
# x- as well as the y-dimension because there are no features, which is
# different for other predictor matrices.
# 2. map_at("snp")
# Ensure that SNP quality checks are applied separately to the genomic data 
# for each heterotic group in order to have only polymorphic markers and no 
# markers in perfect LD.
Z_lst <- pheno %>%
  gather(key = Group, value = ID, -Hybrid) %>%
  select(-Hybrid) %>%
  split(.$Group) %>%
  map(~mutate(., ID = as.factor(ID))) %>%
  map(., ~model.matrix(~ ID - 1, data = .)) %>%
  map(function(x) {
    colnames(x) <- gsub(pattern = "ID", replacement = "", x = colnames(x))
    x
  })

# Separately for each heterotic group, apply SNP quality checks, then generate 
# kernels from the marker data.
M_lst <- pred_lst %>%
  map(.f = split_into_hetgroups, y = hybrid_parents, pedigree = FALSE) %>%
  map_at("snp", .f = ~map(., ~sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE,
    maf_threshold = 0.05, any_missing = FALSE, remove_duplicated = TRUE
    )))

G_lst <- M_lst %>%
  at_depth(.depth = 2, function(M) {
    W <- scale(M, center = TRUE, scale = TRUE)
    G <- tcrossprod(W) / ncol(W)
    diag(G) <- diag(G) + 0.01
    G
  }) 

G_dent <- G_lst[[1]][["Dent"]][levels(pheno$Dent), levels(pheno$Dent)]
G_flint <- G_lst[[1]][["Flint"]][levels(pheno$Flint), levels(pheno$Flint)]


# Define which runs shall be used, i.e., which genotypes shall be included as
# test set members.
if (nchar(runs) == 0) {
  run_length <- seq_along(hybrid_names)
}
param_df <- expand.grid(Trait = traits,
                        Run = run_length,
                        stringsAsFactors = FALSE)



# Compute predictive abilities.
# Repeat the process once so that we can compare the uncertainty within
# BGLR-estimates and within sommer-estimates, respectively.
# Compute the predictive ability for each cross-validation run.
predicted_data <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {

  run <- param_df[i, "Run"]
  trait <- param_df[i, "Trait"]

  geno <- pheno %>%
    slice(i)

  t0_training_set <- pheno %>%
    filter(!Dent %in% (geno %>% .$Dent %>% as.character()),
           !Flint %in% (geno %>% .$Flint %>% as.character())) %>%
    select(Hybrid) %>%
    flatten_chr()

  # Set test-set hybrids to NA so that they will be estimated using training
  # data.
  observed_trait <- pheno %>%
    gather(key = Trait, value = Measurement, -Hybrid, -Dent, -Flint) %>%
    filter(Trait == trait) %>%
    mutate(Measurement = ifelse(Hybrid %in% t0_training_set, 
                                yes = Measurement, no = NA_real_)) %>%
    select(Measurement) %>%
    flatten_dbl()

  ### For the 'sommer' package, no Cholesky decomposition is required in 
  ### advance.
  # For the 'BGLR' package, matrix factorization is necessary in advance.
  ETA <- list(GCA_Dent = list(Z = Z_lst$Dent, K = G_dent),
              GCA_Flint = list(Z = Z_lst$Flint, K = G_flint))

  r <- mmer(
    Y = observed_trait,
    Z = ETA,
    method = "EM",
    silent = TRUE,
    draw = FALSE
    ) %>%
    .$fitted.y %>%
    slot(., "x") %>%
    .[run]
    
  new_column <- paste0("Predicted_", trait)
  geno[, new_column] <- r
  geno
}, mc.cores = use_cores) %>%
  rbindlist()
