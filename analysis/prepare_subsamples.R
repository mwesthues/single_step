# Goal: Sample core sets of genotypes for the Yan-lab data and scenario two of
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("mwesthues/sspredr", update = FALSE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "methods", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# Prediction helper functions (outsourced from this script for improved
# readability).
source("./software/prediction_helper_functions.R")




# load genotypes from the first and the fourth cluster of the pca
cluster14 <- readRDS("./data/derived/maizego/cluster_14_genotypes.RDS")

# for the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic) and which
# belong to either the first or the fourth cluster of the pca
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes[["cluster14"]] <- cluster14
common_genotypes <- common_genotypes %>%
  reduce(intersect)

## --- Scenario 2 SNP preparation.
# Explore the kinship matrix.
genos <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  rownames()

runs <- 10
fractions = seq(from = 0.1, to = 0.9, by = 0.1)



#### Functions
compute_totals_from_fractions <- function(nm, frac) {
  # nm: vector of names
  # frac: fractions of genotypes with both, mRNA and SNP information
  totals <- nm %>%
    length() %>%
    map(function(i) i * frac) %>%
    flatten_dbl() %>%
    ceiling()
}


sample_fraction <- function(nm, frac) {
  # nm: names of genotypes
  # frac: vector of genotypes fractions
  totals <- compute_totals_from_fractions(
    nm = common_genotypes,
    frac = frac
  )
  totals %>%
    map(.f = ~sample(
      x = nm,
      size = .,
      replace = FALSE
    )) %>%
  set_names(., nm = as.character(frac)) %>%
  stack()
}


sample_loo <- function(geno, frac, iter) {
  # fraction: fraction of genotypes used for the training set
  # geno: names of genotypes
  # iter: number of replicates per genotype
  train_geno <- map(geno, .f = ~ setdiff(geno, .))
  names(train_geno) <- geno
  train_lst <- rerun(.n = iter, {
    sub_train_lst <- map(train_geno, .f = ~sample_fraction(., frac = frac))
    sub_train_df <- bind_rows(sub_train_lst, .id = "TST_Geno")
    sub_train <- dplyr::rename(sub_train_df, TRN_Geno = values)
  })
  iter_seq <- seq_len(iter)
  names(train_lst) <-  as.character(iter_seq)
  train_df <- bind_rows(train_lst, .id = "Iter")
  as_data_frame(train_df)
}
###



# Create random LOOCV training set samples for each genotype.
set.seed(9434)
loocv_samples <- sample_loo(
  geno = common_genotypes,
  frac = 0.7,
  iter = runs
)
saveRDS(loocv_samples, "./data/derived/predictor_subsets/loocv_samples.RDS")





# -- PARAMETERS -------------------------------------
init_pred_combi <- c("snp", "mrna", "snp_mrna")
init_pred_model <- "BRR"
init_core_set <- seq(from = 0.1, to = 1, by = 0.1)
init_core_nmb <- seq(from = 1, to = runs, by = 1)

# Combination of all parameters
param_df <- expand.grid(
  Predictor = init_pred_combi,
  Model = init_pred_model,
  Core_Set = init_core_set,
  Random_Sample = init_core_nmb,
  stringsAsFactors = FALSE
  ) %>%
  filter(!(Predictor == "snp_mrna" & Core_Set == 1))


log_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {

  pred_combi <- param_df %>%
    slice(i) %>%
    pull(Predictor)
  pred_model <- param_df %>%
    slice(i) %>%
    pull(Model)
  core_set <- param_df %>%
    slice(i) %>%
    pull(Core_Set)
  core_nmb <- param_df %>%
    slice(i) %>%
    pull(Random_Sample)

  ## -- PREDICTOR INPUT SPECIFICATION -------------------------------------
  pred_sets <- pred_combi %>%
    strsplit(split = "_") %>%
    flatten_chr()
  
  
  
  ## -- DATA ----------------------------------------------------------------
  snp_path <- "./data/processed/maizego/imputed_snp_mat.RDS"
  agro_path <- "./data/derived/maizego/tst_pheno_tibble.RDS"
  mrna_path <- "./data/derived/maizego/mrna.RDS"
  named_list <- create_named_list(snp_path, agro_path, mrna_path)
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
  
  # Make sure that the predictor matrices are in the intended order, which is
  # crucial for the imputation of the predictor matrix that covers fewer
  # genotypes.
  pred_lst <- pred_lst[match(names(pred_lst), pred_sets)]
  
  # Get the random core set of genotypes covering the genetic target space
  # of a pre-specified size (as a fraction of genotypes covered by mRNA data).
  set.seed(9434)
  core_genotypes <- rerun(.n = runs, sample_fraction(
    nm = common_genotypes,
    frac = fractions
    )) %>%
    bind_rows(.id = "Rep") %>%
    as_data_frame() %>%
    mutate(ind = as.character(ind)) %>%
    filter(ind == core_set, Rep == core_nmb) %>%
    pull(values)
  
  # Extract the intersect between genotypes that have data for genomic as well
  # as transcriptomic features.
  common_genotypes <- readRDS(
    "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
  )
  common_genotypes <- common_genotypes %>%
    reduce(intersect)
  
    
  # Reduce the predictor data to match the size of the pre-specified core sets.
  # Ensure that SNP quality checks are applied to the genomic data in order to 
  # have only polymorphic markers and no markers in perfect LD.
  pred_lst <- pred_lst %>%
    modify_at("mrna", .f = function(x) {
      x[rownames(x) %in% core_genotypes, ]
    }) %>%
    modify_at("snp", .f = function(x) {
      x[rownames(x) %in% common_genotypes, ]
    }) %>%
    modify_at("snp",
           .f = ~sspredr::ensure_snp_quality(
             ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
             any_missing = FALSE, remove_duplicated = TRUE
             )
    )
  
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
  covered_genotypes <- pred_genotypes
  saveRDS(
    object = covered_genotypes,
    file = paste0(
      "./data/derived/predictor_subsets/covered_genotypes_",
      i,
      ".RDS"
    )
  )
   
  
  
  pre_eta <- list(Inbred = pred_lst)
  pre_eta[[1]][["geno"]] <- pred_genotypes
  
  
  raw_args <- pre_eta %>%
    map(names) %>%
    reduce(intersect) %>%
    discard(. == "geno")
  pred_number <- raw_args %>%
    .[grep("ped|snp|mrna", x = .)] %>%
    length()
  
  
  # If all predictors are being used, it is known that pedigree data are involved
  # and no measures need to be taken.
  # If only a subset of predigrees will be used, it needs to be determined
  # whether pedigree data are involved or not, since subsequent functions need to
  # 'know' if the predictor matrices need to be scaled or not.
  if (all(grepl("ped", x = raw_args))) {
    pre_eta <- pre_eta %>%
      map(.f = function(x) {
        x$as_kernel <- FALSE
        x$bglr_model <- pred_model
        x
    })
  } else if (any(grepl("ped", x = raw_args))) {
    pre_eta <- pre_eta %>%
      map(.f = function(x) {
        x$is_pedigree <- TRUE
        x$bglr_model <- pred_model
        x$as_kernel <- TRUE
        x
      })
  } else {
    pre_eta <- pre_eta %>%
      map(.f = function(x) {
        x$is_pedigree <- FALSE
        x$bglr_model <- pred_model
        x$as_kernel <- TRUE
        x
      })
  }
  
  # Depending on the number of predictors included in the set, choose the correct
  # function for the set-up of the BGLR kernels.
  eta_fun <- c("complete_eta", "impute_eta") %>%
    .[pred_number]
  
  # Collect the formal arguments from the selected function and then rename the
  # corresponding predictors in the 'pre_eta' object so that the call of the 
  # kernel building function (e.g. "impute_eta") can be generalized.
  eta_fun_args <- eta_fun %>%
    formals() %>%
    names()
  eta_pred_nms <- eta_fun_args %>%
    .[!grepl("as_kernel|is_pedigree|geno|bglr_model", x = .)]
  current_eta_names <- pre_eta %>%
    map(names) %>%
    reduce(intersect)
  current_eta_names[seq_len(pred_number)] <- eta_pred_nms
  pre_eta <- pre_eta %>%
    map(function(x) {
      names(x) <- current_eta_names
      x
    })
  # Build the kernel for BGLR.
  eta <- invoke_map(get(eta_fun), pre_eta) %>%
    flatten()

  ## Save individual files.
  data_location <- paste0("./data/derived/predictor_subsets/eta_", i, ".RDS")
  saveRDS(
    object = eta,
    file = data_location
  )

  # Generate a log file to be able to retrieve all information later on.
  data_frame(
    Data_Location = data_location,
    Predictor = pred_combi,
    Model = pred_model,
    Core_Set = core_set,
    Random_Sample = core_nmb
  )

}, mc.cores = 3, mc.preschedule = FALSE)

# Safe the log file.
log_idx <- log_lst %>%
  map_lgl(., ~all(class(.) != "try-error"))
log_df <- log_lst %>%
  .[log_idx] %>%
  bind_rows()
log_location <- "./data/derived/log_prepare_subsamples.txt"
write.table(
  log_df,
  file = log_location,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  append = FALSE
)

