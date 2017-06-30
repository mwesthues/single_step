## Goal: 
# Sample core sets of genotypes for the Yan-lab data in increments of ten 
# percentage points.

## Topics
# Level 2 sampling
# Level 1 sampling
# ETA setup


## -- PACKAGES AND FUNCTIONS -----------------------------------------------
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("mwesthues/sspredr", update = FALSE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "methods", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# Prediction helper functions (outsourced from this script for improved
# readability).
source("./software/prediction_helper_functions.R")




## -- FUNCTIONS -------------------------------------------------------------
compute_totals_from_fractions <- function(nm, frac) {
  # ---
  # convert a fraction 'frac' * 100% of genotypes, that shall be included in 
  # the training set, into an integer specifying the number of genotypes to be
  # included in the training set.

  # nm: vector of names
  # frac: fractions of genotypes with both, mRNA and SNP information
  #---
  nm %>%
    length() %>%
    map(function(i) i * frac) %>%
    flatten_dbl() %>%
    ceiling()
}





sample_fraction <- function(nm, frac) {
  #---
  # generate a resampled training set provided a vector of genotype names and
  # the fraction of genotypes that shall be included in this training set.

  # nm: names of genotypes
  # frac: vector of genotypes fractions
  #---
  totals <- compute_totals_from_fractions(
    nm = nm,
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





create_t0_trn <- function(hybrids, split_char = "_") {
  #---
  # take a vector of hybrid names ('hybrids'), extract their parental components
  # and, separately for each hybrid, return the corresponding set of hybrids
  # that do not share a parent with the test set hybrid under evaluation

  # hybrids: names of hybrids
  # split_char: separator character for parental components
  #---

  # separate paternal and maternal genotypes and pair them with their offspring
  parent_df <- hybrids %>%
    as_data_frame() %>%
    rename(Hybrid = "value") %>%
    separate(
      Hybrid,
      into = c("Parent_A", "Parent_B"),
      sep = split_char,
      remove = FALSE
    )

  # returns the corresponding t0 hybrids for the test set hybrid under
  # consideration
  select_t0_trn_hybrids <- function(i) {
    parent_vec <- parent_df %>%
      slice(i) %>%
      select(-Hybrid) %>%
      flatten_chr()

    trn_hybrids <- parent_df %>%
      filter(Parent_A != parent_vec[1], Parent_B != parent_vec[2]) %>%
      pull(Hybrid)
    trn_hybrids
  }

  # get the set of t0 hybrid for each test set hybrid
  parent_df %>%
    nrow() %>%
    seq_len() %>%
    map(select_t0_trn_hybrids)
}




sample_loo_sets <- function(geno,
                            frac,
                            iter,
                            material = "Inbred",
                            split_char = "_") {
  #---
  # given a vector of genotype names, that shall be evaluated, generate 'iter'
  # randomly sampled training sets for each genotype, where 'frac' * 100% of the
  # genotypes in 'geno' shall be included in each training set

  # fraction: fraction of genotypes used for the training set
  # geno: names of genotypes
  # iter: number of replicates per genotype
  # material: 'Inbred' or 'Hybrid'
  # split_char: separator for hybrid parents
  #---

  # build the training set by excluding every genotype once from all others
  if (isTRUE(material == "Inbred")) {
    train_geno <- map(geno, .f = ~ setdiff(geno, .))
  } else if (isTRUE(material == "Hybrid")) {
    train_geno <- create_t0_trn(geno, split_char = split_char)
  }
  names(train_geno) <- geno

  # from each training set, sample 'frac' * 100% of genotypes at random and
  # declare them as the resampled training subset. repeat this 'iter' times.
  train_lst <- rerun(.n = iter, {
    sub_train_lst <- map(train_geno, .f = ~sample_fraction(., frac = frac))
    sub_train_df <- sub_train_lst %>% 
      bind_rows(.id = "TST_Geno") %>%
      as_data_frame()
    sub_train <- dplyr::rename(sub_train_df, TRN_Geno = values)
    sub_train
  })

  # concatenate the ten different training sets per test set genotype in one
  # data frame.
  iter_seq <- seq_len(iter)
  names(train_lst) <-  as.character(iter_seq)
  train_df <- bind_rows(train_lst, .id = "Iter")
  
  ### --- Tests ---
  ## test whether the test set genotype is never included in any of the 
  ## training sets
  if (isTRUE(material == "Inbred")) {
    isect <- map2(train_geno, geno, .f = intersect) %>%
      map_int(length) %>%
      unique()
  } else if (isTRUE(material == "Hybrid")) {
    isect <- train_geno %>%
      map(strsplit, split = split_char) %>%
      map(flatten_chr) %>%
      map2(.y = geno, .f = intersect) %>%
      map_int(length) %>%
      unique()
  }
  stopifnot(isect == 0)

  ## test whether the training sets differ in their size, which must be true for
  ## hybrids
  if (isTRUE(material == "Hybrid")) {
    t0_trn_set_size <- train_df %>%
      select(-id) %>%
      group_by(Iter, TST_Geno) %>%
      count() %>%
      pull(n) %>%
      unique() %>%
      length()
      if (isTRUE(t0_trn_set_size == 0)) {
        stop("Hybrid training sets should differ in their sizes")
      }
  }
  ## --- End of tests ---

  # return the data frame with the resampling scheme
  train_df
}






## -- LEVEL2 SAMPLING -------------------------------------------------------
# concatenate all predictor data and phenotypic data in one list
named_list <- create_named_list(
  inb_snp_path = "./data/processed/maizego/imputed_snp_mat.RDS",
  inb_agro_path = "./data/derived/maizego/tst_pheno_tibble.RDS",
  inb_mrna_path = "./data/derived/maizego/mrna.RDS",
  hyb_snp_path = "./data/derived/uhoh/snp_matrix.RDS",
  hyb_agro_path = "./data/derived/uhoh/agro_tibble.RDS",
  hyb_mrna_path = "./data/derived/uhoh/mrna.RDS",
  hyb_ped_path = "./data/derived/uhoh/pedigree_matrix.RDS"
)
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
  select(Path) %>%
  flatten_chr() %>%
  map(readRDS)

pred_lst_names <- named_df %>%
  select(Predictor) %>%
  flatten_chr()
names(pred_lst) <- pred_lst_names

# load genotypes from the first and the fourth cluster of the pca (Scenario 'a')
cluster14 <- readRDS("./data/derived/maizego/cluster_14_genotypes.RDS")



### hybrids
hyb_pheno <- pred_lst %>%
  keep(names(.) == "hyb_agro") %>%
  .[[1]]

full_hyb_genos <- hyb_pheno %>%
  pull(Genotype) %>%
  unique()

core_hyb_parents <- pred_lst %>%
  keep(names(.) == "hyb_mrna") %>%
  .[[1]] %>%
  rownames()

core_hyb_genos <- hyb_pheno %>%
  separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
  filter(Dent %in% core_hyb_parents & Flint %in% core_hyb_parents) %>%
  unite(Hybrid, Dent, Flint, sep = "_") %>%
  mutate(Hybrid = paste0("DF_", Hybrid)) %>%
  select(Hybrid) %>%
  unique() %>%
  flatten_chr()

### inbred lines
full_inb_bc <- pred_lst %>% 
  keep(names(.) == "inb_snp") %>% 
  .[[1]] %>% 
  rownames()

full_inb_a <- intersect(full_inb_bc, cluster14)

core_inb_bc <- pred_lst %>%
  keep(names(.) == "inb_mrna") %>%
  .[[1]] %>%
  rownames()
core_inb_a <- intersect(core_inb_bc, cluster14)

### store all genotype information in one data frame
# note that scenarios 'b' and 'c' are based on the same set of genotypes, and
# thus, will not be considered as separate cases to save memory
geno_lst <- list(
  Full_Hybrid_None = full_hyb_genos,
  Core_Hybrid_None = core_hyb_genos,
  Full_Inbred_A = full_inb_a,
  Full_Inbred_B = full_inb_bc,
  Core_Inbred_A = core_inb_a,
  Core_Inbred_B = core_inb_bc
)

geno_df <- geno_lst %>%
  utils::stack(drop = TRUE) %>%
  separate(ind, into = c("Extent", "Material", "Scenario"), sep = "_") %>%
  rename(G = "values") %>%
  as_data_frame()

# data frame with all combinations that need to be tested
geno_param_df <- geno_df %>%
  select(-G) %>%
  unique()

sample_by_combination <- function(i) {
  param_subset <- slice(geno_param_df, i)
  param_genos <- geno_df %>%
    inner_join(
      y = param_subset,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    pull(G)
  
  loocv_samples <- sample_loo_sets(
  geno = param_genos,
  frac = 0.7,
  iter = 50
  )
  loocv_samples
}

rnd_level2_df <- geno_param_df %>%
  nrow() %>%
  seq_len() %>%
  map(sample_by_combination) %>%
  purrr::set_names(nm = names(geno_lst)) %>%
  bind_rows(.id = "ID") %>%
  separate(ID, into = c("Extent", "Material", "Scenario"), sep = "_") %>%
  rename(Rnd_Level2 = "Iter")
saveRDS(rnd_level2_df, "./data/derived/predictor_subsets/rnd_level2.RDS")








## -- LEVEL 1 SAMPLING ----------------------------------------------------
# level 1 sampling is only required for the core set of inbred lines from the
# maize diversity panel
core_fractions <- seq(from = 0.1, to = 0.9, by = 0.1)
level1_frame <- geno_df %>%
  filter(Extent == "Core", Material == "Inbred") %>%
  pull(Scenario) %>%
  unique()

level1_geno_df <- geno_df %>%
  filter(Extent == "Core", Material == "Inbred") %>%
  select(G, Scenario)

sample_core <- function(x) {
  geno <- pull(x, G)
  smp_fraction <- x %>%
    pull(Core_Fraction) %>%
    unique()
  sample_frac(x, size = smp_fraction)
}

rnd_level1_df <- level1_frame %>%
  expand.grid(
    Extent = "Core",
    Material = "Inbred",
    Scenario = .,
    Core_Fraction = core_fractions,
    Rnd_Level1 = seq_len(50),
    stringsAsFactors = FALSE
  ) %>%
  as_data_frame() %>%
  right_join(., y = level1_geno_df, by = "Scenario") %>%
  split(list(.$Scenario, .$Core_Fraction, .$Rnd_Level1)) %>%
  map(sample_core) %>%
  bind_rows()

# Add all other combinations to the data frame. This will facilitate the setup
# of ETA objects later on.
rnd_level1_df <- geno_df %>%
  filter(!(Extent == "Core" & Material == "Inbred")) %>%
  mutate(Core_Fraction = 0, Rnd_Level1 = 0) %>%
  bind_rows(., rnd_level1_df) %>%
  mutate_all(as.character)
saveRDS(rnd_level1_df, "./data/derived/predictor_subsets/rnd_level1.RDS")




## -- ETA SETUP -----------------------------------------------------------
eta_spec1 <- expand.grid(
  Extent = "Core",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "mrna",
  Scenario = c("A", "B"),
  Rnd_Level1 = seq_len(50) %>% as.character(),
  Core_Fraction = seq(from = 0.1, to = 0.9, by = 0.1) %>% as.character(),
  stringsAsFactors = FALSE
)

eta_spec2 <- expand.grid(
  Extent = "Core",
  Material = "Inbred",
  Pred1 = c("snp", "mrna"),
  Pred2 = "",
  Scenario = c("A", "B"),
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec3 <- expand.grid(
  Extent = "Full",
  Material = "Inbred",
  Pred1 = c("snp", "mrna"),
  Pred2 = "",
  Scenario = "None",
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec4 <- expand.grid(
  Extent = "Full",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "mrna",
  Scenario = "None",
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec5 <- expand.grid(
  Extent = "Core",
  Material = "Hybrid",
  Pred1 = c("snp", "ped", "mrna"),
  Pred2 = "",
  Scenario = "None",
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec6 <- expand.grid(
  Extent = "Full",
  Material = "Hybrid",
  Pred1 = c("snp", "ped"),
  Pred2 = "mrna",
  Scenario = "None",
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec_df <- list(
  eta_spec1,
  eta_spec2,
  eta_spec3,
  eta_spec4,
  eta_spec5,
  eta_spec6
  ) %>%
  bind_rows() %>%
  as_data_frame()


scenario_seq <- eta_spec_df %>%
  nrow() %>%
  seq_len()

lapply(scenario_seq, FUN = function(i) {

  scen_i <- slice(eta_spec_df, i)
  pred1 <- pull(scen_i, Pred1)
  pred2 <- pull(scen_i, Pred2)

  tmp_geno <- inner_join(
    rnd_level1_df,
    scen_i,
    by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
  ) %>%
  pull(G)

  if (isTRUE(nchar(pred2) == 0)) {
    pred1_geno <- tmp_geno
  } else if (all(nchar(c(pred1, pred2)) != 0)) {
    pred2_geno <- tmp_geno
    pred1_geno <- inner_join(
      geno_df,
      scen_i,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    pull(G)
  }
  
  init_pre_model <- "BRR"




}


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


## -- ANALYSIS -------------------------------------------------------------
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
  core_genotypes <- sample_fraction(
    nm = common_genotypes,
    frac = fractions
    ) %>%
    as_data_frame() %>%
    mutate(ind = as.character(ind)) %>%
    filter(ind == core_set) %>%
    pull(values)
  
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

