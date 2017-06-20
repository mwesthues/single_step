# Single-step prediction of hybrid performance based on Fernando et al.: A
# class of Bayesian methods to combine large numbers of genotyped and
# non-genotyped animals for whole-genome analyses. Genetics Selection Evolution
# 2014 46:50
#
# The approach will be adjusted in a sense that we do have complete genomic
# instead of complete pedigree data (i.e. A -> G) and that we have incomplete
# transcriptomic instead of incomplete genomic data (i.e. G -> T).
#
#
#
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("BGLR","data.table", "parallel", "magrittr", "dplyr", "dtplyr",
               "tidyr", "purrr", "methods", "tibble", "lubridate")
#devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")

# Load function to save session info on the software used for the analysis in
# this script.
source("./software/session_info.R")
# Prediction helper functions (outsourced from this script for improved
# readability).
source("./software/prediction_helper_functions.R")



## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "3")
  # Are you using inbred ("Inbred") data from the Yan-lab or hybrid ("Hybrid")
  # data from the UHOH-group?
  Sys.setenv("DATA_TYPE" = "Inbred")
  # Which agronomic trait do you want to evaluate? Leave blank if you want to
  # analyze all traits.
  Sys.setenv("TRAIT" = "100grainweight")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "1000") 
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "BRR")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # Main predictor. If 'Pred2' and 'Pred3' are empty, no imputation will take
  # place.
  Sys.setenv("PRED1" = "snp")
  # If 'Pred3' is empty, 'Pred2' will be imputed via information from 'Pred1'.
  Sys.setenv("PRED2" = "mrna")
  # Fraction of genotypes to be included in the core set.
  Sys.setenv("CORE_SET" = "0.9")
  # Which random core sample should be used (integer)
  Sys.setenv("RANDOM_SAMPLE" = "1")
  # Number of genotypes to predict (only for testing!)
  Sys.setenv("RUNS" = "1-2")
  # Output directory for temporary BGLR files
  Sys.setenv("TMP" = "./tmp")
}

if (!interactive()) {
  job_id <- Sys.getenv("MOAB_JOBID")
} else job_id <- "interactive_00"

# Save session info
write_session_info(directory = "./data/derived/session_info/",
                   job_id = job_id)

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
data_type <- as.character(Sys.getenv("DATA_TYPE"))
traits <- as.character(Sys.getenv("TRAIT"))
iter <- as.integer(Sys.getenv("ITER"))
pred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
pred1 <- as.character(Sys.getenv("PRED1"))
pred2 <- as.character(Sys.getenv("PRED2"))
core_set <- as.character(Sys.getenv("CORE_SET"))
core_nmb <- as.character(Sys.getenv("RANDOM_SAMPLE"))
runs <- as.character(Sys.getenv("RUNS"))
temp <- paste0(as.character(Sys.getenv("TMP")), "/")




# -- PREDICTOR INPUT SPECIFICATION -------------------------------------
pred_combi <- list(pred1, pred2) %>%
  keep(nchar(.) != 0) %>%
  flatten_chr() %>%
  paste(., collapse = "_")

pred_sets <- pred_combi %>%
  strsplit(split = "_") %>%
  flatten_chr()


## -- INPUT CHECKS AND UPDATES ------------------------------------------
if (isTRUE(core_set == "1.0" && length(pred_sets) > 1)) {
  stop("Imputation not possible if all predictors cover the same genotypes")
}

if (isTRUE(nchar(core_set) != 0 &&
           core_set != "1.0" &&
           data_type == "Hybrid")) {
  stop("Core sampling for hybrids is not yet supported")
}
  
if (isTRUE(data_type == "Hybrid")) {
  hybrid <- TRUE
} else hybrid <- FALSE


stopifnot(data_type %in% c("Inbred", "Hybrid"))
if (data_type == "Inbred") {
  if (any(grepl(pred_sets, pattern = "ped"))) {
    stop("No pedigree data for inbred lines")
  }
}

if (!all(pred_sets %in% c("ped", "snp", "mrna"))) {
  stop("Predictors are unknown")
}



## -- LOAD PREDICTOR DATA ------------------------------------------------
if (data_type == "Inbred") {
  snp_path <- "./data/processed/maizego/imputed_snp_mat.RDS"
  agro_path <- "./data/derived/maizego/tst_pheno_tibble.RDS"
  mrna_path <- "./data/derived/maizego/mrna.RDS"
  named_list <- create_named_list(snp_path, agro_path, mrna_path)

} else if (data_type == "Hybrid") {
  snp_path <- "./data/derived/uhoh/snp_matrix.RDS"
  agro_path <- "./data/derived/uhoh/agro_tibble.RDS"
  mrna_path <- "./data/derived/uhoh/mrna.RDS"
  ped_path <- "./data/derived/uhoh/pedigree_matrix.RDS"
  named_list <- create_named_list(snp_path, agro_path, mrna_path, ped_path)
}
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

# If we want to evaluate the predictive ability of genomic data that are being
# imputed via pedigree data, we reduce the set of genotypes covered by SNP data
# to the set of genotypes covered by mrna data. 
# This way, we have a common reference (set of genotypes covered by mrna data)
# for all comparisons.
if (isTRUE(data_type == "Hybrid" &&
           all(c("ped", "snp") %in% pred_sets) ||
           core_set == "1.0")) {
  mrna_genotypes <- readRDS(mrna_path) %>%
    rownames()
  pred_lst <- pred_lst %>%
    map_at("snp", .f = function(x) {
      x[match(mrna_genotypes, rownames(x)), ]
    }) %>%
    map_at("ped", .f = function(x) {
      x[match(mrna_genotypes, rownames(x)),
        match(mrna_genotypes, colnames(x))]
    })
}

# Make sure that the predictor matrices are in the intended order, which is
# crucial for the imputation of the predictor matrix that covers fewer
# genotypes.
pred_lst <- pred_lst[match(names(pred_lst), pred_sets)]





## -- CORE SET SUBSAMPLING OF SNPS -----------------------------------------
if (isTRUE(nchar(core_set) != 0) && data_type == "Inbred") {

  # Get the random core set of genotypes covering the genetic target space
  # of a pre-specified size (as a fraction of genotypes covered by mRNA data).
  core_genotypes <- "./data/derived/predictor_subsets/existing_mrnas.RDS" %>%
    readRDS() %>%
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
    map_at("mrna", .f = function(x) {
      x[rownames(x) %in% core_genotypes, ]
    }) %>%
    map_at("snp", .f = function(x) {
      x[rownames(x) %in% common_genotypes, ]
    }) %>%
    map_at("snp",
           .f = ~sspredr::ensure_snp_quality(
             ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
             any_missing = FALSE, remove_duplicated = TRUE
             )
    )

}




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
if (isTRUE(data_type == "Hybrid")) {
  covered_genotypes <- pheno %>%
    separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
    filter(Dent %in% pred_genotypes & Flint %in% pred_genotypes) %>%
    unite(Hybrid, Dent, Flint, sep = "_") %>%
    mutate(Hybrid = paste0("DF_", Hybrid)) %>%
    select(Hybrid) %>%
    unique() %>%
    flatten_chr()
} else if (isTRUE(data_type == "Inbred")) {
  covered_genotypes <- pred_genotypes
}

pheno <- pheno %>%
  filter(Genotype %in% covered_genotypes)



if (isTRUE(data_type == "Hybrid")) {
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
  pred_lst <- pred_lst %>%
    map_if(.p = names(.) == "ped",
           .f = split_into_hetgroups, y = hybrid_parents, pedigree = TRUE) %>%
    map_if(.p = names(.) != "ped",
           .f = split_into_hetgroups, y = hybrid_parents, pedigree = FALSE) %>%
    map_at("snp", .f = ~map(., ~sspredr::ensure_snp_quality(
      ., callfreq_check = FALSE, maf_check = TRUE,
      maf_threshold = 0.05, any_missing = FALSE, remove_duplicated = TRUE
      )))

  # Add the names of the parental hybrids to the objects.
  # The procedure is necessary for matching predictor data with agronomic data
  # throughout all predictions.
  pre_eta <- pred_lst %>%
    transpose() %>%
    add_parental_hybrid_names()

} else if (isTRUE(data_type == "Inbred")) {

  pre_eta <- list(Inbred = pred_lst)
}

## -- AGRONOMIC DATA PREPARATION ------------------------------------------
# If only one trait was specified, use only that trait, otherwise select all
# traits.
if (nchar(traits) == 0 || length(traits) != 1) {
  traits <- pheno %>% select(Trait) %>% flatten_chr() %>% unique()
} else {
  traits <- traits 
}

# In the case of subsampling based on core samples, we may also have to subset
# the agronomic data afterwards.
# Therefore, we need to determine the largest set of genotypes covered by any
# selected predictor.
if (isTRUE(data_type == "Inbred")) {
  pre_eta[[1]][["geno"]] <- pred_genotypes
} else if (isTRUE(data_type == "Hybrid")) {
  pred_genotypes <- hybrid_names
}

pheno_mat <- pheno %>%
  filter(Trait %in% traits,
         Genotype %in% pred_genotypes) %>%
  spread(key = Trait, value = Value) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "Genotype") %>%
  as.matrix()




## -- GENERATE BGLR() INPUT MATRICES --------------------------------------
# Determine the number of predictors to decide later on, which function to
# call for the set-up of the kernels.
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

# Substitute hybrid names for the current names of the hybrid parents in the
# ETA objects. This is required for the cross-validation function to detect
# correct sorting of genotypes in the variance-covariance matrices with respect
# to the phenotypic traits.
if (isTRUE(data_type == "Hybrid")) {
  eta <- lapply(eta, FUN = function(x) {
    rownames(x$X) <- hybrid_names
    x
  })
}




## -- PREDICTION -----------------------------------------------------------
### Load pre-defined runs
pre_def_runs <- "./data/derived/predictor_subsets/loocv_samples.RDS" %>%
  readRDS()



# Define which runs shall be used, i.e., which genotypes shall be included as 
# test sets.
if (nchar(runs) != 0) {
  run_length <- runs %>%
    strsplit(., split = "-") %>%
    as_vector() %>%
    as.integer() %>%
    invoke(.f = seq, .x = ., by = 1) 
} else {
  run_length <- seq_len(nrow(pheno_mat))
}
param_df <- expand.grid(
  Tst_Geno = pre_def_runs %>% pull(TST_Geno) %>% unique() %>% .[run_length],
  Trait = traits,
  Iter = iter,
  Loocv_Run = pre_def_runs %>% pull(Iter) %>% unique(),
  Runs = length(run_length),
  Core_Number = core_nmb
)
param_df$Trait <- as.character(param_df$Trait)

if (isTRUE(data_type == "Hybrid")) {
  mother_idx <- 2
  father_idx <- 3
  split_char <- "_"
} else if (isTRUE(data_type == "Inbred")) {
  mother_idx <- NULL
  father_idx <- NULL
  split_char <- NULL
}


# Keep track of how long a job is running.
start_time <- Sys.time()
yhat_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {

  # Test genotype
  tst <- param_df %>%
    slice(i) %>%
    pull(Tst_Geno) %>%
    as.character()

  # Miscellaneous parameters
  loocv_run <- param_df %>%
    slice(i) %>%
    pull(Loocv_Run) %>%
    as.character()

  # Training (sub)set
  trn <- pre_def_runs %>%
    filter(
      Iter == loocv_run,
      TST_Geno == tst
    ) %>%
  pull(TRN_Geno)

  # Test and training genotypes
  trn_tst <- c(tst, trn)

  # Operations required for input checks.
  eta <- lapply(eta, FUN = function(x) {
    sub_eta <- x$X
    x$X <- sub_eta[match(trn_tst, rownames(sub_eta)), , drop = FALSE]
    x
  })

  nrow_eta <- unique(vapply(eta, FUN = function(x) {
    nrow(x[["X"]])
  }, FUN.VALUE = integer(1)))

  pheno_mat <- pheno_mat[match(trn_tst, rownames(pheno_mat)), , drop = FALSE]

  # Ensure that the elements of ETA and the phenotyipic values are in the same
  # order.
  eta_rownms <- eta %>%
    map("X") %>%
    map(rownames) %>%
    reduce(intersect) 
  stopifnot(identical(rownames(pheno_mat), eta_rownms))

  # Trait
  trait <- param_df %>%
    slice(i) %>%
    pull(Trait) %>% 
    unique()

  # Set values in the test set as missing.
  y <- pheno_mat[, trait]
  y[tst] <- NA_real_
  stopifnot(sum(is.na(y)) == 1)

  # BGLR implementation
  res <- data.frame(
    Phenotype = NA_character_,
    Geno = NA_character_,
    Mother = NA_character_,
    Father = NA_character_,
    y = NA_complex_,
    yhat = NA_complex_,
    Loocv_run = loocv_run,
    Core_run = core_nmb
  )

  # run the model (GBLUP)
  mod_BGLR <- BGLR::BGLR(
    y = y,
    ETA = eta,
    nIter = iter,
    burnIn = iter / 2,
    saveAt = temp,
    verbose = FALSE
  )

  idx <- which(names(mod_BGLR$yHat) == tst)
  yhat <- mod_BGLR$yHat
  sd_yhat <- mod_BGLR$SD.yHat

  # store results
  res$Phenotype <- trait
  res$Geno <- tst
  res$y <- pheno_mat[rownames(pheno_mat) == tst, trait]
  res$yhat <- yhat[idx]
  res$var_yhat <- sd_yhat[idx] ^ 2
  res
}, mc.cores = use_cores)
res <- rbindlist(yhat_lst)


# If not all three predictors were used, declare the missing ones as "none" and
# add them, which will facilitate further data analyses.
max_pred_number <- 2L
pred_nms <- c(pred_sets,
              rep("none", times = max_pred_number - length(pred_sets)))
res <- res %>%
  dplyr::rename(
    Trait = Phenotype,
    Dent = Mother,
    Flint = Father) %>%
  mutate(
    Data_Type = data_type,
    Pred1 = pred_nms[1],
    Pred2 = pred_nms[2],
    Core_Set = core_set,
    Job_ID = job_id,
    Runs = runs,
    Elapsed_Time = get_elapsed_time(start_time),
    Start_Time = as.character(start_time),
    Cores = use_cores
  ) %>%
  as.data.table()


log_file <- unique(res[, .(Job_ID, Pred1, Pred2, Data_Type, 
                           Core_Set, Elapsed_Time, Iter, CV, Start_Time, 
                           Runs, Cores), ])

# Reduce the size of the prediction object to the minimum possible size.
res[, `:=` (Iter = NULL, 
            Pred1 = NULL,
            Pred2 = NULL,
            Data_Type = NULL,
            Core_Set = NULL,
            CV = NULL,
            Dent = NULL,
            Flint = NULL,
            Runs = NULL,
            Elapsed_Time = NULL,
            Start_Time = NULL,
            Cores = NULL),
  ]


# Prediction results file
saveRDS(res, 
        file = paste0("./data/derived/predictions/", job_id, ".RDS"),
        compress = FALSE)
                     
# Log file
log_location <- "./data/derived/uhoh_maizego_prediction_log.txt"
write.table(log_file,
            file = log_location,
            sep = "\t", row.names = FALSE,
            col.names = ifelse(file.exists(log_location), yes = FALSE,
                               no = TRUE),
            append = ifelse(file.exists(log_location), yes = TRUE, no = FALSE))

