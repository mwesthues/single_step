## Goal:
# Sample core sets of genotypes for the Yan-lab data in increments of ten
# percentage points.

## Topics
# Level 2 sampling
# Level 1 sampling
# ETA setup


## -- PACKAGES AND FUNCTIONS -----------------------------------------------
if (isTRUE(interactive())) {
  .libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.4/"))
}
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("mwesthues/sspredr", update = FALSE)
pacman::p_load("tidyverse", "methods", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# Prediction helper functions (outsourced from this script for improved
# readability).
source("./software/prediction_helper_functions.R")



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
  tibble::as_data_frame() %>%
  dplyr::mutate(Path = as.character(Path))

pred_lst <- named_df %>%
  tibble::as_data_frame() %>%
  dplyr::select(Path) %>%
  purrr::flatten_chr() %>%
  purrr::map(readRDS)

pred_lst_names <- named_df %>%
  dplyr::select(Predictor) %>%
  purrr::flatten_chr()
names(pred_lst) <- pred_lst_names
saveRDS(pred_lst, "./data/derived/predictor_subsets/pred_lst.RDS")

# load genotypes from the first and the fourth cluster of the pca (Scenario 'a')
cluster14 <- readRDS("./data/derived/maizego/cluster_14_genotypes.RDS")



### hybrids
hyb_pheno <- pred_lst %>%
  purrr::keep(names(.) == "hyb_agro") %>%
  .[[1]]

full_hyb_genos <- hyb_pheno %>%
  dplyr::pull(Genotype) %>%
  base::unique()

core_hyb_parents <- pred_lst %>%
  purrr::keep(names(.) == "hyb_mrna") %>%
  .[[1]] %>%
  rownames()

core_hyb_genos <- hyb_pheno %>%
  tidyr::separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
  dplyr::filter(Dent %in% core_hyb_parents & Flint %in% core_hyb_parents) %>%
  tidyr::unite(Hybrid, Dent, Flint, sep = "_") %>%
  dplyr::mutate(Hybrid = paste0("DF_", Hybrid)) %>%
  dplyr::select(Hybrid) %>%
  base::unique() %>%
  purrr::flatten_chr()

### inbred lines
full_inb_bc <- pred_lst %>%
  purrr::keep(names(.) == "inb_snp") %>%
  .[[1]] %>%
  rownames()

full_inb_a <- intersect(full_inb_bc, cluster14)

core_inb_bc <- pred_lst %>%
  purrr::keep(names(.) == "inb_mrna") %>%
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
  tidyr::separate(
    ind,
    into = c("Extent", "Material", "Scenario"),
    sep = "_"
  ) %>%
  dplyr::rename(G = "values") %>%
  tibble::as_data_frame()
saveRDS(geno_df, "./data/derived/predictor_subsets/geno_df.RDS")

# data frame with all combinations that need to be tested
geno_param_df <- geno_df %>%
  dplyr::select(-G) %>%
  base::unique()



sample_by_combination <- function(i) {

  param_subset <- dplyr::slice(geno_param_df, i)

  material_i <- geno_param_df %>%
    dplyr::slice(i) %>%
    dplyr::pull(Material)

  param_genos <- geno_df %>%
    dplyr::inner_join(
      y = param_subset,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    dplyr::pull(G) %>%
    gsub(pattern = "DF_", x = ., replacement = "")

  loocv_samples <- sample_loo_sets(
  geno = param_genos,
  frac = 0.7,
  material = material_i,
  iter = 50
  )
  if (isTRUE(material_i == "Hybrid")) {
    loocv_samples %>%
    dplyr::mutate(TST_Geno = paste0("DF_", TST_Geno))
  } else if (isTRUE(material_i == "Inbred")) {
    loocv_samples
  }
}



rnd_level2_df <- geno_param_df %>%
  nrow() %>%
  seq_len() %>%
  purrr::map(sample_by_combination) %>%
  purrr::set_names(nm = names(geno_lst)) %>%
  dplyr::bind_rows(.id = "ID") %>%
  tidyr::separate(
    ID,
    into = c("Extent", "Material", "Scenario"),
    sep = "_"
  ) %>%
  dplyr::rename(Rnd_Level2 = "Iter")

# Save the result in two separate data frames to reduce memory requirements when
# loading only a subset.
rnd_level2_df %>%
  dplyr::filter(Rnd_Level2 %in% seq_len(25)) %>%
  saveRDS(., "./data/derived/predictor_subsets/rnd_level2_1-25.RDS")
rnd_level2_df %>%
  dplyr::filter(Rnd_Level2 %in% seq(from = 26, to = 50, by = 1)) %>%
  saveRDS(.,  "./data/derived/predictor_subsets/rnd_level2_26-50.RDS")








## -- LEVEL 1 SAMPLING ----------------------------------------------------
# level 1 sampling is only required for the core set of inbred lines from the
# maize diversity panel
level2_geno_df <- rnd_level2_df %>%
  dplyr::filter(Extent == "Core", Material == "Inbred") %>%
  dplyr::select(-Extent, -Material)

level1_frame <- geno_df %>%
  filter(Extent == "Core", Material == "Inbred") %>%
  pull(Scenario) %>%
  unique()

# Sample a fraction of training set individuals based on the core set size,
# which is specified as a fraction of 1.
sample_core <- function(x) {
  smp_fraction <- x %>%
    dplyr::pull(Core_Fraction) %>%
    base::unique()
  sample_frac(x, size = smp_fraction)
}


core_fractions <- seq(from = 0.1, to = 0.9, by = 0.1)
lapply(core_fractions, FUN = function(i) {
  print(i)
  level1_frame %>%
    expand.grid(
      Extent = "Core",
      Material = "Inbred",
      Scenario = .,
      Core_Fraction = i,
      Rnd_Level1 = seq_len(50),
      stringsAsFactors = FALSE
    ) %>%
    tibble::as_data_frame() %>%
    dplyr::right_join(
      x = .,
      y = level2_geno_df,
      by = "Scenario"
    ) %>%
    base::split(list(.$Scenario, .$Rnd_Level1, .$Rnd_Level2, .$TST_Geno), drop = TRUE) %>%
    purrr::map(sample_core) %>%
    dplyr::bind_rows() %>%
    saveRDS(
      object = .,
      file = paste0("./data/derived/predictor_subsets/core_set", i, ".RDS")
    )
})



## Add all other combinations to the data frame. This will facilitate the setup
## of ETA objects later on.
#rnd_level1_df <- geno_df %>%
#  dplyr::mutate(Core_Fraction = 0, Rnd_Level1 = 0) %>%
#  dplyr::bind_rows(., rnd_level1_df) %>%
#  dplyr::mutate_all(as.character)
#saveRDS(rnd_level1_df, "./data/derived/predictor_subsets/rnd_level1.RDS")




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
  Pred1 = "snp",
  Pred2 = "",
  Scenario = c("A", "B"),
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

eta_spec4 <- expand.grid(
  Extent = "Full",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "mrna",
  Scenario = c("A", "B"),
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

eta_spec7 <- expand.grid(
  Extent = "Full",
  Material = "Hybrid",
  Pred1 = c("snp", "ped"),
  Pred2 = "",
  Scenario = "None",
  Rnd_Level1 = "0",
  Core_Fraction = "0",
  stringsAsFactors = FALSE
)

spec_eta_df <- list(
  eta_spec1,
  eta_spec2,
  eta_spec3,
  eta_spec4,
  eta_spec5,
  eta_spec6,
  eta_spec7
  ) %>%
  dplyr::bind_rows() %>%
  tibble::as_data_frame()
saveRDS(spec_eta_df, "./data/derived/predictor_subsets/spec_eta_df.RDS")



pred_lst <- readRDS("./data/derived/predictor_subsets/pred_lst.RDS")
rnd_level1_df <- readRDS("./data/derived/predictor_subsets/rnd_level1.RDS")
rnd_level2_df <- readRDS("./data/derived/predictor_subsets/rnd_level2.RDS")
geno_df <- readRDS("./data/derived/predictor_subsets/geno_df.RDS")
spec_eta_df <- readRDS("./data/derived/predictor_subsets/spec_eta_df.RDS")

scenario_seq <- spec_eta_df %>%
  nrow() %>%
  seq_len()

lapply(scenario_seq, FUN = function(i) {
  scen_i <- dplyr::slice(spec_eta_df, i)
  pred1 <- dplyr::pull(scen_i, Pred1)
  pred2 <- dplyr::pull(scen_i, Pred2)
  extent <- dplyr::pull(scen_i, Extent)
  pred_model <- "BRR"
  if (dplyr::pull(scen_i, Material) == "Hybrid") {
    material_type <- "hyb"
  } else {
    material_type <- "inb"
  }
  
  pred_combi <- list(pred1, pred2) %>%
     purrr::keep(nchar(.) != 0) %>%
     purrr::flatten_chr() %>%
     paste(., collapse = "_")
   
   pred_sets <- pred_combi %>%
     strsplit(split = "_") %>%
     purrr::flatten_chr()

  pred_lst_selector <- paste(material_type, pred_sets, sep = "_")

  # Get the random core set of genotypes covering the genetic target space
  # of a pre-specified size (as a fraction of genotypes covered by mRNA data).
  tmp_geno <- dplyr::inner_join(
    rnd_level1_df,
    scen_i,
    by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
  ) %>%
  dplyr::pull(G)

  if (isTRUE(nchar(pred2) == 0 && pred1 != "mrna")) {
    full_geno <- tmp_geno
    core_geno <- tmp_geno
  } else if (pred1 == "mrna" && pred2 == "" && extent == "Full") {
    full_geno <- tmp_geno
    core_geno <- scen_i %>%
      dplyr::mutate(Extent = "Core") %>%
      dplyr::inner_join(
        rnd_level1_df,
        by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
      ) %>%
      dplyr::pull(G)
  } else if (all(nchar(pred_sets) != 0) && extent == "Core") {
    core_geno <- tmp_geno
    full_geno <- dplyr::inner_join(
      geno_df,
      scen_i,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    dplyr::pull(G)
  } else if (all(nchar(pred_sets) != 0) && extent == "Full") {
    full_geno <- tmp_geno
    core_geno <- scen_i %>%
      dplyr::mutate(Extent = "Core") %>%
      dplyr::inner_join(
        rnd_level1_df,
        by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
      ) %>%
      dplyr::pull(G)
  } 


  # Reduce the predictor data to match the size of the pre-specified core sets.
  # Ensure that SNP quality checks are applied to the genomic data in order to 
  # have only polymorphic markers and no markers in perfect LD.
  if (isTRUE(pull(scen_i, Material) == "Inbred")) {
    curr_pred_lst <- pred_lst %>%
      purrr::keep(names(.) %in% pred_lst_selector) %>%
      purrr::set_names(
        nm = gsub(pattern = "^[^_]*_", replacement = "", x = names(.))
      ) %>%
      purrr::modify_at("mrna", .f = function(x) {
        x[rownames(x) %in% core_geno, ]
      }) %>%
      purrr::modify_at("snp", .f = function(x) {
        x[rownames(x) %in% full_geno, ]
      }) %>%
      purrr::modify_at("snp",
        .f = ~sspredr::ensure_snp_quality(
          ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
          any_missing = FALSE, remove_duplicated = TRUE
        )
      )  
  }


  if (isTRUE(dplyr::pull(scen_i, Material) == "Hybrid")) {
    hybrid_parents <- full_geno %>%
      tibble::as_data_frame() %>%
      dplyr::rename(Genotype = "value") %>%
      tidyr::separate(
        Genotype,
        into = c("Prefix", "Dent", "Flint"),
        sep = "_"
      ) %>%
      dplyr::select(Dent, Flint) %>%
      unique() %>%
      tidyr::gather(key = Pool, value = Genotype) %>%
      base::split(.$Pool) %>%
      purrr::map(., ~select(., Genotype)) %>%
      purrr::map(flatten_chr)
  
    hybrid_names <- full_geno

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
    curr_pred_lst <- pred_lst %>%
      purrr::keep(names(.) %in% pred_lst_selector) %>%
      purrr::set_names(
        nm = gsub(pattern = "^[^_]*_", replacement = "", x = names(.))
      ) %>%
      purrr::modify_if(.p = names(.) == "ped",
             .f = split_into_hetgroups, y = hybrid_parents, pedigree = TRUE) %>%
      purrr::modify_if(.p = names(.) != "ped",
             .f = split_into_hetgroups, y = hybrid_parents, pedigree = FALSE) %>%
      purrr::modify_at("snp", .f = ~map(., ~sspredr::ensure_snp_quality(
        ., callfreq_check = FALSE, maf_check = TRUE,
        maf_threshold = 0.05, any_missing = FALSE, remove_duplicated = TRUE
        )
      )
    )
  
   # Make sure that the predictor matrices are in the intended order, which is
   # crucial for the imputation of the predictor matrix that covers fewer
   # genotypes.
   curr_pred_lst <- curr_pred_lst[match(names(curr_pred_lst), pred_sets)]


    # Add the names of the parental hybrids to the objects.
    # The procedure is necessary for matching predictor data with agronomic data
    # throughout all predictions.
    pre_eta <- curr_pred_lst %>%
      purrr::transpose() %>%
      add_parental_hybrid_names()
  
  } else if (isTRUE(dplyr::pull(scen_i, Material) == "Inbred")) {

    pre_eta <- list(Inbred = curr_pred_lst)
    pre_eta[[1]][["geno"]] <- full_geno
  }
  
  
  # Determine the number of predictors to decide later on, which function to
  # call for the set-up of the kernels.
  raw_args <- pre_eta %>%
    purrr::map(names) %>%
    purrr::reduce(intersect) %>%
    purrr::discard(. == "geno")
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
      purrr::map(.f = function(x) {
        x$as_kernel <- FALSE
        x$bglr_model <- pred_model
        x
    })
  } else if (any(grepl("ped", x = raw_args))) {
    pre_eta <- pre_eta %>%
      purrr::map(.f = function(x) {
        x$is_pedigree <- TRUE
        x$bglr_model <- pred_model
        x$as_kernel <- TRUE
        x
      })
  } else {
    pre_eta <- pre_eta %>%
      purrr::map(.f = function(x) {
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
    purrr::map(names) %>%
    purrr::reduce(intersect)
  current_eta_names[seq_len(pred_number)] <- eta_pred_nms
  pre_eta <- pre_eta %>%
    purrr::map(function(x) {
      names(x) <- current_eta_names
      x
    })
  # Build the kernel for BGLR.
  eta <- purrr::invoke_map(get(eta_fun), pre_eta) %>%
    purrr::flatten()

  ## Save individual files.
  data_location <- paste0("./data/derived/predictor_subsets/eta_", i, ".RDS")
  saveRDS(
    object = eta,
    file = data_location
  )
})

# Safe the log file.
log_df <- scenario_seq %>%
  purrr::map(function(i) {
    dplyr::slice(spec_eta_df, i)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(Data_Location = paste0(
    "./data/derived/predictor_subsets/eta_",
    scenario_seq,
    ".RDS"
  ))

log_location <- "./data/derived/log_prepare_subsamples.txt"
write.table(
  log_df,
  file = log_location,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  append = FALSE
)

