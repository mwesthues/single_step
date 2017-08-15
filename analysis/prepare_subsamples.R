## Goal:
# Sample core sets of genotypes for the Yan-lab data in increments of ten
# percentage points.

## Topics
# Level 2 sampling
# Level 1 sampling
# ETA setup


## -- PACKAGES AND FUNCTIONS -----------------------------------------------
options(max.print = 2000)
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("mwesthues/sspredr", update = FALSE)
#devtools::install_github("tidyverse/dplyr")
pacman::p_load(
  "tidyverse", "methods", "parallel", "dtplyr", "data.table", "uuid"
)
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



## -- GENOTYPES --------------------------------------------------------------
# Define which genotypes belong to which material, extent and scenario.
## hybrids
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

## inbred lines
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




## -- LEVEL 2 SAMPLING -------------------------------------------------------
# For each combination of 'Extent', 'Material' and 'Scenario' sample 70% of
# the available genotypes at random.
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
  iter = 20
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
  saveRDS(., "./data/derived/predictor_subsets/rnd_level2_1-20.RDS")






## -- LEVEL 1 SAMPLING ----------------------------------------------------
# level 1 sampling is only required for the core set of inbred lines from the
# maize diversity panel
level2_geno_df <- rnd_level2_df %>%
  dplyr::filter(Extent == "Core", Material == "Inbred") %>%
  dplyr::select(-Extent, -Material)

level1_frame <- geno_df %>%
  dplyr::filter(Extent == "Core", Material == "Inbred") %>%
  dplyr::distinct(Scenario)

# Sample a fraction of training set individuals based on the core set size,
# which is specified as a fraction of 1.
sample_core <- function(x) {
  smp_fraction <- x %>%
    dplyr::pull(Core_Fraction) %>%
    base::unique()
  sample_frac(x, size = smp_fraction)
}


core_location <- "./data/derived/predictor_subsets/core_set"
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
    base::split(
      list(.$Scenario, .$Rnd_Level1, .$Rnd_Level2, .$TST_Geno),
      drop = TRUE
    ) %>%
    purrr::map(sample_core) %>%
    dplyr::bind_rows() %>%
    saveRDS(
      object = .,
      file = paste0(core_location, i, ".RDS")
    )
})




## -- LEVEL 1 AND LEVEL 2 CONCATENATION --------------------------------------
# Specify all possible scenarios that need to be considered. This is crucial
# for augmenting previously generated data frames with information on which
# genotypes are part of the training or test set (or none of the two) for each
# combination of 'Material', 'Extent', 'Scenario', 'Rnd_Level1' and
# 'Rnd_Level2', respectively.
# This is crucial because ETA objects will differ based on which predictor in
# 'Pred1' and 'Pred2', respectively, are involved.
scen_spec1 <- expand.grid(
  Extent = "Core",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "mrna",
  Scenario = c("A", "B"),
  Core_Fraction = seq(from = 0.1, to = 0.9, by = 0.1),
  stringsAsFactors = FALSE
)

scen_spec2 <- expand.grid(
  Extent = "Core",
  Material = "Inbred",
  Pred1 = c("snp", "mrna"),
  Pred2 = "",
  Scenario = c("A", "B"),
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec3 <- expand.grid(
  Extent = "Full",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "",
  Scenario = c("A", "B"),
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec4 <- expand.grid(
  Extent = "Full",
  Material = "Inbred",
  Pred1 = "snp",
  Pred2 = "mrna",
  Scenario = c("A", "B"),
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec5 <- expand.grid(
  Extent = "Core",
  Material = "Hybrid",
  Pred1 = c("snp", "ped", "mrna"),
  Pred2 = "",
  Scenario = "None",
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec6 <- expand.grid(
  Extent = "Full",
  Material = "Hybrid",
  Pred1 = c("snp", "ped"),
  Pred2 = "mrna",
  Scenario = "None",
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec7 <- expand.grid(
  Extent = "Full",
  Material = "Hybrid",
  Pred1 = c("snp", "ped"),
  Pred2 = "",
  Scenario = "None",
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)

scen_spec8 <- expand.grid(
  Extent = "Full",
  Material = "Hybrid",
  Pred1 = "ped",
  Pred2 = "snp",
  Scenario = "None",
  Core_Fraction = 1,
  stringsAsFactors = FALSE
)
scen_df <- list(
  scen_spec1,
  scen_spec2,
  scen_spec3,
  scen_spec4,
  scen_spec5,
  scen_spec6,
  scen_spec7,
  scen_spec8
  ) %>%
  dplyr::bind_rows() %>%
  data.table::as.data.table() %>%
  dplyr::mutate(Core_Fraction = as.numeric(Core_Fraction))
saveRDS(scen_df, file = "./data/derived/predictor_subsets/scen_df.RDS")

# Collect the data from the level 2 randomization, the level 1 randomization
# and combine them in a single data frame.
filter_files <- function(file_loc) {
  file_loc %>%
    readRDS() %>%
    dplyr::mutate(Rnd_Level2 = as.integer(Rnd_Level2)) %>%
    dplyr::filter(
      Rnd_Level1 %in% seq_len(20),
      Rnd_Level2 %in% seq_len(20)
    )
}

rnd_level1_df <- core_fractions %>%
  purrr::map(.f = ~paste0(core_location, ., ".RDS")) %>%
  purrr::map(filter_files) %>%
  bind_rows()


# Combine the correct randomizations for all materials, extents and scenarios.
rnd_df <- rnd_level2_df %>%
  dplyr::select(-ind) %>%
  dplyr::mutate(
    Core_Fraction = 1,
    Rnd_Level1 = 0,
    Rnd_Level2 = as.integer(Rnd_Level2)
  ) %>%
  dplyr::bind_rows(rnd_level1_df %>% dplyr::select(-ind))
saveRDS(rnd_df, "./data/derived/predictor_subsets/rnd_df.RDS")






## -- THREE MAJOR SCENARIOS ---------------------------------------------------
# 1) Inbred, Core_Fraction == 1
inb_level2 <- rnd_df %>%
  dplyr::filter(Core_Fraction == 1, Material == "Inbred") %>%
  dplyr::full_join(
    y = dplyr::filter(scen_df, Material == "Inbred", Core_Fraction == 1),
    by = c("Extent", "Material", "Scenario", "Core_Fraction")
  )
inb_level2 %>%
  saveRDS(file = "./data/derived/predictor_subsets/inbred_level2_table.RDS")
rm(inb_level2)

# 2) Inbred, Core_Fraction != 1
scen_inb_level1 <- scen_df %>%
  dplyr::filter(Extent == "Core", Material == "Inbred", Core_Fraction != 1)
inb_level1 <- rnd_df %>%
  dplyr::filter(Extent == "Core", Material == "Inbred", Core_Fraction != 1) %>%
  dplyr::full_join(
    y = scen_inb_level1,
    by = c("Extent", "Material", "Scenario", "Core_Fraction")
  ) %>%
  dplyr::filter(complete.cases(.))
inb_level1 %>%
  saveRDS(file = "./data/derived/predictor_subsets/inbred_level1_table.RDS")
rm(inb_level1)

# 3) Hybrid
hyb <-  rnd_df %>%
  dplyr::filter(Material == "Hybrid") %>%
  dplyr::full_join(
    y = dplyr::filter(scen_df, Material == "Hybrid"),
    by = c("Extent", "Material", "Scenario", "Core_Fraction")
  )
  dplyr::filter(complete.cases(.))
hyb %>%
  saveRDS(file = "./data/derived/predictor_subsets/hybrid_table.RDS")
rm(hyb)






## -- HYBRID ETA SETUP -----------------------------------------------------
pred_lst <- readRDS("./data/derived/predictor_subsets/pred_lst.RDS")
geno_df <- readRDS( "./data/derived/predictor_subsets/geno_df.RDS")
scen_df <- readRDS("./data/derived/predictor_subsets/scen_df.RDS")

hyb_pre_eta <- geno_df %>%
  dplyr::filter(Material == "Hybrid") %>%
  dplyr::full_join(
    y = scen_df %>%
          dplyr::filter(Material == "Hybrid") %>%
          dplyr::select(-Core_Fraction) %>%
          tibble::as_data_frame(),
    by = c("Extent", "Material", "Scenario")
  ) %>%
  tidyr::unite(col = "Predictor", c("Pred1", "Pred2"), remove = TRUE) %>%
  dplyr::mutate(Predictor = gsub("_$", x = Predictor, replacement = ""))



# Generate a sequence of all possible scenarios for easy looping.
# Additionally, generate a unique ID for each combination of parameters for
# which a separate ETA object will be created.
hyb_unique_eta_combis <- hyb_pre_eta %>%
  dplyr::distinct(Material, Extent, Predictor, .keep_all = FALSE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(UUID = uuid::UUIDgenerate())
saveRDS(
  hyb_unique_eta_combis,
  file = "./data/derived/predictor_subsets/hyb_unique_eta_combis.RDS"
)
hyb_scenario_seq <- hyb_unique_eta_combis %>%
  nrow() %>%
  seq_len()

# Specify which BGLR algorithm shall be used.
pred_model <- "BRR"


lapply(hyb_scenario_seq, FUN = function(i) {
  # For setting up ETA objects it is crucial to know which predictors are in
  # use.
  current_scenario <- dplyr::slice(hyb_unique_eta_combis, i)
  pred_combi <- dplyr::pull(current_scenario, Predictor)
  pred_sets <- pred_combi %>%
    base::strsplit(split = "_") %>%
    purrr::flatten_chr()
  pred_lst_selector <-  current_scenario %>%
    dplyr::mutate(Material = dplyr::if_else(
      Material == "Hybrid",
      true = "hyb",
      false = "inb"
    )) %>%
    dplyr::pull(Material) %>%
    paste(., pred_sets, sep = "_")


  # Get the names of the hybrids and their parental components to properly
  # augment the ETA objects to match their phenotypic data.
  current_hyb_df <- hyb_pre_eta %>%
    dplyr::inner_join(
      y = current_scenario,
      by = c("Extent", "Predictor", "Material")
  )
  hybrid_names <- dplyr::pull(current_hyb_df, G)
  hybrid_parents <- current_hyb_df %>%
    tidyr::separate(
      G,
      into = c("Prefix", "Dent", "Flint"),
      sep = "_"
    ) %>%
    dplyr::select(Dent, Flint) %>%
    unique() %>%
    tidyr::gather(key = Pool, value = Genotype) %>%
    base::split(.$Pool) %>%
    purrr::map(., ~select(., Genotype)) %>%
    purrr::map(.f = ~purrr::flatten_chr(.))

  if (isTRUE(length(pred_sets) == 2 && pred_sets[1] == "ped")) {
    snp_hybrid_parents <- hyb_pre_eta %>%
      dplyr::inner_join(
        y = hyb_unique_eta_combis %>%
              dplyr::filter(
                Material == "Hybrid",
                Extent == "Core",
                Predictor == "ped"
              ),
        by = c("Extent", "Predictor", "Material")
      ) %>%
      tidyr::separate(
        G,
        into = c("Prefix", "Dent", "Flint"),
        sep = "_"
      ) %>%
      dplyr::select(Dent, Flint) %>%
      tidyr::gather(key = Pool, value = Genotype) %>%
      base::split(.$Pool) %>%
      purrr::map(., ~select(., Genotype)) %>%
      purrr::map(.f = ~purrr::flatten_chr(.))
  } else if (isTRUE(!"ped" %in% pred_sets)) {
    snp_hybrid_parents <- hybrid_parents
  }

  #if (exists("hybrid_parents")) {
  #  print(paste("'hybrid_parents' exists for ", i))
  #} else if (!exists("hybrid_parents")) {
  #  print(paste("'hybrid_parents' does not exist! for", i))
  #}
  # In the case of hybrid data, we still need to split all predictor matrices
  # into Flint and Dent components first.
  # 1. modify_if
  # In the case of pedigree data, we need to split the data by genotypes in the
  # x- as well as the y-dimension because there are no features, which is
  # different for other predictor matrices.
  # 2. modify_at("snp")
  # Ensure that SNP quality checks are applied separately to the genomic data
  # for each heterotic group in order to have only polymorphic markers and no
  # markers in perfect LD.
  curr_pred_lst <- pred_lst %>%
    purrr::keep(names(.) %in% pred_lst_selector) %>%
    purrr::set_names(
      nm = gsub(pattern = "^[^_]*_", replacement = "", x = names(.))
    ) %>%
    purrr::modify_if(.p = names(.) == "ped",
      .f = split_into_hetgroups, y = hybrid_parents, pedigree = TRUE
    ) %>%
    purrr::modify_if(.p = !names(.) %in% c("ped", "snp"),
      .f = split_into_hetgroups, y = hybrid_parents, pedigree = FALSE
    ) %>%
    purrr::modify_if(.p = names(.) == "snp",
      .f = split_into_hetgroups, y = snp_hybrid_parents, pedigree = FALSE
    ) %>%
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
    purrr::modify_if(.p = names(.) == "Dent", .f = function(x) {
      append(x, list(geno = hybrid_parents$Dent))
    }) %>%
    purrr::modify_if(.p = names(.) == "Flint", .f = function(x) {
      append(x, list(geno = hybrid_parents$Flint))
    })

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

  current_uuid <- dplyr::pull(current_scenario, UUID)
  saveRDS(
    eta,
    file = paste0(
      "./data/derived/predictor_subsets/ETA_",
      current_uuid,
      ".RDS"
    )
  )
  print(paste("Successfully processed element", i))
})





## -- INBRED ETA SETUP ------------------------------------------------------
pred_lst <- readRDS("./data/derived/predictor_subsets/pred_lst.RDS")
geno_df <- readRDS( "./data/derived/predictor_subsets/geno_df.RDS")
scen_df <- readRDS("./data/derived/predictor_subsets/scen_df.RDS")

inb_level2_pre_eta <- geno_df %>%
  dplyr::filter(Material == "Inbred") %>%
  dplyr::full_join(
    y = scen_df %>%
          dplyr::filter(Material == "Inbred", Core_Fraction == 1) %>%
          dplyr::select(-Core_Fraction) %>%
          tibble::as_data_frame(),
    by = c("Extent", "Material", "Scenario")
  ) %>%
  tidyr::unite(col = "Predictor", c("Pred1", "Pred2"), remove = TRUE) %>%
  dplyr::mutate(Predictor = gsub("_$", x = Predictor, replacement = ""))



# Generate a sequence of all possible scenarios for easy looping.
# Additionally, generate a unique ID for each combination of parameters for
# which a separate ETA object will be created.
inb_level2_unique_eta_combis <- inb_level2_pre_eta %>%
  dplyr::distinct(Material, Extent, Predictor, Scenario, .keep_all = FALSE) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(UUID = uuid::UUIDgenerate())
saveRDS(
  inb_level2_unique_eta_combis,
  file = "./data/derived/predictor_subsets/inb_level2_unique_eta_combis.RDS"
)
inb_level2_scenario_seq <- inb_level2_unique_eta_combis %>%
  nrow() %>%
  seq_len()

# Specify which BGLR algorithm shall be used.
pred_model <- "BRR"


lapply(inb_level2_scenario_seq, FUN = function(i) {
  # For setting up ETA objects it is crucial to know which predictors are in
  # use.
  current_scenario <- dplyr::slice(inb_level2_unique_eta_combis, i)
  pred_combi <- dplyr::pull(current_scenario, Predictor)
  pred_sets <- pred_combi %>%
    base::strsplit(split = "_") %>%
    purrr::flatten_chr()
  pred_lst_selector <-  current_scenario %>%
    dplyr::mutate(Material = dplyr::if_else(
      Material == "Hybrid",
      true = "hyb",
      false = "inb"
    )) %>%
    dplyr::pull(Material) %>%
    paste(., pred_sets, sep = "_")

  # Get the names of the hybrids and their parental components to properly
  # augment the ETA objects to match their phenotypic data.
  current_inb_df <- inb_level2_pre_eta %>%
    dplyr::inner_join(
      y = current_scenario,
      by = c("Extent", "Predictor", "Material", "Scenario")
  )
  full_geno <- dplyr::pull(current_inb_df, G)

  if (length(pred_sets) != 1 || "mrna" %in% pred_sets) {
    core_geno <- current_scenario %>%
      dplyr::mutate(
        Extent = "Core",
        Predictor = "mrna"
      ) %>%
      dplyr::inner_join(
        y = inb_level2_pre_eta,
        by = c("Extent", "Predictor", "Material", "Scenario")
    ) %>%
    dplyr::pull(G)
  }

  # In the case of hybrid data, we still need to split all predictor matrices
  # into Flint and Dent components first.
  # 1. modify_if
  # In the case of pedigree data, we need to split the data by genotypes in the
  # x- as well as the y-dimension because there are no features, which is
  # different for other predictor matrices.
  # 2. modify_at("snp")
  # Ensure that SNP quality checks are applied separately to the genomic data
  # for each heterotic group in order to have only polymorphic markers and no
  # markers in perfect LD.
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
    purrr::modify_at("snp", .f = ~sspredr::ensure_snp_quality(
      ., callfreq_check = FALSE, maf_check = TRUE,
      maf_threshold = 0.05, any_missing = FALSE, remove_duplicated = TRUE
      )
    )

  # Make sure that the predictor matrices are in the intended order, which is
  # crucial for the imputation of the predictor matrix that covers fewer
  # genotypes.
  curr_pred_lst <- curr_pred_lst[match(names(curr_pred_lst), pred_sets)]

  # The additionn of names is necessary for matching predictor data with
  # agronomic data throughout all predictions.
  pre_eta <- list(Inbred = curr_pred_lst)
  pre_eta[[1]][["geno"]] <- full_geno

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
  pre_eta <- pre_eta %>%
   purrr::map(.f = function(x) {
     x$is_pedigree <- FALSE
     x$bglr_model <- pred_model
     x$as_kernel <- TRUE
     x
  })


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

  current_uuid <- dplyr::pull(current_scenario, UUID)
  saveRDS(
    eta,
    file = paste0(
      "./data/derived/predictor_subsets/ETA_",
      current_uuid,
      ".RDS"
    )
  )
  print(i)
})









## -- INBRED CORE_SET != 1 ETA SETUP ------------------------------------------
pred_lst <- readRDS("./data/derived/predictor_subsets/pred_lst.RDS")
scen_df <- readRDS("./data/derived/predictor_subsets/scen_df.RDS")
geno_df <-  "./data/derived/predictor_subsets/geno_df.RDS" %>%
  readRDS() %>%
  data.table::as.data.table() %>%
  data.table::setkey()

# Define how the data should be sorted. This will considerably accelerate the
# inner join of the 'current' scenario with the entire data set in order to
# select the corresponding genotype names.
keycols <- c(
  "Extent",
  "Predictor",
  "Scenario",
  "Core_Fraction",
  "TST_Geno",
  "Rnd_Level1",
  "Rnd_Level2"
)

inb_level1_pre_eta <- readRDS(
  "./data/derived/predictor_subsets/inbred_level1_table.RDS"
  ) %>%
  tidyr::unite(Predictor, c("Pred1", "Pred2"), sep = "_") %>%
  data.table::as.data.table() %>%
  data.table::setkeyv(cols = keycols)

inb_level1_unique_eta_combis <- inb_level1_pre_eta %>%
  dplyr::distinct(
    Material,
    Extent,
    Scenario,
    Rnd_Level1,
    Rnd_Level2,
    Core_Fraction,
    TST_Geno,
    Predictor
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(UUID = uuid::UUIDgenerate()) %>%
  dplyr::ungroup()

inb_level1_scenario_seq <- inb_level1_unique_eta_combis %>%
  nrow() %>%
  seq_len()

# Specify which BGLR algorithm shall be used.
pred_model <- "BRR"



system.time(
lapply(inb_level1_scenario_seq, FUN = function(i) {
  # For setting up ETA objects it is crucial to know which predictors are in
  # use.
  current_scenario <- dplyr::slice(inb_level1_unique_eta_combis, i) %>%
    data.table::as.data.table() %>%
    data.table::setkeyv(cols = keycols)
  pred_combi <- dplyr::pull(current_scenario, Predictor)
  pred_sets <- pred_combi %>%
    base::strsplit(split = "_") %>%
    purrr::flatten_chr()
  pred_lst_selector <- paste("inb", pred_sets, sep = "_")

  # Get the names of the hybrids and their parental components to properly
  # augment the ETA objects to match their phenotypic data.
  # This implementation usies data.table because its use of keys and the size
  # of the two data sets to be merged make it a lot faster than the dplyr
  # solution.
  core_geno <- inb_level1_pre_eta[current_scenario, nomatch = 0] %>%
    .[, .(TST_Geno, TRN_Geno), ] %>%
    tidyr::gather(key = Set, value = Genotype) %>%
    dplyr::distinct(Genotype) %>%
    dplyr::pull(Genotype)

  full_geno <- current_scenario %>%
    dplyr::select(Extent, Material, Scenario) %>%
    dplyr::inner_join(
      y = geno_df,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    dplyr::pull(G)

  # In the case of hybrid data, we still need to split all predictor matrices
  # into Flint and Dent components first.
  # 1. modify_if
  # In the case of pedigree data, we need to split the data by genotypes in the
  # x- as well as the y-dimension because there are no features, which is
  # different for other predictor matrices.
  # 2. modify_at("snp")
  # Ensure that SNP quality checks are applied separately to the genomic data
  # for each heterotic group in order to have only polymorphic markers and no
  # markers in perfect LD.
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
    purrr::modify_at("snp", .f = ~sspredr::ensure_snp_quality(
      ., callfreq_check = FALSE, maf_check = TRUE,
      maf_threshold = 0.05, any_missing = FALSE, remove_duplicated = TRUE
      )
    )

  # Make sure that the predictor matrices are in the intended order, which is
  # crucial for the imputation of the predictor matrix that covers fewer
  # genotypes.
  curr_pred_lst <- curr_pred_lst[match(names(curr_pred_lst), pred_sets)]

  # The additionn of names is necessary for matching predictor data with
  # agronomic data throughout all predictions.
  pre_eta <- list(Inbred = curr_pred_lst)
  pre_eta[[1]][["geno"]] <- full_geno

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
  pre_eta <- pre_eta %>%
   purrr::map(.f = function(x) {
     x$is_pedigree <- FALSE
     x$bglr_model <- pred_model
     x$as_kernel <- TRUE
     x
  })


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

  current_uuid <- dplyr::pull(current_scenario, UUID)
  saveRDS(
    eta,
    file = paste0(
      "./data/derived/predictor_subsets/ETA_",
      current_uuid,
      ".RDS"
    ),
    compress = FALSE
  )
})


