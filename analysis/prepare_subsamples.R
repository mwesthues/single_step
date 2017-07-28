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
saveRDS(pred_lst, "./data/derived/predictor_subsets/pred_lst.RDS")

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
saveRDS(geno_df, "./data/derived/predictor_subsets/geno_df.RDS")

# data frame with all combinations that need to be tested
geno_param_df <- geno_df %>%
  select(-G) %>%
  unique()



sample_by_combination <- function(i) {

  param_subset <- slice(geno_param_df, i)

  material_i <- geno_param_df %>% 
    slice(i) %>%
    pull(Material)

  param_genos <- geno_df %>%
    inner_join(
      y = param_subset,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    pull(G) %>%
    gsub(pattern = "DF_", x = ., replacement = "")
  
  loocv_samples <- sample_loo_sets(
  geno = param_genos,
  frac = 0.7,
  material = material_i,
  iter = 50
  ) %>%
  mutate(TST_Geno = paste0("DF_", TST_Geno))
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

# Save the result in two separate data frames to reduce memory requirements when
# loading only a subset.
rnd_level2_df %>%
  filter(Rnd_Level2 %in% seq_len(25)) %>%
  saveRDS(., "./data/derived/predictor_subsets/rnd_level2_1-25.RDS")
rnd_level2_df %>%
  filter(Rnd_Level2 %in% seq(from = 26, to = 50, by = 1)) %>%
  saveRDS(.,  "./data/derived/predictor_subsets/rnd_level2_26-50.RDS")








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

spec_eta_df <- list(
  eta_spec1,
  eta_spec2,
  eta_spec3,
  eta_spec4,
  eta_spec5,
  eta_spec6
  ) %>%
  bind_rows() %>%
  as_data_frame()
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
  scen_i <- slice(spec_eta_df, i)
  pred1 <- pull(scen_i, Pred1)
  pred2 <- pull(scen_i, Pred2)
  extent <- pull(scen_i, Extent)
  pred_model <- "BRR"
  if (pull(scen_i, Material) == "Hybrid") {
    material_type <- "hyb"
  } else {
    material_type <- "inb"
  }
  
  pred_combi <- list(pred1, pred2) %>%
     keep(nchar(.) != 0) %>%
     flatten_chr() %>%
     paste(., collapse = "_")
   
   pred_sets <- pred_combi %>%
     strsplit(split = "_") %>%
     flatten_chr()

  pred_lst_selector <- paste(material_type, pred_sets, sep = "_")

  # Get the random core set of genotypes covering the genetic target space
  # of a pre-specified size (as a fraction of genotypes covered by mRNA data).
  tmp_geno <- inner_join(
    rnd_level1_df,
    scen_i,
    by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
  ) %>%
  pull(G)

  if (isTRUE(nchar(pred2) == 0 && pred1 != "mrna")) {
    full_geno <- tmp_geno
    core_geno <- tmp_geno
  } else if (pred1 == "mrna" && pred2 == "" && extent == "Full") {
    full_geno <- tmp_geno
    core_geno <- scen_i %>%
      mutate(Extent = "Core") %>%
      inner_join(
        rnd_level1_df,
        by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
      ) %>%
      pull(G)
  } else if (all(nchar(pred_sets) != 0) && extent == "Core") {
    core_geno <- tmp_geno
    full_geno <- inner_join(
      geno_df,
      scen_i,
      by = c("Extent", "Material", "Scenario")
    ) %>%
    pull(G)
  } else if (all(nchar(pred_sets) != 0) && extent == "Full") {
    full_geno <- tmp_geno
    core_geno <- scen_i %>%
      mutate(Extent = "Core") %>%
      inner_join(
        rnd_level1_df,
        by = c("Extent", "Material", "Scenario", "Core_Fraction", "Rnd_Level1")
      ) %>%
      pull(G)
  } 


  # Reduce the predictor data to match the size of the pre-specified core sets.
  # Ensure that SNP quality checks are applied to the genomic data in order to 
  # have only polymorphic markers and no markers in perfect LD.
  if (isTRUE(pull(scen_i, Material) == "Inbred")) {
    curr_pred_lst <- pred_lst %>%
      keep(names(.) %in% pred_lst_selector) %>%
      set_names(
        nm = gsub(pattern = "^[^_]*_", replacement = "", x = names(.))
      ) %>%
      modify_at("mrna", .f = function(x) {
        x[rownames(x) %in% core_geno, ]
      }) %>%
      modify_at("snp", .f = function(x) {
        x[rownames(x) %in% full_geno, ]
      }) %>%
      modify_at("snp",
        .f = ~sspredr::ensure_snp_quality(
          ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
          any_missing = FALSE, remove_duplicated = TRUE
        )
      )  
  }


  if (isTRUE(pull(scen_i, Material) == "Hybrid")) {
    hybrid_parents <- full_geno %>%
      as_data_frame() %>%
      rename(Genotype = "value") %>%
      separate(Genotype, into = c("Prefix", "Dent", "Flint"), sep = "_") %>%
      select(Dent, Flint) %>%
      unique() %>%
      gather(key = Pool, value = Genotype) %>%
      split(.$Pool) %>%
      map(., ~select(., Genotype)) %>%
      map(flatten_chr)
  
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
      keep(names(.) %in% pred_lst_selector) %>%
      set_names(
        nm = gsub(pattern = "^[^_]*_", replacement = "", x = names(.))
      ) %>%
      modify_if(.p = names(.) == "ped",
             .f = split_into_hetgroups, y = hybrid_parents, pedigree = TRUE) %>%
      modify_if(.p = names(.) != "ped",
             .f = split_into_hetgroups, y = hybrid_parents, pedigree = FALSE) %>%
      modify_at("snp", .f = ~map(., ~sspredr::ensure_snp_quality(
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
  
  } else if (isTRUE(pull(scen_i, Material) == "Inbred")) {

    pre_eta <- list(Inbred = curr_pred_lst)
    pre_eta[[1]][["geno"]] <- full_geno
  }
  
  
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

  ## Save individual files.
  data_location <- paste0("./data/derived/predictor_subsets/eta_", i, ".RDS")
  saveRDS(
    object = eta,
    file = data_location
  )
})

# Safe the log file.
log_df <- scenario_seq %>%
  map(function(i) {
    slice(spec_eta_df, i)
  }) %>%
  bind_rows() %>%
  mutate(Data_Location = paste0(
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

