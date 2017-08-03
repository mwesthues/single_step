if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse")

# Load data
rnd_level2_df <- "./data/derived/predictor_subsets/rnd_level2_1-25.RDS" %>%
  readRDS() %>%
  dplyr::filter(Rnd_Level2 %in% seq_len(20))

# Check, if training set hybrids are only T0 hybrids
check_if_t0 <- function(df, tst, trn) {
  # --------------------------------------------------------------------------
  # Concepts from:
  # https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html
  # http://www.win-vector.com/blog/2015/03/using-closures-as-objects-in-r/
  # --------------------------------------------------------------------------

  trn <- rlang::enquo(trn)
  tst <- rlang::enquo(tst)

  # Define a simple function doing the replacement
  gsub_hyb <- function(var_in) {
    gsub("DF_", x = var_in, replacement = "")
  }

  # Extract the parents
  extract_parents <- function(df, genotype) {
    genotype <- rlang::enquo(genotype)
    df %>%
      dplyr::pull(!!genotype) %>%
      tibble::as_data_frame() %>%
      dplyr::rename(Hybrid = "value") %>%
      tidyr::separate(Hybrid, into = c("Dent", "Flint"), sep = "_") %>%
      unique()
  }

  # Split the data frame into separate ones for the test and training set and
  # split the hybrid names into components for Dent and Flint, respectively.
  tst_trn_lst <- df %>%
    dplyr::select(!!trn, !!tst) %>%
    dplyr::mutate(!!rlang::quo_name(tst) := gsub_hyb((!!tst))) %>%
    tidyr::gather(key = Set, value = G) %>%
    dplyr::mutate(Set = gsub("_Geno", x = Set, replacement = "")) %>%
    split(.$Set) %>%
    purrr::map(~extract_parents(., genotype = G))

  # Determine the intersect between training and test set parents, separately
  # for Dent and Flint material.
  tst_trn_intersect <- tst_trn_lst %>%
    purrr::map(~tidyr::gather(., key = Group, value = G)) %>%
    dplyr::bind_rows(.id = "Set") %>%
    split(., list(.$Group)) %>%
    purrr::map(~split(., list(.$Set))) %>%
    purrr::modify_depth(.depth = 2, .f = ~dplyr::pull(., G)) %>%
    purrr::map(~reduce(., intersect))

  isect_length <- tst_trn_intersect %>%
    purrr::flatten_chr() %>%
    length()

  isect_length
#  if (isTRUE(isect_length == 0)) {
#    message("The object 'df' contains only T0 hybrids")
#  } else {
#    tst_trn_intersect
#    warning("The object 'df' contains not only T0 hybrids", .call = FALSE)
#  }
}

rnd_level2_df %>%
  split(list(.$Rnd_Level2, .$TST_Geno)) %>%
  purrr::map(~check_if_t0(., tst = TST_Geno, trn = TRN_Geno)) %>%
  reduce(unique)

tst_trn_lst <- check_if_t0(smp_df, tst = TST_Geno, trn = TRN_Geno)


tst_df <- tst_trn_lst[["TST"]]
trn_df <- tst_trn_lst[["TRN"]]
dplyr::intersect(tst_df, trn_df %>% bind_rows(., data_frame(Dent = "3205", Flint = "A")))
