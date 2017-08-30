# Goal: Bootstrap prediction results.

## -- PACKAGES ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  "dplyr", "tibble", "tidyr", "purrr", "dtplyr", "data.table", "boot", "broom"
)


## -- INPUT PARAMETERS --------------------------------------------------------
# Set the number of cores to use in parallel.
if (interactive()) {

  use_cores <- 2L

} else {

  use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
}


## -- DATA --------------------------------------------------------------------
dat <- "./data/processed/predictions/predictions.RDS" %>%
  readRDS() %>%
  dplyr::mutate(Combi = substr(Interval, start = 1, stop = 3)) %>%
  data.table::as.data.table()


## -- USER FUNCTIONS ----------------------------------------------------------
# Compute the Pearson correlation coefficient for (y, yhat) using a data.table
# object.
# data.table bootstrap: https://stackoverflow.com/a/36731848/2323832
boot_cor <- function(x, i) {

  x[i, cor(y, yhat), ]
}



## -- ANALYSIS ----------------------------------------------------------------
# Get the names for all unique combinations that will be bootstrapped to
# assign them to the list of resuls from the 'boot()' function.
boot_group_nms <- dat %>%
  dplyr::distinct(
    Combi,
    Trait,
    Core_Fraction,
    Predictor,
    Rnd_Level1,
    Rnd_Level2
  ) %>%
  dplyr::mutate(
    Core_Fraction = dplyr::if_else(
      is.na(Core_Fraction),
      true = 1,
      false = Core_Fraction
    )
  ) %>%
  dplyr::mutate(
    Rnd_Level1 = dplyr::if_else(
      is.na(Rnd_Level1),
      true = 1,
      false = Rnd_Level1
      )
  ) %>%
  tidyr::unite(
    col = Boot_Group,
    Combi,
    Trait,
    Core_Fraction,
    Predictor,
    Rnd_Level1,
    Rnd_Level2,
    sep = "-"
  ) %>%
  dplyr::pull(Boot_Group)

# Use these variables to group the prediction results into distinct clusters.
boot_group_vars <- c(
  "Combi",
  "Trait",
  "Core_Fraction",
  "Predictor",
  "Rnd_Level1",
  "Rnd_Level2"
  )

# Run the bootstrap.
set.seed(34094)
bootstrap_result <- dat[,
  list(list(boot::boot(
    .SD,
    statistic = boot_cor,
    R = 1e4,
    parallel = "multicore",
    ncpus = use_cores
    ))),
  by = boot_group_vars
 ]$V1


# Name and tidy up the bootstrap results.
tidy_boot_df <- bootstrap_result %>%
  purrr::map(.f = ~broom::tidy(.)) %>%
  purrr::set_names(nm = boot_group_nms) %>%
  dplyr::bind_rows(.id = "Boot_Group") %>%
  tibble::as_data_frame() %>%
  tidyr::separate(
    col = Boot_Group,
    into = boot_group_vars,
    sep = "-"
  ) %>%
  dplyr::rename(r = "statistic")

saveRDS(
  tidy_boot_df,
  file =  "./data/processed/predictions/bootstrapped_predictions.RDS"
)
