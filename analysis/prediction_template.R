if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("tidyverse", "data.table", "dtplyr")




## -- SCENARIO INFORMATION -----------------------------------------------
main_dir <- "./data/derived/predictor_subsets/"

# Use this function to extract the first letter of 'Extent', 'Scenario'
# and 'Material' to concatenate them.
# This will save memory and reduce the number of variables that need to
# be specified when picking a cluster for predictions.
abbreviate_combi <- function(x) {
  substr(x, start = 1, stop = 1)
}

# Sort all data.tables using this common key, which is essential for
# concatenating them.
common_key <- c(
  "UUID",
  "Combi",
  "Predictor",
  "Core_Fraction",
  "Rnd_Level1",
  "Rnd_Level2",
  "TST_Geno",
  "TRN_Geno"
)

# Data for all genotypes by Material, Extent, Scenario.
geno_df <- readRDS(paste0(main_dir, "geno_df.RDS"))


### Material == "Inbred" && Core_Fraction == 1
inbred_aug <- readRDS(paste0(main_dir, "inbred_trn_df.RDS"))
inbred_pre_eta_df <- readRDS(paste0(main_dir, "inbred_pre_eta_df.RDS"))

inbred_df <- inbred_aug %>%
  tidyr::unite(Predictor, c("Pred1", "Pred2"), remove = TRUE) %>%
  dplyr::mutate(Predictor = gsub("_$", x = Predictor, replacement = "")) %>%
  dplyr::full_join(
    y = inbred_pre_eta_df,
    by = c("Material", "Extent", "Scenario", "Predictor", "G")
   ) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::rename(TRN_Geno = "G") %>%
  data.table::as.data.table()

# Concatenate the initials of 'Extent' (E), 'Material' (M) and 'Scenario' (S)
# for alleviated referencing of combinations and add some dummy variables for
# consistency with data from other combinations.
cols_to_del <- c("Extent", "Material", "Scenario", "E", "M", "S", "Pred1", "Pred2")
inbred_df[
  ,
  c("E", "M", "S") := lapply(.SD, FUN = abbreviate_combi),
  .SDcols = c("Extent", "Material", "Scenario")
 ][
  ,
  Combi := do.call(paste0, args = .SD),
  .SDcols = c("E", "M", "S")
  ][
  ,
  `:=` (TST_Geno = NA_character_,
        TRN_Geno = NA_character_
        ),
  ][
  ,
  (cols_to_del) := NULL,
  ]
data.table::setkeyv(inbred_df, cols = common_key)




### Material == "Hybrid"
keycols = c("Material", "Extent", "Scenario", "TST_Geno")
hyb_lvl2_df <- readRDS(paste0(main_dir, "hybrid_lvl2_df.RDS")) %>%
  data.table::as.data.table()
data.table::setkeyv(hyb_lvl2_df, cols = keycols)

hybrid_pre_eta_df <- readRDS(paste0(main_dir, "hybrid_pre_eta_df.RDS")) %>%
  dplyr::select(-Pred1, -Pred2) %>%
  data.table::as.data.table() %>%
  data.table::setkeyv(cols = keycols)

hybrid_df <- merge(
  x = hyb_lvl2_df,
  y = hybrid_pre_eta_df,
  all = TRUE,
  allow.cartesian = TRUE
)

# Concatenate the initials of 'Extent' (E), 'Material' (M) and 'Scenario' (S)
# for alleviated referencing of combinations and add some dummy variables for
# consistency with data from other combinations.
cols_to_del <- c("Extent", "Material", "Scenario", "E", "M", "S")
hybrid_df[
  ,
  c("E", "M", "S") := lapply(.SD, FUN = abbreviate_combi),
  .SDcols = c("Extent", "Material", "Scenario")
 ][
  ,
  Combi := do.call(paste0, args = .SD),
  .SDcols = c("E", "M", "S")
  ][
  ,
  `:=` (Core_Fraction = NA_real_,
        Rnd_Level1 = NA_real_
        ),
  ][
  ,
  (cols_to_del) := NULL,
  ]
data.table::setkeyv(hybrid_df, cols = common_key)



### Material == "Inbred" && Core_Fraction != 1
cols_to_del <- c("Extent", "Material", "Scenario", "E", "M", "S")
core_df <- readRDS(paste0(main_dir, "core_pre_eta_df.RDS")) %>%
  tidyr::unite(Predictor, c("Pred1", "Pred2"), remove = TRUE) %>%
  dplyr::select(-G) %>%
  dplyr::distinct() %>%
  data.table::as.data.table()

# Concatenate the initials of 'Extent' (E), 'Material' (M) and 'Scenario' (S)
# for alleviated referencing of combinations and add some dummy variables for
# consistency with data from other combinations.
core_df[
  ,
  c("E", "M", "S") := lapply(.SD, FUN = abbreviate_combi),
  .SDcols = c("Extent", "Material", "Scenario")
  ][
  ,
  Combi := do.call(paste0, args = .SD),
  .SDcols = c("E", "M", "S")
  ][
  ,
  `:=` (TST_Geno = NA_character_,
        TRN_Geno = NA_character_
        ),
  ][
  ,
  (cols_to_del) := NULL,
  ]
data.table::setkeyv(core_df, cols = common_key)





## -- COMBINE UNIQUE SCENARIOS
# Keep only scenarios not including training or test set hybrids, respectively.
hybrid_scen_df <- hybrid_df %>%
  .[
  ,
  `:=` (TST_Geno = NA_character_,
        TRN_Geno = NA_character_
        ),
  ] %>%
  data.table::setkey(., NULL) %>%
  unique()

#rm(hybrid_df)
#gc();gc();gc();gc();gc();gc()

core_scen_df <- core_df %>%
  data.table::setkey(., NULL) %>%
  unique()

inbred_scen_df <- inbred_df %>%
  data.table::setkey(., NULL) %>%
  unique()





## -- PHENOTYPIC DATA ----------------------------------------------------
# Augment the scenario data.table by the phenotypic traits, separately for
# inbred lines and hybrids, respectively.
inb_traits <- "./data/derived/maizego/tst_pheno_tibble.RDS" %>%
  readRDS() %>%
  dplyr::distinct(Trait) %>%
  dplyr::mutate(Material = "I")

hyb_traits <- "./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  dplyr::distinct(Trait) %>%
  dplyr::mutate(Material = "H")

trait_df <- dplyr::bind_rows(inb_traits, hyb_traits) %>%
  data.table::as.data.table() %>%
  data.table::setkeyv(cols = "Material")


full_scen_df <- rbindlist(
  list(core_scen_df, inbred_scen_df, hybrid_scen_df),
  use.names = TRUE
)
full_scen_df[, Material := substr(Combi, start = 2, stop = 2), ]
data.table::setkeyv(full_scen_df, cols = "Material")
full_scen_df <- merge(
  trait_df,
  full_scen_df,
  all = TRUE,
  allow.cartesian = TRUE
)
full_scen_df[, Material := NULL, ]

# Rename the UUID-variable to clarify that it pertains to unique ETA objects
# and not to unique rows in this data.table.
data.table::setnames(
  x = full_scen_df,
  old = "UUID",
  new = "ETA_UUID"
)

saveRDS(
  full_scen_df,
  file = "./data/derived/prediction_runs/prediction_template.RDS",
  compress = FALSE
)



## -- SPLIT HYBRID DATA ---------------------------------------------------
# The pile of hybrid data has become excessively large and we do not want
# to load everything (several GB) into memory every time we are making
# predictions using merely a small subset of the data.
# Hence, split the hybrid data based on unique combinations defined above.

# The function `uuid::UUIDgenerate()` generates random character vectors with
# 36 characters, each.
# The first test ensures that the hybrid data have not yet been created and
# only proceeds otherwise.
if (length(list.files(main_dir, pattern = "hybrid_.{36}\\.RDS")) == 0) {

  hybrid_df[
    ,
    Hybrid_UUID := uuid::UUIDgenerate(),
    by = .(UUID, Combi, Predictor, Rnd_Level2)
    ]

  data.table::setnames(
    hybrid_df,
    old = "UUID",
    new = "ETA_UUID"
  )


  hybrid_lst <- split(
    x = hybrid_df,
    by = "Hybrid_UUID",
    drop = TRUE
  )

  rm(hybrid_df)
  gc();gc();gc();gc();gc();gc();gc()


  lapply(seq_along(hybrid_lst), FUN = function(i) {
    name_i <- names(hybrid_lst)[i]
    object_i <- hybrid_lst[[i]]
    saveRDS(
      object_i,
      file = paste0(
        "./data/derived/predictor_subsets/hybrid_",
        name_i,
        ".RDS"
      )
    )
  })

  # Now that we have assigned unique IDs to each hybrid scenario, spit those out
  # so that we can use them as a reference later.
  hybrid_scen_df <- hybrid_lst %>%
    purrr::map(function(x) {
    x[, `:=` (TST_Geno = NA_character_, TRN_Geno = NA_character_), ]
    data.table::setkey(x, NULL)
    unique(x)
    }) %>%
    dplyr::bind_rows()

  saveRDS(
    hybrid_scen_df,
    file = paste0(main_dir, "hybrid_scenario_df.RDS")
  )

} else {
  print("The hybrid subsets have already been created. Don't overwrite them!")
}







## -- ALL COMBINATIONS ----------------------------------------------------
template_precursor <- rnd_level2_df %>%
  select(-TST_Geno, -TRN_Geno, -ind) %>%
  unique() %>%
  full_join(y = spec_eta_df, by = c("Extent", "Material", "Scenario"))

rm(rnd_level2_df, spec_eta_df)

# Get the number of entries per hybrid and per inbred, respectively.
# Then, augment each material group so that every scenario contains the multiple
# traits.
n_hyb <- template_precursor %>%
  filter(Material == "Hybrid") %>%
  nrow()

n_inb <- template_precursor %>%
  filter(Material == "Inbred") %>%
  nrow()

hyb_template <- template_precursor %>%
  filter(Material == "Hybrid") %>%
  slice(rep(1:n(), each = length(hyb_traits))) %>%
  mutate(Trait = rep(hyb_traits, times = n_hyb))

inb_template <- template_precursor %>%
  filter(Material == "Inbred") %>%
  slice(rep(1:n(), each = length(inb_traits))) %>%
  mutate(Trait = rep(inb_traits, times = n_inb))


template_lst <- inb_template %>%
  bind_rows(., hyb_template) %>%
  mutate_at(c("Rnd_Level1", "Rnd_Level2"), as.numeric) %>%
  arrange(
    Material,
    Extent,
    Scenario,
    Pred1,
    Pred2,
    Core_Fraction,
    Rnd_Level1,
    Rnd_Level2,
    Trait
  ) %>%
  rowid_to_column() %>%
  filter(
    Rnd_Level1 %in% seq(from = 0, to = 25, by = 1),
    Rnd_Level2 %in% seq(from = 0, to = 25, by = 1)
  ) %>%
  mutate_at(.vars = c("Extent", "Material", "Scenario"), .funs = as.factor) %>%
  split(list(.$Extent, .$Material, .$Scenario), drop = TRUE) %>%
  keep(~ nrow(.x) != 0)




## -- FUNCTIONS -------------------------------------------------------------
cut_data <- function(job_df, minutes_per_job, minutes_per_node, tolerance) {
  #---
  # Goal: Split a data frame of jobs 'job_df' into intervals corresponding to
  # the number of jobs that can be run on a single node given the amount of time
  # that a single job takes to finish ('minutes_per_job'), the number of minutes
  # that is available on a given node ('minutes_per_node') and a tolerance
  # factor to ensure that a chain of jobs successfully finishes on time.
  #
  ## Parameters
  # job_df: data frame with one job per row
  #
  # minutes_per_job: integer or double scalar with the number of minutes per job
  #
  # minutes_per_node: integer or double scalar with the number of minutes per
  # node
  #
  # tolerance: double scalar to adjust the estimated run time of a single job
  #---

  job_number <- nrow(job_df)
  jobs_per_node <- minutes_per_node / (minutes_per_job * tolerance)
  # The 'node_number' is equivalent to the number of nodes.
  node_number <- ceiling(job_number / jobs_per_node)
  if (isTRUE(node_number != 1)) {
    intervals <- cut(
      x = seq_len(job_number),
      breaks = node_number,
      labels = FALSE
    )
  } else {
    intervals <- rep(1, times = job_number)
  }
  job_df$Interval <- as.character(intervals)
  job_df
}







## -- GENERATE PREDICTION JOB INTERVALS -------------------------------------
# Time per node (hours * minutes/hour * cores)
node_time_per_day <- 24 * 60 * 16
# Time per job in minutes
job_minutes <- c(18, 25, 25, 35, 105, 3200)
names(job_minutes) <- names(template_lst)
# Fractions of the day times number of minutes per node and day
template_node_time <- c(0.5, 1, 1, 1, 1, 2.5) * node_time_per_day
names(template_node_time) <- names(template_lst)

pre_map_lst <- list(
    job_df = template_lst,
    minutes_per_job = job_minutes,
    minutes_per_node = template_node_time,
    tolerance = rep(1.3, times = length(template_lst))
)

# Create short identifiers for each of the six blocks that need to be analyzed
# to ensure that the intervals are unique.
combi_names <- template_lst %>%
  names() %>%
  strsplit(., split = "[.]") %>%
  map(~ substr(., start = 1, stop = 1)) %>%
  map_chr(~ paste(., collapse = ""))

prediction_template <- pmap(. = pre_map_lst, .f = cut_data) %>%
  set_names(., nm = combi_names) %>%
  bind_rows(.id = "Combi") %>%
  mutate(Interval = paste(Combi, Interval, sep = "_"))

# Check how many nodes will be requested for each of the six blocks.
prediction_template %>%
  mutate(Combi = as.factor(Combi)) %>%
  split(.$Combi) %>%
  map(~ pull(., Interval)) %>%
  map(unique) %>%
  map_int(length)


pred_loc <- "./data/derived/prediction_runs/"
saveRDS(prediction_template, paste0(pred_loc, "prediction_template.RDS"))

