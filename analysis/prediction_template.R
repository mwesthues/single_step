if (isTRUE(interactive())) {
  .libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.4/"))
}

if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("tidyverse")


rnd_level2_df <- readRDS("./data/derived/predictor_subsets/rnd_level2.RDS")
spec_eta_df <- readRDS("./data/derived/predictor_subsets/spec_eta_df.RDS")

## -- PHENOTYPIC DATA ----------------------------------------------------
inb_traits <- "./data/derived/maizego/tst_pheno_tibble.RDS" %>%
  readRDS() %>%
  pull(Trait) %>%
  unique()
hyb_traits <- "./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  pull(Trait) %>%
  unique()



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

