# Goal: Collect all individual prediction results in a single data frame.
if(!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse")

pred_dir <- "./data/derived/predictions"
template_names <- list.files(pred_dir, pattern = "template*", full.names = TRUE)

predictions <- template_names %>%
  purrr::map(readRDS) %>%
  dplyr::bind_rows()

saveRDS(predictions, "./data/processed/predictions/predictions.RDS")
