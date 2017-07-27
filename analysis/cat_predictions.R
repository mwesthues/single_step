# Goal: Collect all individual prediction results in a single data frame.
.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")
if(!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse")

pred_dir <- "./data/derived/predictions"
template_names <- list.files(pred_dir, pattern = "template*", full.names = TRUE)

predictions <- template_names %>%
  purrr::map(readRDS) %>%
  dplyr::bind_rows()

saveRDS(predictions, "./data/processed/predictions/predictions.RDS")
