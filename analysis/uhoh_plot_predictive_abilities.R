# Goal: Plot results of hybrid-based predictive abilities.
# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "stringr", "forcats",
               "gsubfn", "viridis", "ggthemes")

source("./software/prediction_helper_functions.R")


log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

# Specify a unified predictor order.
pred_order <- c("T", "P", "G", "PG", "PT", "GT")

# Specify a unique trait order.
trait_order <- c("DMY", "DMC", "FAT", "PRO", "STA", "SUG")

raw_pred_df <- log_file %>%
  filter(Runs == "",
         Data_Type == "Hybrid") %>%
  select(Job_ID) %>%
  flatten_chr() %>%
  map(~readRDS(paste0(prediction_path, ., ".RDS"))) %>%
  bind_rows(.id = "id") %>%
  left_join(
    y = log_file %>% select(Job_ID, Core_Fraction, Pred1, Pred2), by = "Job_ID"
    ) %>%
  unite(col = Predictor, Pred1, Pred2, sep = "-") %>%
  filter(Trait != "ADL") %>%
  droplevels() %>%
  mutate(Trait = fct_recode(Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  mutate(
    Trait = fct_relevel(Trait, trait_order)
  ) %>%
  select(-Job_ID)

# Replace long predictor names by short, one-character versions.
raw_pred_df <- raw_pred_df %>%
  mutate(
    Predictor = gsubfn(
      pattern = "\\S+", 
      replacement = list(
        "mrna-none" = "T",
        "ped-none" = "P",
        "snp-none" = "G",
        "ped-snp" = "PG",
        "ped-mrna" = "PT",
        "snp-mrna" = "GT"
      ),
      x = Predictor
    ),
    Predictor = fct_relevel(Predictor, pred_order)
  )
 
 
pred_df <- raw_pred_df %>%
  mutate(
    Group = ifelse(Predictor == "T", yes = "Core", no = "Full"),
    Group = ifelse(Core_Fraction == "1.0", yes = "Core", no = Group) 
  )
  
# Get the names of genotypes that were imputed based on core genotypes.
imputed_genotypes <- pred_df %>% 
  split(.$Group) %>% 
  map(., ~select(., Geno)) %>% 
  map(flatten_chr) %>% 
  .[c("Full", "Core")] %>% 
  reduce(setdiff)
  
compute_r_and_cv <- function(x) {
  x %>% 
    group_by(Trait, Predictor) %>% 
    summarize(
      r = cor(y, yhat),
      CV = coefficient_of_variation(var_yhat, yhat)
    ) %>%
    ungroup()
}

ext_pred_order <- c(
  "P_Imputed", "G_Imputed", "PG_Imputed", "PT_Imputed", "GT_Imputed",
  "P_Core", "G_Core", "T_Core"
)

plot_data <- pred_df %>% 
  split(.$Group) %>% 
  map_at("Full", .f = ~filter(., Geno %in% imputed_genotypes)) %>% 
  map(compute_r_and_cv) %>% 
  set_names(c("Core", "Imputed")) %>% 
  bind_rows(.id = "Group") %>% 
  unite(col = Exact_Predictor, Predictor, Group, sep = "_", remove = FALSE) %>% 
  mutate(Pred_Order = match(Exact_Predictor, ext_pred_order)) %>% 
  arrange(Pred_Order)

# Define the top and bottom of the errorbars
limits <- aes(ymax = r + CV, ymin = r - CV)
# Because the bars and errorbars have different widths
# we need to specify how wide the objects we are dodging are
dodge <- position_dodge(width = 0.9)

plot_data %>% 
  ggplot(aes(x = Trait, y = r, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  facet_grid(. ~ Group) +
  scale_fill_viridis(discrete = TRUE) +
  theme_base()
  






