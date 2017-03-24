# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "xtable", "stringr",
               "forcats", "gsubfn")
# xtable options
options(width = 60)

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
  ) %>%
  select(-Core_Fraction) %>%
  group_by(Trait, Predictor, Group) %>%
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>%
  ungroup()

ext_pred_order <- c(
  "P_Full", "G_Full", "PG_Full", "PT_Full", "GT_Full",
  "P_Core", "G_Core", "T_Core"
)

xtable_pred_df <- pred_df %>%
  mutate(
    r = round(r, digits = 2),
    CV = round(CV, digits = 3),
    Value = paste0(
      r, " (", CV, ")"
    )
  ) %>%
  select(-r, -CV) %>%
  spread(key = Trait, value = Value) %>%
  mutate(
    Predictor = paste0(Predictor, "_", Group),
    Pred_Order = match(Predictor, ext_pred_order)
  ) %>% 
  arrange(Pred_Order) %>% 
  mutate(Predictor = gsub("_", replacement = " ", x = Predictor)) %>% 
  select(-Group, -Pred_Order)
 

hybrid_caption <- paste(
  "Predictive abilities and corresponding coefficients of variation for the",
  "set of maize hybrids. As predictors, pedigree (P), genomic (G),",
  "transcriptomic (T) data and combinations thereof we used."
)
xtable_pred_df <- xtable(
  xtable_pred_df,
  caption = hybrid_caption
)
# Function for printing labels in bold font.
bold <- function(x) {
  paste0('{\\bfseries ', x, '}')
}
print.xtable(
  xtable_pred_df, 
  include.rownames = FALSE,
  label = "tbl:hybrid_list",
  booktabs = TRUE,
  sanitize.subheadings.function = bold,
  caption.placement = "top"
)

