# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "stringr", "forcats",
               "xtable", "gsubfn")
source("./software/prediction_helper_functions.R")

# xtable options
options(width = 60)

log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

# Specify which prediction results shall be loaded.
filtered_log_file <- log_file %>%
  filter(Runs == "",
         Data_Type == "Inbred") %>%
  mutate(Core_Fraction = if_else(nchar(Core_Fraction) == 0, 
                                 true = "99", false = Core_Fraction)
  ) %>%
  mutate(Core_Fraction = as.numeric(Core_Fraction)) %>%
  filter(Core_Fraction >= 1) %>%
  mutate(Reduced = if_else(Core_Fraction == 1, true = TRUE, false = FALSE),
         Reduced = if_else(Pred1 == "mrna", true = TRUE, false = Reduced)
  )

# Read the results from the prediction models based on the log files and 
# concatenate them in a common data frame.
pred_df <- filtered_log_file %>%
  select(Job_ID) %>%
  flatten_chr() %>%
  map(~readRDS(paste0(prediction_path, ., ".RDS"))) %>%
  bind_rows(.id = "id") %>%
  left_join(
    y = filtered_log_file %>% select(Job_ID, Reduced, Pred1, Pred2),
    by = "Job_ID"
  ) %>%
  unite(col = Predictor, Pred1, Pred2, sep = "-") %>%
  droplevels() %>%
  select(-Job_ID)

# Replace long predictor names by short, one-character versions.
pred_df <- pred_df %>%
  mutate(
    Predictor = gsubfn(
      pattern = "\\S+", 
      replacement = list(
        "snp-none" = "G",
        "mrna-none" = "T",
        "snp-mrna" = "GT"
      ),
      x = Predictor
    )
  )


# For each combination of trait, genotype coverage and predictor, compute the 
# predictive ability and the coefficient of variation of predictions.
# Combine the two statistics in a single variable where the CV-values are 
# printed in parentheses.
# Then split the data into a list with a separate data frame for each level of
# genotype coverage (i.e., only genotypes covered by both predictors vs. all
# genotypes).
pred_lst <- pred_df %>%
  group_by(Trait, Reduced, Predictor) %>%
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>% 
  ungroup() %>%
  mutate(
    r = round(r, digits = 2),
    CV = round(CV, digits = 2),
    Value = paste0(
      r, " (", CV, ")"
  )) %>%
  select(-r, -CV) %>%
  spread(key = Trait, value = Value) %>%
  split(.$Reduced) %>%
  map(., ~select(., -Reduced)) %>%
  map(., ~column_to_rownames(., var = "Predictor"))

# Prepare the split data for being displayed inside a LaTeX table.
names(pred_lst)
coverage_subheading <- "Only incomplete predictor:"
attr(pred_lst, "subheadings") <- paste(coverage_subheading, names(pred_lst))
inbred_caption <- paste(
  "Predictive abilities and corresponding coefficients of variation for the",
  "set of 211 tropical/subtropical lines from the maize diversity panel."
)
xtable_lst <- xtableList(pred_lst,
                         caption = inbred_caption)
# Function for printing labels in bold font.
bold <- function(x) {
  paste0('{\\bfseries ', x, '}')
}
print.xtableList(xtable_lst, 
                 label = "tbl:inbred_list",
                 booktabs = TRUE,
                 sanitize.subheadings.function = bold,
                 caption.placement = "top")
