# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "xtable", "stringr",
               "forcats", "ggthemes", "viridis")
source("./software/prediction_helper_functions.R")

log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

filtered_log_file <- log_file %>%
  filter(Runs == "",
         Data_Type == "Inbred"
  ) %>%
  mutate(Core_Fraction = if_else(nchar(Core_Fraction) == 0,
                                 true = "Full", false = Core_Fraction)
  ) %>%
  mutate(Core_Fraction = if_else(
    (Pred1 == "mrna" & Pred2 == "none"), true = "1.0", false = Core_Fraction
  ))

pred_df <- filtered_log_file %>%
  select(Job_ID) %>%
  flatten_chr() %>%
  map(~readRDS(paste0(prediction_path, ., ".RDS"))) %>%
  bind_rows(.id = "id") %>%
  left_join(
    y = filtered_log_file %>% select(Job_ID, Core_Fraction, Pred1, Pred2),
    by = "Job_ID"
  ) %>%
  unite(col = Predictor, Pred1, Pred2, sep = "-") %>%
  droplevels() %>%
  select(-Job_ID)
  
  

pred_lst <- pred_df %>%
  mutate(
    Reduced = if_else(Core_Fraction == "Full", true = FALSE, false = TRUE)
  ) %>%
  group_by(Trait, Predictor, Core_Fraction, Reduced) %>%
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>%
  ungroup() %>%
  mutate(
    r = round(r, digits = 2),
    CV = round(CV, digits = 4),
    Value = paste0(
      r, " (", CV, ")"
    )
  ) %>%
  select(-r, -CV) %>%
  spread(key = Trait, value = Value) %>%
  split(.$Reduced) %>%
  map(., ~select(., -Reduced)) %>%
  map_at("TRUE", .f = ~mutate(
    ., Predictor = paste0(Predictor, "_", Core_Fraction)
  )) %>%
  map(., ~select(., -Core_Fraction)) %>%
  map(., ~column_to_rownames(., var = "Predictor"))

# Prepare the split data for being displayed inside a LaTeX table.
names(pred_lst)
coverage_subheading <- "Only incomplete predictor:"
attr(pred_lst, "subheadings") <- paste(coverage_subheading, names(pred_lst))
core_caption <- paste(
  "Predictive abilities and corresponding coefficients of variation for",
  "different core sets from the tropical/subtropical material of the maize",
  "diversity panel."
)
xtable_lst <- xtableList(pred_lst,
                         caption = core_caption)
# Function for printing labels in bold font.
bold <- function(x) {
  paste0('{\\bfseries ', x, '}')
}
print.xtableList(xtable_lst, 
                 label = "tbl:core_list",
                 booktabs = TRUE,
                 sanitize.subheadings.function = bold,
                 caption.placement = "top")



# PLOT THE TREND IN R OVER CORE SETS --------------------------------------
g1 <- pred_df %>%
  filter(!Core_Fraction %in% c("Full", "1.0")) %>%
  group_by(Trait, Core_Fraction) %>%
  summarize(r = cor(y, yhat)) %>%
  rename(`Core Fraction` = Core_Fraction) %>%
  ggplot(aes(x = `Core Fraction`, y = r, color = Trait, group = Trait)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  theme_bw(base_size = 10)
plot_name <- "./paper/tables_figures/core_fraction_predictive_ability_trend.pdf"
ggsave(plot = g1, 
       filename = plot_name,
       width = 4,
       height = 4,
       units = "in")
  

