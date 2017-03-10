# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "ggthemes", "viridis",
               "stringr", "forcats")

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
  
pred_df %>%
  group_by(Predictor, Core_Fraction) %>%
  count()
  

## Predictive ability for seven agronomic traits (rows) and nine predictor sets 
# (columns). The predictor sets are split into the categories 'Reduced' 
# (i.e. 685 hybrids) and 'Full' (i.e. 1,521 hybrids)."
g1 <- pred_df %>%
  group_by(Trait, Predictor, Core_Fraction) %>%
  summarize(r = cor(y, yhat)) %>%
  ggplot(aes(x = Predictor, y = r)) +
  geom_bar(stat = "identity", aes(fill = Predictor), color = "black") +
  geom_text(aes(label = round(r, digits = 3)), vjust = -0.2) +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_grid(Trait ~ Core_Fraction, space = "free", scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_rect(color = "white")
  ) +
  guides(fill = guide_legend(nrow = 1,
                             byrow = TRUE,
                             keywidth = 0.25,
                             keyheight = 0.25,
                             default.unit = "inch",
                             title.position = "top")
  ) +
  ylim(c(0, 1))
g1



