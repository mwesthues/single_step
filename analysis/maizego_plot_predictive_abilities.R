# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "ggthemes", "viridis",
               "stringr", "forcats")

log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

filtered_log_file <- log_file %>%
  filter(Runs == "",
         Data_Type == "Inbred") %>%
  mutate(Core_Fraction = if_else(nchar(Core_Fraction) == 0, 
                                 true = "99", false = Core_Fraction)
  ) %>%
  mutate(Core_Fraction = as.numeric(Core_Fraction)) %>%
  filter(Core_Fraction >= 1) %>%
  mutate(Reduced = if_else(Core_Fraction == 1, true = "yes", false = "no"),
         Reduced = if_else(Pred1 == "mrna", true = "yes", false = Reduced)
  )

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
  
  

## Distribution of predicted and observed values
legend_labs <- expression(y, hat(y))
g1 <- pred_df %>%
  unite(col = "Unique_Predictor", Reduced, Predictor, sep = "_") %>%
  ggplot(aes(x = y, y = yhat)) +
  geom_point() +
  facet_wrap(Unique_Predictor ~ Trait, scales = "free") +
  theme_pander() 

## Predictive ability for seven agronomic traits (rows) and nine predictor sets 
# (columns). The predictor sets are split into the categories 'Reduced' 
# (i.e. 685 hybrids) and 'Full' (i.e. 1,521 hybrids)."
g2 <- pred_df %>%
  group_by(Trait, Predictor, Reduced) %>%
  summarize(r = cor(y, yhat)) %>%
  ggplot(aes(x = Predictor, y = r)) +
  geom_bar(stat = "identity", aes(fill = Predictor), color = "black") +
  geom_text(aes(label = round(r, digits = 3)), vjust = -0.2) +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_grid(Trait ~ Reduced, space = "free", scales = "free_x") +
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
g2



coefficient_of_variation <- function(x, y) {
  sqrt(mean(x)) / mean(y)
}


pred_df %>%
  group_by(Trait, Predictor, Reduced) %>%
  summarize(CV = coefficient_of_variation(var_yhat, yhat)) %>%
  ungroup() %>%
  ggplot(aes(x = CV, y = Predictor)) +
  geom_segment(aes(yend = Predictor), xend = 0, color = "gray50") +
  geom_point(aes(color = Predictor), size = 3) +
  scale_color_viridis(discrete = TRUE, option = "inferno") +
  facet_grid(Reduced ~ Trait, scales = "free") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top"
  ) +
  guides(color = guide_legend(ncol = 2,
                              reverse = FALSE)) +
  xlab("Coefficient of Variation") +
  ylab("Predictor Set") +
  xlim(0, NA)
