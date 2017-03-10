# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "ggthemes", "viridis",
               "stringr", "forcats")

log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

# Specify a unified predictor order.
pred_order <- c("mrna-none", "ped-none", "snp-none",
                "ped-snp", "ped-mrna", "snp-mrna")

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
  mutate(Predictor = fct_relevel(Predictor, pred_order)) %>%
  droplevels() %>%
  mutate(Trait = fct_recode(Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  mutate(Trait = fct_relevel(Trait, trait_order),
         Predictor = as.factor(Predictor)
  ) %>%
  select(-Job_ID)
  
  
pred_df <- raw_pred_df %>%
  mutate(
    Reduced = ifelse(Predictor == "mrna-none", yes = TRUE, no = FALSE),
    Reduced = ifelse(Core_Fraction == "1.0", yes = TRUE, no = Reduced) 
  ) %>%
  select(-Core_Fraction)



## Distribution of predicted and observed values
legend_labs <- expression(y, hat(y))
g1 <- pred_df %>%
  ggplot(aes(x = y, y = yhat)) +
  geom_point() +
  facet_wrap(Predictor ~ Trait, scales = "free") +
  theme_pander() 

## Predictive ability for seven agronomic traits (rows) and nine predictor sets 
# (columns). The predictor sets are split into the categories 'Reduced' 
# (i.e. 685 hybrids) and 'Full' (i.e. 1,521 hybrids)."
g2 <- pred_df %>%
  group_by(Trait, Predictor, Reduced) %>%
  summarize(r = cor(y, yhat)) %>%
  mutate(Predictor = fct_relevel(Predictor, pred_order)
  ) %>%
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
    legend.position = "bottom",
    legend.key = element_rect(color = "white")
  ) +
  guides(fill = guide_legend(nrow = 3,
                             byrow = TRUE,
                             keywidth = 0.25,
                             keyheight = 0.25,
                             default.unit = "inch",
                             title.position = "top")
  ) +
  ylim(c(0, 1))



coefficient_of_variation <- function(x, y) {
  sqrt(mean(x)) / mean(y)
}


pred_df %>%
  group_by(Trait, Predictor, Reduced) %>%
  summarize(CV = coefficient_of_variation(var_yhat, yhat)) %>%
  ungroup() %>%
  mutate(Predictor = fct_relevel(Predictor, pred_order),
         Trait = fct_relevel(Trait, trait_order)) %>%
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
