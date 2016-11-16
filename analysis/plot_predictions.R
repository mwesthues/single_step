if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "ggthemes", "viridis",
               "stringr", "corrplot", "forcats")

log_file <- fread("./data/derived/pred_log.txt")

pred_df <- log_file %>%
  .$Job_ID %>%
  as.list() %>%
  map(~readRDS(paste0("./data/derived/predictions/", ., ".RDS"))) %>%
  bind_rows(.id = "id") %>%
  left_join(y = log_file %>% select(Job_ID, Pred1, Pred2, Pred3, Trait),
            by = "Job_ID") %>%
  unite(col = Predictor, Pred1, Pred2, Pred3, sep = "-") %>%
  mutate(Trait = fct_recode(Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  mutate(Trait = fct_relevel(Trait,
                             "DMY", "DMC", "ADF", "FAT", "PRO", "STA", "SUG"),
         Predictor = as.factor(Predictor)) %>%
  select(-Job_ID)


# Update the variable 'Trait' for subsequent joins with 'pred_df'.
log_file <- log_file %>%
  mutate(Trait = fct_recode(Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  mutate(Trait = fct_relevel(Trait,
                             "DMY", "DMC", "ADF", "FAT", "PRO",
                             "STA", "SUG")) %>%
  unite(col = Predictor, Pred1, Pred2, Pred3, sep = "-")


# Compute the average standard deviation for each combination of predictors and 
# traits -> http://stats.stackexchange.com/a/26647
coefficient_of_variation <- function(x, y) {
  sqrt(mean(x)) / mean(y)
}
cv_df <- pred_df %>%
  split(list(.$Trait, .$Predictor)) %>%
  map(~coefficient_of_variation(x = .$var_yhat, y = .$yhat)) %>%
  stack() %>%
  separate(col = ind, into = c("Trait", "Predictor"), sep = "[.]") %>%
  rename(CV = values) %>%
  tbl_dt() %>%
  mutate(Predictor = as.factor(Predictor),
         Trait = as.factor(Trait))


# Combine predictive abilities and coefficients of variation for all predictor
# and trait combinations for plotting.
overview_df <- pred_df %>%
  group_by(Predictor, Trait) %>%
  summarize(r = cor(y, yhat)) %>%
  ungroup() %>%
  as.data.frame() %>%
  left_join(y = cv_df, by = c("Predictor", "Trait")) %>%
  as_tibble() %>%
  mutate(Predictor = as.factor(Predictor),
         Predictor = forcats::fct_reorder(Predictor, CV),
         Trait = fct_relevel(Trait,
                             "DMY", "DMC", "ADF", "FAT", "PRO",
                             "STA", "SUG"),
         Group = ifelse(Predictor %in% c("snp42-none-none",
                                         "ped42-none-none",
                                         "mrna42-none-none"),
                        yes = "Reduced", no = "Full"),
         Group = factor(Group),
  )


# Specify a unified predictor order.
pred_order <- c("ped42-none-none", "snp42-none-none", "mrna42-none-none",
                "ped100-none-none", "snp100-none-none", "ped100-snp77-none",
                "ped100-mrna42-none", "snp100-mrna42-none",
                "ped100-snp77-mrna42")

# Distribution of predicted and observed values.
# Bear in mind that the distribution of the observed values will be different
# for different predictors because some predictors "xx42" and "xx77", where
# 'xx' denotes a predictor, contain less hybrids than there are available 
# (i.e. "xx100").
legend_labs <- expression(y, hat(y))
pred_df %>%
  filter(Predictor != "snp77-none-none") %>%
  gather(key = Value, value = `Phenotypic Value`, y, yhat) %>%
  mutate(Predictor = fct_relevel(Predictor, pred_order)) %>%
  ggplot(aes(x = `Phenotypic Value`, color = Value)) +
  geom_freqpoly() +
  scale_color_manual(values = viridis(length(legend_labs)),
                     labels = legend_labs) +
  facet_grid(Predictor ~ Trait, scales = "free") +
  theme(legend.position = "top") +
  ylab("Frequency")



# Compare the coefficients of variation for all predictor-trait combinations.
overview_df %>%
  filter(Predictor != "snp77-none-none") %>%
  mutate(Predictor = fct_relevel(Predictor, pred_order)) %>%
  ggplot(aes(x = CV, y = Predictor)) +
  geom_segment(aes(yend = Predictor), xend = 0, color = "gray50") +
  geom_point(aes(color = Predictor)) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap(~ Trait, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.7, 0.13)
  ) +
  guides(color = guide_legend(ncol = 2,
                              reverse = FALSE)) +
  xlab("Coefficient of Variation") +
  ylab("Predictor Set") +
  xlim(0, NA)



# Plot the predictive abilities for all predictor-trait combinations 
# separately for the reduced and the full data set, respectively.
overview_df %>%
  as_tibble() %>%
  filter(Predictor != "snp77-none-none") %>%
  mutate(Predictor = fct_relevel(Predictor, pred_order),
         Group = fct_relevel(Group, "Reduced", "Full")
  ) %>%
  ggplot(aes(x = Predictor, y = r)) +
  geom_bar(stat = "identity", aes(fill = Predictor), color = "black") +
  geom_text(aes(label = round(r, digits = 2)), vjust = -0.2) +
  scale_fill_viridis(discrete = TRUE) +
  facet_grid(Trait ~ Group, space = "free", scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank()
  ) +
  ylim(c(0, 1))
