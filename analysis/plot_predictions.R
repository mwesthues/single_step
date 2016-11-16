if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "ggthemes", "viridis",
               "stringr", "corrplot", "forcats")

log_file <- fread("./data/derived/pred_log.txt")

#log_file %>%
#  filter(inverted_grepl(pattern = "2016-11-14", x = Date),
#         Runs == 1521,
#         Cores != 16) %>%
#  .$Job_ID %>%
#  as.list() %>%
#  map(~file.remove(paste0("./data/derived/predictions/", ., ".RDS")))

inverted_grepl <- compose(`!`, grepl)
pred_df <- log_file %>%
  filter(inverted_grepl(pattern = "2016-11-14", x = Date),
         Runs == 1521,
         Cores == 16,
         Job_ID != "10085976") %>%
  distinct(Pred1, Pred2, Pred3, Trait, .keep_all = TRUE) %>%
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
                             "DMY", "DMC", "ADF", "FAT", "PRO", "STA", "SUG"))


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
                             "DMY", "DMC", "ADF", "FAT", "PRO", "STA", "SUG"))


# Compute the average standard deviation for each combination of predictors and 
# traits -> http://stats.stackexchange.com/a/26647
coefficient_of_variation <- function(x, y) {
  sqrt(mean(x)) / mean(y)
}
cv_df <- pred_df %>%
  split(.$Job_ID) %>%
  map(~coefficient_of_variation(x = .$var_yhat, y = .$yhat)) %>%
  stack() %>%
  rename(Job_ID = ind,
         CV = values) %>%
  mutate(Job_ID = as.character(Job_ID)) %>%
  left_join(y = log_file %>% select(Job_ID, Pred1, Pred2, Pred3, Trait),
            by = "Job_ID") %>%
  unite(col = Predictor, Pred1, Pred2, Pred3, sep = "-")

# Combine predictive abilities and coefficients of variation for all predictor
# and trait combinations for plotting.
overview_df <- pred_df %>%
  group_by(Predictor, Trait) %>%
  summarize(r = cor(y, yhat)) %>%
  ungroup() %>%
  as.data.frame() %>%
  left_join(y = cv_df, by = c("Predictor", "Trait"))

# Compare the coefficients of variation for all predictor-trait combinations.
overview_df %>%
  as_tibble() %>%
  mutate(Predictor = as.factor(Predictor),
         Predictor = forcats::fct_reorder(Predictor, CV, .desc = TRUE)) %>%
  ggplot(aes(x = CV, y = Predictor)) +
  geom_segment(aes(yend = Predictor), xend = 0, color = "gray50") +
  geom_point(aes(color = Trait)) +
  scale_color_viridis(guide = FALSE, discrete = TRUE) +
  facet_wrap(~ Trait, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  xlab("Coefficient of Variation") +
  ylab("Predictor Set") +
  xlim(0, NA)



# Distribution of predicted and observed values.
# Bear in mind that the distribution of the observed values will be different
# for different predictors because some predictors "xx42" and "xx77", where
# 'xx' denotes a predictor, contain less hybrids than there are available 
# (i.e. "xx100").
legend_labs <- expression(y, hat(y))
pred_df %>%
  gather(key = Value, value = `Phenotypic Value`, y, yhat) %>%
  ggplot(aes(x = `Phenotypic Value`, color = Value)) +
  geom_freqpoly() +
  scale_color_manual(values = viridis(length(legend_labs)),
                     labels = legend_labs) +
  facet_grid(Predictor ~ Trait, scales = "free") +
  theme(legend.position = "top") +
  ylab("Frequency")


#log_file <- log_file %>%
#  mutate(Runs = str_replace(Runs, pattern = "^$", replacement = "1521")) %>%
#  filter(!Runs %in% c("1-500", "501-1000"))
#write.table(x = log_file,
#            file = "./data/derived/pred_log.txt",
#            sep = "\t",
#            row.names = FALSE)



