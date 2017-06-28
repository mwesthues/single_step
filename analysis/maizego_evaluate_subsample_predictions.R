if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "broom")

# load the log file
log_file <- "./data/derived/uhoh_maizego_prediction_log.txt" %>%
  read_tsv() %>%
  mutate(Job_ID = as.character(Job_ID)) %>%
  filter(grepl("2017-06-27", x = Start_Time))

# concatenate log file with prediction results for complete information
pred_data <- log_file %>%
  pull(Job_ID) %>%
  map(., .f = ~readRDS(paste0("./data/derived/predictions/", ., ".RDS"))) %>%
  map(as_data_frame) %>%
  bind_rows() %>%
  left_join(x = ., y = log_file, by = "Job_ID")

boot_pred <- pred_data %>%
  nest(-Core_Set, -Core_run, -Loocv_run, -Trait) %>%
  mutate(cors_boot = map(data, ~broom::bootstrap(., 100) %>%
                       do(tidy(cor(.$y, .$yhat))))) %>%
  unnest(cors_boot) %>%
  group_by(Core_Set, Trait) %>%
  summarize(r = mean(x), se = sd(x))


saveRDS(boot_pred, "./data/derived/maizego/bootstrap_110_genotypes.RDS")


# bootstrap plot
boot_pred %>%
  ungroup() %>%
  mutate(upper = r + se, lower = r - se) %>%
  ggplot(aes(Core_Set %>% as.character() %>% as.factor(), y = r)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymax = upper, ymin = lower)) +
  facet_wrap(~ Trait)


pred_data %>%
  group_by(Core_Set, Trait) %>%
  summarize(r = cor(y, yhat)) %>%
  ggplot(aes(Core_Set %>% as.character() %>% as.factor(), y = r)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Trait)

