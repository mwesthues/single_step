pacman::p_load("data.table", "tidyverse", "dtplyr", "lubridate", "forcats")
log_file <- fread("./data/derived/pred_log.txt")


log_file <- log_file %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S")) %>%
  filter(Date >= as.POSIXct("2016-11-27"))

pred_df <- log_file %>%
  .$Job_ID %>%
  as.list() %>%
  map(~readRDS(paste0("./data/derived/predictions/", ., ".RDS"))) %>%
  bind_rows(.id = "id")

pred_df %>%
  left_join(y = log_file %>% select(Job_ID, Pred1, Pred2, Pred3, Trait, Iter),
            by = "Job_ID") %>%
  unite(col = Predictor, Pred1, Pred2, Pred3, sep = "-") %>%
  select(-Job_ID, -Iter) %>%
  group_by(Predictor, Trait) %>%
  summarize(PredAbility = cor(y, yhat)) %>%
  as.data.frame()



readRDS("./data/derived/predictions/10094280.RDS") %>%
  plot(.$y, .$yhat)
  .[, .(cor(y, yhat)), ]
