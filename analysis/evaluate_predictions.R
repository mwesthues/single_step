if (isTRUE(interactive())) {
  .libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")
}


if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "broom")

# load the data
dat <- "./data/processed/predictions/predictions.RDS" %>%
  readRDS() %>%
  filter(complete.cases(.))


boot_iter <- 2

boot_pred <- dat %>%
  dplyr::rename(y = "Observed") %>%
  tidyr::unite(PredCombi, Pred1, Pred2, sep = "_") %>%
  tidyr::nest(
    -Extent,
    -Material, 
    -Scenario,
    -PredCombi,
    -Rnd_Level1,
    -Rnd_Level2,
    -Core_Fraction,
    -Trait
  ) %>%
  dplyr::mutate(cors_boot = purrr::map(
    data,
    ~broom::bootstrap(., boot_iter) %>% do(broom::tidy(cor(.$y, .$yhat))))
  ) %>%
  tidyr::unnest(cors_boot) %>%
  dplyr::group_by(
    Extent,
    Material,
    Scenario,
    PredCombi,
    Core_Fraction,
    Trait
  ) %>%
  dplyr::summarize(r = mean(x), se = sd(x))

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

