# Goal: Compare single step prediction with individual predictor results.

## -- PACKAGES ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "dtplyr", "data.table", "broom")


## -- DATA --------------------------------------------------------------------
dat <- "./data/processed/predictions/predictions.RDS" %>%
  readRDS() %>%
  dplyr::mutate(Combi = substr(Interval, start = 1, stop = 3)) %>%
  data.table::as.data.table()


# Find the complement of the full set and the reduced set.
geno_df <- "./data/derived/predictor_subsets/geno_df.RDS" %>%
  readRDS()

hybrid_imp_genos <- geno_df %>%
  dplyr::filter(Material == "Hybrid") %>%
  base::split(x = ., f = list(.$Extent)) %>%
  purrr::map(~dplyr::pull(., G)) %>%
  .[c("Full", "Core")] %>%
  purrr::reduce(setdiff)

# Specify a unified predictor order.
pred_order <- c("P", "G", "T", "PG", "PT", "GT")

# Specify a unique trait order.
trait_order <- c("DMY", "DMC", "ADL", "FAT", "PRO", "STA", "SUG")

# Predefine the color scheme.
mycols <- scales::brewer_pal(type = "div", palette = "Spectral")(n = 6)
names(mycols) <- c("P", "G", "T", "PG", "PT", "GT")

imputed_hybrid_plot <- dat %>%
  tibble::as_data_frame() %>%
  dplyr::filter(Combi == "FHN", TST_Geno %in% hybrid_imp_genos) %>%
  dplyr::group_by(Trait, Predictor) %>%
  dplyr::summarize(r = cor(y, yhat)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Trait = forcats::fct_recode(
    Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  dplyr::mutate(Trait = forcats::fct_relevel(Trait, trait_order)) %>%
  dplyr::mutate(
    Predictor = gsubfn::gsubfn(
      pattern = "\\S+",
      replacement = list(
        "mrna" = "T",
        "ped" = "P",
        "snp" = "G",
        "ped_snp" = "PG",
        "ped_mrna" = "PT",
        "snp_mrna" = "GT"
      ),
      x = Predictor
    ),
    Predictor = forcats::fct_relevel(Predictor, pred_order)
  ) %>%
  dplyr::rename(`Predictive Ability` = "r") %>%
  ggplot(aes(x = Trait, y = `Predictive Ability`, fill = Predictor)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  scale_fill_manual(values = mycols) +
  ggthemes::theme_pander(base_size = 10)


ggsave(
  plot = imputed_hybrid_plot,
  filename = "./tabs_figs/pred_ability_hybrid.pdf",
  width = 7,
  height = 4,
  units = "in"
  )


