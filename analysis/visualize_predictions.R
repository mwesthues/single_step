if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "dtplyr",
  "data.table",
  "gsubfn",
  "ggthemes",
  "cowplot"
  )

boot_df <- "./data/processed/predictions/bootstrapped_predictions.RDS" %>%
  readRDS()


## -- HYBRIDS -----------------------------------------------------------------
# Specify a unified predictor order.
pred_order <- c("P", "G", "T", "PG", "PT", "GT")

# Specify a unique trait order.
trait_order <- c("DMY", "DMC", "ADL", "FAT", "PRO", "STA", "SUG")


hyb_df <- boot_df %>%
  dplyr::filter(Combi %in% c("CHN", "FHN")) %>%
  dplyr::mutate(Trait = forcats::fct_recode(
    Trait,
    "DMY" = "GTM",
    "DMC" = "GTS",
    "FAT" = "FETT",
    "PRO" = "RPR",
    "SUG" = "XZ"
  )) %>%
  dplyr::mutate(Combi = forcats::fct_recode(
    Combi,
    "Core" = "CHN",
    "Full" = "FHN"
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
  )

# Because the bars and errorbars have different widths
# we need to specify how wide the objects we are dodging are
dodge <- position_dodge(width = 0.9)

# Predefine the color scheme.
mycols <- scales::brewer_pal(type = "div", palette = "Spectral")(n = 6)
names(mycols) <- c("P", "G", "T", "PG", "PT", "GT")

hybrid_plot <- hyb_df %>%
  ggplot(aes(x = Trait, y = r, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge, color = "black") +
  facet_grid(. ~ Combi) +
  scale_fill_manual(values = mycols) +
  ggthemes::theme_pander(base_size = 10)

ggsave(
  plot = hybrid_plot,
  filename = "./paper/tables_figures/pred_ability_hybrid.pdf",
  width = 7,
  height = 4,
  units = "in"
  )






## -- INBREDS CORE_FRACTION == 1 ----------------------------------------------
inbred_df <- boot_df %>%
  dplyr::filter(!Combi %in% c("CHN", "FHN")) %>%
  dplyr::mutate(
    Predictor = gsubfn::gsubfn(
      pattern = "\\S+",
      replacement = list(
        "snp" = "G",
        "mrna" = "T",
        "snp_mrna" = "GT"
      ),
      x = Predictor
    ),
    Trait = gsubfn::gsubfn(
      pattern = "\\S+",
      replacement = list(
        "Plantheight" = "PH",
        "Eardiameter" = "ED",
        "100grainweight" = "GW",
        "cobweight" = "CW",
        "Kernelwidth" = "KW",
        "Silkingtime" = "DS"
      ),
      x = Trait
    )
  ) %>%
  dplyr::mutate(Combi = forcats::fct_recode(
    Combi,
    "Core" = "CHN",
    "Full" = "FHN"
  ))

# Predefine the color scheme.
# Use also the levels of predictors in the hybrid data set to ensure that,
# regardless of the material, the same color is picked for the same predictor.
a_colors <- scales::brewer_pal(type = "div", palette = "Spectral")(n = 6) %>% 
  set_names(c("P", "G", "T", "PG", "PT", "GT")) %>% 
  .[names(.) %in% c("G", "T", "GT")]


inbred_plot <- inbred_df %>%
  dplyr::filter(Core_Fraction == "1") %>%
  ggplot(aes(x = Trait, y = r, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge) +
  facet_grid(. ~ Combi) +
  scale_fill_manual(values = a_colors) +
  theme_pander(base_size = 10)

ggsave(
  plot = inbred_plot,
  filename = "./paper/tables_figures/pred_ability_inbred.pdf",
  width = 7,
  height = 6,
  units = "in"
  )





## -- INBREDS CORE_FRACTION != 1 ----------------------------------------------
core_plot <- inbred_df %>%
  dplyr::filter(Predictor == "GT", Combi %in% c("CIA", "CIB")) %>%
  dplyr::mutate_at(vars(Predictor), funs(as.factor)) %>%
  dplyr::rename(`Core Fraction` = "Core_Fraction") %>%
  ggplot(aes(x = `Core Fraction`, y = r, color = Trait, group = Trait)) +
  geom_line(stat = "identity", position = position_dodge(width = 0.3)) +
  facet_grid(~ Combi) +
  ggthemes::scale_color_tableau() +
  ggthemes::theme_pander(base_size = 10) +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(
  plot = core_plot,
  filename = "./paper/tables_figures/pred_ability_inbred_core_fraction.pdf",
  width = 7,
  height = 6,
  units = "in"
  )


