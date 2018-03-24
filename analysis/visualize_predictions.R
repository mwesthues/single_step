if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "dtplyr",
  "data.table",
  "gsubfn",
  "ggthemes",
  "cowplot",
  "svglite"
  )

boot_df <- "./data/processed/predictions/bootstrapped_predictions.RDS" %>%
  readRDS() %>%
  dplyr::group_by(Combi, Trait, Core_Fraction, Predictor) %>%
  dplyr::summarize(
    avg_r = mean(r),
    se = mean(std.error)
  ) %>%
  dplyr::ungroup()


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
limits <- aes(
  ymax = `Predictive Ability` + se,
  ymin = `Predictive Ability` - se
  )

hybrid_plot <- hyb_df %>%
  dplyr::rename(`Predictive Ability` = "avg_r") %>%
  dplyr::mutate(Combi = forcats::fct_recode(
    Combi,
    `Reduced Set` = "Core",
    `Full Set` = "Full"
  )) %>%
  ggplot(aes(x = Trait, y = `Predictive Ability`, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge, color = "black", width = 0.8) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  facet_grid(. ~ Combi) +
  scale_fill_manual(values = mycols) +
  ggthemes::theme_pander(base_size = 10) +
  theme(
    legend.position = "top",
    panel.spacing = unit(2, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(
  plot = hybrid_plot,
  filename = "./paper/tables_figures/pred_ability_hybrid.svg",
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

inbred_scen_a_df <- inbred_df %>%
  dplyr::filter(Core_Fraction == "1") %>%
  dplyr::rename(`Predictive Ability` = "avg_r") %>%
  dplyr::mutate(Combi = forcats::fct_recode(
    Combi,
    `Reduced Set` = "CIA",
    `Full Set` = "FIA"
  )) 

inbred_plot <- inbred_scen_a_df %>%
  ggplot(aes(x = Trait, y = `Predictive Ability`, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge, color = "black") +
  geom_errorbar(limits, position = dodge, width = 0.15) +
  facet_grid(. ~ Combi) +
  scale_fill_manual(values = a_colors) +
  theme_pander(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top"
  )

ggsave(
  plot = inbred_plot,
  filename = "./paper/tables_figures/pred_ability_inbred.svg",
  width = 7,
  height = 4,
  units = "in"
  )





## -- INBREDS CORE_FRACTION != 1 ----------------------------------------------
line_type_df <- inbred_df %>%
  dplyr::filter(Predictor != "GT", Combi == "CIA") %>%
  dplyr::mutate_at(vars(Predictor), funs(as.factor)) %>%
  dplyr::select(-se) %>%
  tidyr::spread(key = Predictor, value = avg_r) %>%
  dplyr::rename(Transcriptome = "T") %>%
  dplyr::rename(Genome = "G") %>%
  dplyr::mutate(Diff = Genome - Transcriptome) %>%
  dplyr::mutate(BetterPredictor = dplyr::case_when(
    Diff > 0.01 ~ "G",
    Diff < -0.01 ~ "T",
    TRUE ~ "None"
    )
  ) %>%
  dplyr::mutate(LineType = dplyr::case_when(
    BetterPredictor == "G" ~ "longdash",
    BetterPredictor == "T" ~ "dotted",
    BetterPredictor == "None" ~ "solid"
    )
  ) %>%
  dplyr::select(-Combi, -Core_Fraction)


core_df <- inbred_df %>%
  dplyr::filter(Combi == "CIA") %>%
  dplyr::mutate(Core_Fraction = dplyr::case_when(
    Predictor == "T" ~ "1",
    Predictor == "G" ~ "0",
    Predictor == "GT" ~ Core_Fraction
    )
  ) %>%
  dplyr::mutate_at(vars(Predictor), funs(as.factor)) %>%
  dplyr::inner_join(y = line_type_df, by = "Trait")

core_plot <- core_df %>%
  dplyr::rename(`Core Fraction (%)` = "Core_Fraction") %>%
  dplyr::rename(`Predictive Ability` = "avg_r") %>%
  ggplot(aes(
    x = `Core Fraction (%)`,
    y = `Predictive Ability`,
    color = Trait,
    group = Trait
    )
  ) +
  geom_line(aes(linetype = LineType), size = 1) +
  geom_point(size = 3, shape = 20, color = "black") +
  ggthemes::scale_color_tableau() +
  ggthemes::theme_pander(base_size = 10) +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_linetype(guide = "none") +
  scale_shape(guide = "none") +
  scale_x_discrete(
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    labels = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  )

ggsave(
  plot = core_plot,
  filename = "./paper/tables_figures/pred_ability_inbred_core_fraction.svg",
  width = 7,
  height = 6,
  units = "in"
  )


