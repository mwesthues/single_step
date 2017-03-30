# LOOCV-Prediction
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "dtplyr", "xtable", "stringr",
               "forcats", "ggthemes", "viridis", "gsubfn", "cowplot")
source("./software/prediction_helper_functions.R")

log_path <- "./data/derived/uhoh_maizego_prediction_log.txt"
prediction_path <- "./data/derived/predictions/"
log_file <- fread(log_path)

filtered_log_file <- log_file %>%
  filter(Runs == "",
         Data_Type == "Inbred"
  ) %>%
  mutate(Core_Fraction = if_else(nchar(Core_Fraction) == 0,
                                 true = "Full", false = Core_Fraction)
  ) %>%
  mutate(Core_Fraction = if_else(
    (Pred1 == "mrna" & Pred2 == "none"), true = "1.0", false = Core_Fraction
  ))

raw_pred_df <- filtered_log_file %>%
  select(Job_ID) %>%
  flatten_chr() %>%
  map(~readRDS(paste0(prediction_path, ., ".RDS"))) %>%
  bind_rows(.id = "id") %>%
  left_join(
    y = filtered_log_file %>% select(Job_ID, Core_Fraction, Pred1, Pred2),
    by = "Job_ID"
  ) %>%
  unite(col = Predictor, Pred1, Pred2, sep = "-") %>%
  droplevels() %>%
  select(-Job_ID)
  
pred_df <- raw_pred_df %>%
  mutate(
    Predictor = gsubfn(
      pattern = "\\S+",
      replacement = list(
        "snp-none" = "G",
        "mrna-none" = "T",
        "snp-mrna" = "GT"
      ),
      x = Predictor
    ),
    Trait = gsubfn(
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
  )

# Get the names of genotypes that were imputed based on core genotypes.
#imputed_genotypes <- pred_df %>% 
#  filter(Core_Fraction %in% c("Full", "1.0")) %>% 
#  split(.$Core_Fraction) %>% 
#  map(., ~select(., Geno)) %>% 
#  map(flatten_chr) %>% 
#  .[c("Full", "1.0")] %>% 
#  reduce(setdiff)
 


# PLOT THE TREND IN R OVER CORE SETS --------------------------------------
# Define the top and bottom of the errorbars
limits <- aes(ymax = r + CV, ymin = r - CV)
# Because the bars and errorbars have different widths
# we need to specify how wide the objects we are dodging are
dodge <- position_dodge(width = 0.9)

plot_a_data <- pred_df %>% 
  filter(Core_Fraction %in% c("Full", "1.0")) %>% 
  split(.$Core_Fraction) %>% 
  set_names(c("Core", "Full")) %>% 
#  map_at("Imputed", .f = ~filter(., Geno %in% imputed_genotypes)) %>% 
  bind_rows(.id = "Group") %>% 
  group_by(Trait, Predictor, Group) %>% 
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>%
  ungroup() %>% 
  mutate(
    Predictor = factor(Predictor, levels = c("G", "T", "GT"))
  )

# Predefine the color scheme.
a_colors <- scales::brewer_pal(type = "div", palette = "Spectral")(n = 6) %>% 
  set_names(c("P", "G", "T", "PG", "PT", "GT")) %>% 
  .[names(.) %in% c("G", "T", "GT")]

plot_a <- plot_a_data %>% 
  ggplot(aes(x = Trait, y = r, fill = Predictor)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  facet_grid(. ~ Group) +
  scale_fill_manual(values = a_colors) +
  theme_pander(base_size = 10)
  

plot_b_data <- pred_df %>% 
  filter(Core_Fraction %in% seq(from = 0.1, to = 0.9, by = 0.1)) %>% 
  group_by(Trait, Core_Fraction, Predictor) %>% 
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>%
  ungroup() %>% 
  rename(`Core Fraction` = Core_Fraction)
  
# Because the bars and errorbars have different widths
# we need to specify how wide the objects we are dodging are
dodge_b <- position_dodge(width = 0.3)
plot_b <- plot_b_data %>% 
  ggplot(aes(x = `Core Fraction`, y = r, color = Trait, group = Trait)) +
  geom_line(stat = "identity", position = dodge_b) +
  geom_errorbar(limits, position = dodge_b, width = 0.25) +
  scale_color_tableau() +
  theme_pander(base_size = 10) +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot_c_data <- pred_df %>% 
  filter(Core_Fraction %in% paste0("A", seq_len(3))) %>% 
  group_by(Trait, Core_Fraction, Predictor) %>% 
  summarize(
    r = cor(y, yhat),
    CV = coefficient_of_variation(var_yhat, yhat)
  ) %>%
  ungroup() %>% 
  rename(`Imputed Ancestral Population` = Core_Fraction)

plot_c <- plot_c_data %>% 
  ggplot(aes(x = Trait, y = r, fill = `Imputed Ancestral Population`)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  scale_fill_manual(values = c("#258039", "#F5BE41", "#31A9B8")) +
  theme_pander(base_size = 10) +
  theme(legend.position = "top")


bottom_row <- plot_grid(
  plot_b, plot_c, labels = c("B", "C"), align = 'h', rel_widths = c(1, 1)
)
combi_plot <- plot_grid(
  plot_a, bottom_row, labels = c("A", ""), ncol = 1, rel_heights = c(1, 1.2)
)

ggsave(plot = combi_plot, 
       filename = "./paper/tables_figures/inbred_line_combi_plot.pdf",
       width = 7,
       height = 6,
       units = "in")
  

