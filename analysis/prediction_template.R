if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("tidyverse")


rnd_level2_df <- readRDS("./data/derived/predictor_subsets/rnd_level2.RDS")
spec_eta_df <- readRDS("./data/derived/predictor_subsets/spec_eta_df.RDS")

## -- PHENOTYPIC DATA ----------------------------------------------------
inb_traits <- "./data/derived/maizego/tst_pheno_tibble.RDS" %>%
  readRDS() %>%
  pull(Trait) %>%
  unique()
hyb_traits <- "./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  pull(Trait) %>%
  unique()



## -- ALL COMBINATIONS ----------------------------------------------------
template_precursor <- rnd_level2_df %>%
  select(-TST_Geno, -TRN_Geno, -ind) %>%
  unique() %>%
  full_join(y = spec_eta_df, by = c("Extent", "Material", "Scenario"))

rm(rnd_level2_df, spec_eta_df)

# Get the number of entries per hybrid and per inbred, respectively.
# Then, augment each material group so that every scenario contains the multiple
# traits.
n_hyb <- template_precursor %>%
  filter(Material == "Hybrid") %>%
  nrow()

n_inb <- template_precursor %>%
  filter(Material == "Inbred") %>%
  nrow()
  
hyb_template <- template_precursor %>%
  filter(Material == "Hybrid") %>%
  slice(rep(1:n(), each = length(hyb_traits))) %>%
  mutate(Trait = rep(hyb_traits, times = n_hyb))

inb_template <- template_precursor %>%
  filter(Material == "Inbred") %>%
  slice(rep(1:n(), each = length(inb_traits))) %>%
  mutate(Trait = rep(inb_traits, times = n_inb))


prediction_template <- inb_template %>%
  bind_rows(., hyb_template) %>%
  mutate_at(c("Rnd_Level1", "Rnd_Level2"), as.numeric) %>%
  arrange(
    Material,
    Extent,
    Scenario,
    Pred1,
    Pred2,
    Core_Fraction,
    Rnd_Level1,
    Rnd_Level2,
    Trait
)
  
# Split the template into multiple intervals so that each scenario can be run on
# a separate node on a server.
cut_idx <- prediction_template %>%
  nrow() %>%
  seq_len() %>%
  cut(breaks = 200, labels = FALSE)

prediction_template <- prediction_template %>%
  mutate(Interval = as.character(cut_idx))

pred_loc <- "./data/derived/prediction_runs/"
saveRDS(prediction_template, paste0(pred_loc, "prediction_template.RDS"))
