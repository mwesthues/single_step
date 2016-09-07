# Goal: Check what the predictive abilities are when mRNA values were missing 
# for no, one or both parents.

if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes", "dplyr", "forcats", 
               "viridis")
run_log <- fread("./data/processed/prediction_log/log_list.txt")
geno <- readRDS("./data/processed/common_genotypes.RDS")
dent_imp_nms <- setdiff(geno$Dent$snp, geno$Dent$mrna)
flint_imp_nms <- setdiff(geno$Flint$snp, geno$Flint$mrna)


# LOAD DATA ---------------------------------------------------------------
# CV1000
cv1000_nms <- run_log[CV != "LOOCV", Job_ID, ]
cv1000_lst <- lapply(seq_along(cv1000_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", cv1000_nms[i], ".RDS"))
})
cv1000 <- rbindlist(cv1000_lst)

# LOOCV
loocv_nms <- run_log[CV == "LOOCV", Job_ID, ]
loocv <- lapply(seq_along(loocv_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", loocv_nms[i], ".RDS"))
})
loocv <- rbindlist(loocv)


# COMPUTE PREDICTIVE ABILITY ----------------------------------------------
# CV1000
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  rename(y = EST) %>%
  select(G, y, Trait)
cv1000_avg <- cv1000 %>%
  filter(Imputation == "TRUE", Set == "T0") %>%
  select(-Imputation) %>%
  # Get the hybrid names from the phenotypic data, which will be necessary to 
  # determine the parents for which mRNA data where imputed.
  left_join(pheno, by = c("y", "Trait")) %>%
  separate(G, into = c("DF", "Dent", "Flint"), sep = "_") %>%
  select(-DF) %>%
  # Determine for which parents mRNA data were imputed.
  mutate(Dent_Imputed = ifelse(Dent %in% dent_imp_nms, yes = 1, no = 0),
         Flint_Imputed = ifelse(Flint %in% flint_imp_nms, yes = 1, no = 0),
         Number_Imputed = Dent_Imputed + Flint_Imputed) %>%
  group_by(Trait, Number_Imputed) %>%
  summarize(Pred_Ability = cor(y, yHat)) %>%
  ungroup() %>%
  mutate(CV_Mode = "CV1000")

# LOOCV
loocv_avg <- loocv %>%
  filter(Imputation == "TRUE") %>%
  select(-Imputation) %>%
  rename(yHat = yhat) %>%
  mutate(Dent_Imputed = ifelse(Mother %in% dent_imp_nms, 
                               yes = TRUE, no = FALSE),
         Flint_Imputed = ifelse(Father %in% flint_imp_nms, 
                                yes = TRUE, no = FALSE),
         Number_Imputed = Dent_Imputed + Flint_Imputed) %>%
  group_by(Trait, Number_Imputed) %>%
  summarize(Pred_Ability = cor(y, yHat)) %>%
  ungroup() %>%
  mutate(CV_Mode = "LOOCV")
 

# ORDER DATA --------------------------------------------------------------
preds <- rbind(cv1000_avg, loocv_avg)
preds <- preds %>%
  mutate(Trait = as.factor(Trait),
         Number_Imputed = as.factor(as.character(Number_Imputed))) %>%
  mutate(Trait = fct_recode(Trait,
                            "FAT" = "FETT",
                            "PRO" = "RPR",
                            "DMY" = "GTM",
                            "DMC" = "GTS")) %>%
  mutate(Trait = fct_relevel(Trait, 
                             c("DMY", "DMC", "ADF", "FAT", "PRO",
                               "STA", "XZ")),
         Number_Imputed = fct_relevel(Number_Imputed,
                                      c(0, 1, 2)))
# Plot predictive abilities for CV1000 depending on the number of parents for 
# which mRNA data were imputed.
g1 <- preds %>%
  ggplot(aes(x = Trait, y = Pred_Ability, color = Number_Imputed)) +
  geom_point(size = 3) +
  facet_wrap(~ CV_Mode) +
  theme_light() +
  scale_color_viridis(discrete = TRUE,
                      guide = guide_legend(
    title = "Number of parents for which mRNA values were imputed"
  )) +
  labs(title = "Predictive abilities higher for imputed parents") +
  ylab("Predictive ability (r)") +
  theme(legend.position = "bottom")
g1
ggsave(plot = g1, 
       file = "./tabs_figs/pred_ability_imputed_parents.pdf",
       width = 6, height = 6)
