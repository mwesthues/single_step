if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes", "dplyr", "forcats", 
               "viridis")
run_log <- fread("./data/processed/prediction_log/log_list.txt")
geno <- readRDS("./data/processed/common_genotypes.RDS")
dent_imp_nms <- setdiff(geno$Dent$snp, geno$Dent$mrna)
flint_imp_nms <- setdiff(geno$Flint$snp, geno$Flint$mrna)

# CV1000
cv1000_nms <- run_log[CV != "LOOCV", Job_ID, ]
cv1000_lst <- lapply(seq_along(cv1000_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", cv1000_nms[i], ".RDS"))
})
cv1000 <- rbindlist(cv1000_lst)
cv1000_avg <- cv1000[Set == "T0", 
                     .(Pred_Ability = cor(y, yHat)),
                     by = .(Imputation, Trait)]
cv1000_avg[, CV_Mode := "CV1000", ]


# LOOCV 
loocv_nms <- run_log[CV == "LOOCV", Job_ID, ]
loocv <- lapply(seq_along(loocv_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", loocv_nms[i], ".RDS"))
})
loocv <- rbindlist(loocv)
loocv_avg <- loocv[, .(Pred_Ability = cor(y, yhat)), 
                   by = .(Imputation, Trait)]
loocv_avg[, CV_Mode := "LOOCV", ]

# Concatenation and transformation
cv_avg <- rbind(loocv_avg, cv1000_avg)
cv_avg <- cv_avg %>%
  mutate(Trait = fct_recode(Trait, 
    "FAT" = "FETT",
    "PRO" = "RPR",
    "DMY" = "GTM",
    "DMC" = "GTS"
  )) %>%
  mutate(Trait = fct_relevel(Trait,
    c("DMY", "DMC", "ADF", "FAT", "PRO", "STA", "XZ")
  )) %>%
  rename(ssBLUP = Imputation)

# Plot
g1 <- cv_avg %>%
  ggplot(aes(x = Trait, y = Pred_Ability, fill = ssBLUP)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_grid(CV_Mode ~ .) +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "bottom") +
  ylab("Predictive ability (r)") +
  labs(title = "single-step BLUP works best with LOOCV")
ggsave(plot = g1, 
       file = "./tabs_figs/Pred_Ability.pdf", width = 6, height = 6)


