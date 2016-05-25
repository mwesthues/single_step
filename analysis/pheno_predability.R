if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes")
fl_nms <- list.files("./data/derived/", pattern = "pheno_prediction")
fl_lst <- lapply(seq_along(fl_nms), FUN = function(i) {
  dat <- readRDS(paste0("./data/derived/", fl_nms[i]))
  dat
})
DT <- rbindlist(fl_lst)
pred_ability <- DT[, .(PredAbility = cor(y, yhat)), 
                   by = .(Phenotype, Predictor)]

ggplot(pred_ability, aes(x = Phenotype, y = PredAbility, fill = Predictor)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.position = "top") +
  xlab("Trait") +
  ylab(bquote(r^{2}))
