if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes")

fl_nms <- list.files("./data/derived/", pattern = "pheno_prediction")
fl_lst <- lapply(seq_along(fl_nms), FUN = function(i) {
  fl_nm <- fl_nms[i]
  mod_pos <- regexpr("(?<=Model=).*(?=_VCOV)", text = fl_nm, perl = TRUE)
  mod_nm <- substring(text = fl_nm, first = mod_pos,
                      last = mod_pos + attr(mod_pos, "match.length") - 1)
  vcov_pos <- regexpr("(?<=VCOV=).*(?=_Iter)", text = fl_nm, perl = TRUE)
  vcov_nm <- substring(text = fl_nm, first = vcov_pos,
                       last = vcov_pos + attr(vcov_pos, "match.length") - 1)
  iter_pos <- regexpr("(?<=Iter=).*(?=_SnpFilter)", text = fl_nm, perl = TRUE)
  iter_nm <- substring(text = fl_nm, first = iter_pos,
                       last = iter_pos + attr(iter_pos, "match.length") - 1)
  fltr_pos <- regexpr("(?<=SnpFilter=).*(?=_Imputed)", text = fl_nm, 
                      perl = TRUE)
  fltr_nm <- substring(text = fl_nm, first = fltr_pos,
                       last = fltr_pos + attr(fltr_pos, "match.length") - 1)
  imp_pos <- regexpr("(?<=Imputed=).*(?=_Comparison)", text = fl_nm, 
                      perl = TRUE)
  imp_nm <- substring(text = fl_nm, first = imp_pos,
                      last = imp_pos + attr(imp_pos, "match.length") - 1)
  cmp_pos <- regexpr("(?<=Comparison=).*(?=.RDS)", text = fl_nm, 
                      perl = TRUE)
  cmp_nm <- substring(text = fl_nm, first = cmp_pos,
                      last = cmp_pos + attr(cmp_pos, "match.length") - 1)
  dat <- readRDS(paste0("./data/derived/", fl_nm))
  dat[, `:=`(Model = mod_nm,
             VCOV = vcov_nm,
             Iter = iter_nm,
             SNP_Filter = fltr_nm,
             Imputed = imp_nm,
             Comparison = cmp_nm), ]

  dat
})
DT <- rbindlist(fl_lst)
pred_ability <- DT[, .(PredAbility = cor(y, yhat)), 
                   by = .(Phenotype, Predictor, Model, VCOV, Iter, SNP_Filter,
                          Imputed, Comparison)]
pred_ability <- pred_ability[!(Predictor == "mrna" & Comparison == "mrna"), , ]


p1 <- ggplot(pred_ability, aes(x = Phenotype, y = PredAbility, 
                               fill = Predictor)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_pander(base_size = 8) +
  facet_grid(. ~ Imputed) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  xlab("Trait") +
  ylab(bquote(r))
ggsave("./tabs_figs/pheno_predability-barplot.pdf", width = 4, height = 4, 
       plot = p1)

