if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "ggthemes")
run_log <- fread("./data/processed/prediction_log/log_list.txt")


# CV1000
cv1000_nms <- run_log[CV != "LOOCV", Job_ID, ]
cv1000 <- lapply(seq_along(cv1000_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", cv1000_nms[i], ".RDS"))
})
cv1000 <- rbindlist(cv1000)
cv1000[, mean(Pred_Ability), by = .(Imputation)]


# LOOCV 
loocv_nms <- run_log[CV == "LOOCV", Job_ID, ]
loocv <- lapply(seq_along(loocv_nms), FUN = function(i) {
  readRDS(paste0("./data/processed/predictions/", loocv_nms[i], ".RDS"))
})
loocv <- rbindlist(loocv)
loocv[, mean(cor(y, yhat)), by = .(Imputation)]

