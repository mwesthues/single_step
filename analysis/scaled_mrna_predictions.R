if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("BGLR", "data.table", "tidyverse", "dtplyr", "parallel", 
               "caret", "e1071", "viridis")
pacman::p_load_gh("mwesthues/sspredr")


# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
mrna_inbreds <- genos %>%
  filter(Pool != "Hybrid", Data_Type == "mrna") %>%
  split(.$Pool) %>%
  map("G")

# mRNA
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna_df <- mrna %>%
  as_data_frame() %>%
  mutate(G = rownames(mrna)) %>%
  mutate(Group = ifelse(G %in% mrna_inbreds$Dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group))

## Data pre-processing
# Pre-process the data in the following order:
# 1.    Near-zero-variance filtering.
# 2.    Box-Cox tranformation
# 3.    Centering
# 4.    Scaling
mrna_lst <- mrna_df %>% 
  split(.$Group) %>%
  map(as.data.frame) %>%
  map(~column_to_rownames(., var = "G")) %>%
  map(.f = function(x) {
    x$Group <- NULL
    x
  }) %>%
  map(~as.matrix(.))

# Data transformation
trans_lst <- mrna_lst %>%
  map(preProcess,
      method = c("BoxCox", "center", "scale")) %>%
  map2(.y = mrna_lst, .f = predict)

mrna_cor_idx <- trans_lst %>%
  map(cor) %>%
  map(~findCorrelation(x = ., cutoff = 0.999))

trans_lst <- lapply(seq_along(trans_lst), FUN = function(i) {
  mat <- trans_lst[[i]]
  high_cor_idx <- mrna_cor_idx[[i]]
  if (isTRUE(length(high_cor_idx) != 0)) {
    mat <- mat[, -high_cor_idx]
  }
  mat
})


pred_sub_lst <- readRDS("./data/derived/pred_sub_list.RDS")
snp77 <- pred_sub_lst %>% 
  transpose() %>%
  .[names(.) == "ped100_snp77_none"] %>%
  at_depth(.depth = 2, .f = ~.[c(2, 3)]) %>%
  .[[1]]
bare_snp77 <- snp77 %>%
  map(.f = 1)
snp77_cor <- bare_snp77 %>%
  map(cor)
snp77_high_ld <- snp77_cor %>%
  map(~findCorrelation(x = ., cutoff = 0.999))
thinned_snp77 <- lapply(seq_along(bare_snp77), FUN = function(i) {
  mat <- bare_snp77[[i]]
  high_cor_idx <- snp77_high_ld[[i]]
  if (isTRUE(length(high_cor_idx) != 0)) {
    mat <- mat[, -high_cor_idx]
  }
  mat
})
names(thinned_snp77) <- c("Dent", "Flint")
thinned_snp77 <- list(
  list(thinned_snp77[[1]], geno = snp77[[1]][["geno"]]),
  list(thinned_snp77[[2]], geno = snp77[[2]][["geno"]])
  )


snp77_mrna42 <- lapply(seq_along(thinned_snp77), FUN = function(i) {
  snp_grp <- thinned_snp77[[i]] 
  names(snp_grp)[1] <- "snp"
  mrna_grp <- trans_lst[[i]]
  snp_grp$mrna <- mrna_grp
  snp_grp
})
names(snp77_mrna42) <- c("Dent", "Flint")

ped100 <- pred_sub_lst %>% 
  transpose() %>%
  .[names(.) == "ped100_snp77_none"] %>%
  at_depth(.depth = 2, .f = 1) %>%
  .[[1]]

ped100_snp77_mrna42 <- snp77_mrna42
ped100_snp77_mrna42[["Dent"]][["ped"]] <- ped100[[1]]
ped100_snp77_mrna42[["Flint"]][["ped"]] <- ped100[[2]]

ped100_snp77_mrna42_eta <- ped100_snp77_mrna42 %>%
  map(., ~impute2(ped = .$ped,
                  snp = .$snp,
                  mrna = .$mrna,
                  geno = .$geno,
                  as_kernel = TRUE,
                  bglr_model = "BRR"
    )
  )

ped100_snp77_mrna42_eta %>%
  map(2) %>%
  map("X") %>%
  map(~.[match(colnames(.), rownames(.)), ]) %>%
  map(~.[lower.tri(., diag = FALSE)]) %>%
  stack() %>%
  rename(Values = values,
         Group = ind) %>%
  ggplot(aes(x = Values, fill = Group)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#669933", "#FFCC66"))


# Select hybrids for whose parent lines at least one predictor has records.
hybrid <- paste0("DF_", 
                 snp77[[1]][["geno"]],
                 "_",
                 snp77[[2]][["geno"]])

# Agronomic data
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  as_data_frame() %>%
  filter(G %in% hybrid) %>%
  dplyr::select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  dplyr::select(-ADL) %>%
  as.data.frame %>%
  remove_rownames() %>%
  column_to_rownames(var = "G") %>%
  as.matrix

eta <- ped100_snp77_mrna42_eta %>%
  flatten() %>%
  map(function(x) {
    rownames(x$X) <- hybrid
    x
  })

param_df <- expand.grid(Trait = "GTS",
                        Iter = 30000,
                        Run = seq_len(nrow(pheno)))
param_df$Trait <- as.character(param_df$Trait)
set.seed(34923)
param_df <- param_df[sample(rownames(param_df), size = 60), ]
use_cores <- 3L

keep_objs <- c("use_cores", "param_df", "eta", "pheno")
rm(list = ls()[!ls() %in% keep_objs])
gc()

# Determine how much time (hh:mm:ss format) has elapsed since script initiation.
get_elapsed_time <- function(start_time, tz = "CEST") {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt, tz = tz), "%H:%M:%S")
}

start_time <- Sys.time()
# Keep track of how long a job is running.
yhat_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {
  run <- param_df[i, "Run"]
  trait <- param_df[i, "Trait"]
  iter <- as.integer(param_df[i, "Iter"])
  pred <- run_loocv(Pheno = pheno,
                    ETA = eta,
                    hybrid = TRUE,
                    mother_idx = 2,
                    father_idx = 3, 
                    split_char = "_",
                    trait = trait,
                    iter = iter,
                    speed_tst = FALSE,
                    run = run,
                    verbose = FALSE,
                    out_loc = "./tmp/")
  cbind(pred, Iter = iter)
}, mc.cores = use_cores)
res <- rbindlist(yhat_lst)
elapsed_time <- get_elapsed_time(start_time)

res <- res %>%
  rename(Trait = Phenotype,
         Dent = Mother,
         Flint = Father) %>%
  mutate(CV = "LOOCV",
         Elapsed_Time = elapsed_time,
         Date = as.character(Sys.time()),
         Cores = use_cores) %>%
  as.data.table()
res
saveRDS(res, file = "~/playground/new_gtm_snp77_mrna42_none.RDS")
