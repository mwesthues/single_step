# Goal: Sample cross-validation folds for prediction with subset of genotypes
# for which records of any endophenotype (incl. root metabolites) are
# available.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("data.table", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

genos <- readRDS("./data/processed/common_genotypes.RDS")
hybrids <- gsub("DF_", replacement = "", x = genos$Hybrid)
n_flint <- ceiling(length(genos$Flint$agro) * 0.8)
n_dent <- ceiling(length(genos$Dent$agro) * 0.8)
cv_options <- expand.grid(n_trn = 500, 
                          min_size = c(40, 50, 60, 70))
cv_lst <- mclapply(seq_len(nrow(cv_options)), FUN = function(i) {
  cv <- sample_cv(hybrids, n_father = n_flint, n_mother = n_dent,
                  n_hyb_trn = cv_options[i, "n_trn"], 
                  rounds = 1000, 
                  min_size = cv_options[i, "min_size"],
                  hybrid_split = "_", 
                  progress = FALSE)
  stopifnot(check_cv(cv) == "success")
  cv
}, mc.cores = 4)
cv_lst <- unlist(cv_lst, recursive = FALSE)

# Plot the distributions of each scenario and each set.
lapply(seq_along(cv_lst), FUN = function(i) {
  dat <- cv_lst[[i]]
  t0 <- colSums(dat == "T0")
  t1 <- colSums(dat == "T1")
  t2 <- colSums(dat == "T2")
  pdf(file = paste0("./data/derived/cv/",
                    names(cv_lst)[i], ".pdf"),
      width = 7, height = 5)
  opar <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  plot.new()
  text(x = 0.5, y = 0.5,
       labels = names(cv_lst)[i], cex = 0.8)
  hist(t0, xlim = c(0, max(t0)), main = "T0", xlab = "")
  abline(v = mean(t0), col = "tomato")
  hist(t1, xlim = c(0, max(t1)), main = "T1", xlab = "")
  abline(v = mean(t1), col = "tomato")
  hist(t2, xlim = c(0, max(t2)), main = "T2", xlab = "")
  abline(v = mean(t2), col = "tomato")
  par(opar)
  dev.off()
})

# Convert CV-scenarios into long-format for compatbility with subsequent
# analyses.
long_cv_lst <- vector(mode = "list", length = length(cv_lst))
names(long_cv_lst) <- names(cv_lst)
lapply(seq_along(cv_lst), FUN = function(i) {
  dat <- cv_lst[[i]]
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  colnames(dat) <- seq_len(ncol(dat))
  dat$Sample_ID <- paste0("DF_", rownames(dat))
  m_dat <- melt(dat, 
                id.vars = "Sample_ID",
                variable.name = "Run",
                value.name = "Set")
  m_dat$Run <- as.numeric(as.character(m_dat$Run))
  m_dat
  saveRDS(m_dat, paste0("./data/processed/cv1000_", names(cv_lst[i]), ".RDS"))
})

