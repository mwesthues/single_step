# This script contains a cross-validation function based on the method of
# Technow et al. (2014). "Genome Propeties and Prospects of Genomic Prediction
# of Hybrid Performance in a Breeding Program of Maize."
# Genetics, 197: 1343-1355
if (!require("pacman")) install.packages("pacman")
pacman::p_load("reshape2", "data.table", "parallel")

## -------------------------------------------------------------------------
### FUNCTIONS
sample_cv <- function(x, n_dent, n_flint, ts_size) {
  #-- PURPOSE
  # Sampling function for generating CV1000-scenarios.

  # -- INPUT
  # x: data.frame with three columns: 1) hybrid names, 2) parental Dent names,
  #                                   3) parental Flint names
  # n_dent: vector with number of parental Dent lines to sample.
  # n_flint: vector with number of parental Flint lines to sample.
  # ts_size: scalar (vector) specifying the size of the training set
  # rnd_seed: set a random seed
  stopifnot(class(x) == "data.frame")
  stopifnot(all(sapply(list(n_dent, n_flint, ts_size), FUN = is.vector)))
  if (!(all(colnames(x) == c("Hybrid", "HybDent", "HybFlint")))) {
    stop("Rename data frame columns as 'Hybrid', 'HybDent' and 'HybFlint'")
  }

  # -- COMPUTATION
  # Get the original random seed, which shall be stored upon exit.
  # Sample parental Dent and Flint lines, respectively.
  smp_dent_nms <- sample(unique(x$HybDent), size = n_dent, replace = FALSE)
  smp_flint_nms <- sample(unique(x$HybFlint), size = n_flint, replace = FALSE)

  # Check whether actual parental inbred lines match the sampled parental
  # inbred lines or not.
  x$Dent_Match <-  as.numeric(x$HybDent %in% smp_dent_nms)
  x$Flint_Match <- as.numeric(x$HybFlint %in% smp_flint_nms)
  x$Set_Count <- x$Dent_Match + x$Flint_Match
  # Get indices of all hybrids, which are progeny of two sampled inbred lines.
  ts_candidates <- which(x$Set_Count == 2)

  # Shuffle the hybrid parents
  rnd_x <- x[ts_candidates, ]
  rnd_x <- rnd_x[sample(seq_len(nrow(rnd_x)), replace = FALSE), ]

  # Match hybrids and inbreds.
  dent_hyb <- sapply(seq_along(smp_dent_nms), FUN = function(i) {
    rnd_x[smp_dent_nms[i] == rnd_x$HybDent, "Hybrid"][1]
  })
  flint_hyb <- sapply(seq_along(smp_flint_nms), FUN = function(i) {
    rnd_x[smp_flint_nms[i] == rnd_x$HybFlint, "Hybrid"][1]
  })

  # Remove duplicated, sampled TS-hybrids.
  ts_cond_hyb <- union(dent_hyb, flint_hyb)
  if (all(smp_dent_nms %in% x[x$Hybrid %in% ts_cond_hyb, "HybDent"]) &&
      all(smp_flint_nms %in% x[x$Hybrid %in% ts_cond_hyb, "HybFlint"])) {
    
    # Add the difference between already sampled TS-hybrids and the specified
    # TS-size.
    ts_add_hyb <- sample(rnd_x[!rnd_x$Hybrid %in% ts_cond_hyb, "Hybrid"],
                         size = ts_size - length(ts_cond_hyb),
                         replace = FALSE)
    ts_hyb <- c(ts_cond_hyb, ts_add_hyb)
    stopifnot(isTRUE(length(ts_hyb) == ts_size))
    if (any(duplicated(ts_hyb))) {
      stop("Duplicated TS hybrids were sampled")
    }
  
    # Assign labels to the hybrids.
    x[x$Set_Count == 0, "Set_Label"] <- "T0"
    x[x$Set_Count == 1, "Set_Label"] <- "T1"
    x[x$Set_Count == 2, "Set_Label"] <- "T2"
    x[x$Hybrid %in% ts_hyb, "Set_Label"] <- "TS"
  
  } else {
    x$Set_Label <- "error"
  }
    x$Set_Label
}

detect_small_sets <- function(x, min_set_size) {
  #-- Purpose: 
  # Detect any column in which not all sets (i.e. T0, T1, T2) contain at least
  # 'min_set_size' elements and remove duplicated CV-runs.

  #-- INPUT
  # x: a matrix with set labels for each CV-run. CV-runs are stored in columns.
  # min_set_size: the minimum frequency of each set (T0, T1, T2) per CV-run.
  if (class(x) != "matrix") {
    stop("x has to be of class 'matrix'")
  }
  if (!is.numeric(min_set_size) || length(min_set_size) != 1) {
    stop("min_set_size has to be a numeric scalar")
  }

  #-- COMPUTATIONS
  colSums(x == "T0") >= min_set_size &
  colSums(x == "T1") >= min_set_size &
  colSums(x == "T2") >= min_set_size &
  !duplicated(x, MARGIN = 2)
}


# Unique combinations of data.frames.
expand.grid.df <- function(...) {
  Reduce(function(...) merge(..., by = NULL), list(...))
}
## -------------------------------------------------------------------------
### ANALYSIS
# Load common samples.
pred_info <- readRDS("./data/processed/legarra_inverses.RDS")
dent <- pred_info$NamesDent
flint <- pred_info$NamesFlint
hybrid <- pred_info$NamesHybrid
use_cores <- 4


# Extract the names of the parental inbred lines, which constitute the hybrids.
hd <- sapply(strsplit(hybrid, split = "_"), FUN = "[[", 1)
all(dent %in% hd)
hf <- sapply(strsplit(hybrid, split = "_"), FUN = "[[", 2)
all(flint %in% hf)

# Extract a baseline data frame, which contains the names of the hybrids and
# the names of the corresponding parental Dent and Flint inbred lines,
# respectively.
base_df <- data.frame(Hybrid = hybrid,
                      HybDent = hd,
                      HybFlint = hf,
                      stringsAsFactors = FALSE)


# Explore a plethora of CV-set compositions.
flint_sizes <- data.frame(Flint_Sizes = seq(from = 30, to = 130, by = 10))
dent_sizes <- data.frame(Dent_Sizes = seq(from = 20, to = 100, by = 10))
ts_sizes <- data.frame(TS_Sizes = seq(from = 150, to = 800, by = 50))
tst_cv_size <- 1500 # allows for errors
fin_cv_size <- 800 # final CV-size
min_set_size <- 10
cv_params <- expand.grid.df(flint_sizes, dent_sizes, ts_sizes)
cv_params <- cv_params[cv_params$Flint_Sizes < cv_params$Dent_Sizes, ]
cv_lst <- mclapply(seq_len(nrow(cv_params)), FUN = function(i) {
  scen_lst <- do.call(cbind, lapply(seq_len(tst_cv_size), FUN = function(j) {
  try(sample_cv(x = base_df, 
                n_dent = cv_params$Dent_Sizes[i],
                n_flint = cv_params$Flint_Sizes[i],
                ts_size = cv_params$TS_Sizes[i]),
      silent = TRUE)
  }))
#  scen_lst
#})
#lapply(cv_lst, FUN = function(x) {
#  if (isTRUE(all(c(x) %in% c("T0", "T1", "T2", "TS", "error")))) {
#    cv_mat <- x[, as.logical(colMeans(x != "error"))]
#  } else cv_mat <- NULL
#  if (!is.null(cv_mat) && ncol(cv_mat) >= fin_cv_size) {
#    large_sets <- cv_mat[, detect_small_sets(cv_mat, 
#                                             min_set_size = min_set_size)]
#  }
#  if (exists("large_sets")) {
#    dim(large_sets)
#  } else NULL
#})
  succ_cv <- all(c(scen_lst) %in% c("T0", "T1", "T2", "TS", "error"))
  if (isTRUE(succ_cv)) {
    cv_mat <- scen_lst
    cv_mat <- cv_mat[, as.logical(colMeans(cv_mat != "error"))]
    if (ncol(cv_mat) >= fin_cv_size) {
      large_sets <- cv_mat[, detect_small_sets(cv_mat, 
                                               min_set_size = min_set_size)]
    }
    if (exists("large_sets") && ncol(large_sets) >= fin_cv_size) {
      res <- large_sets[, sample(seq_len(ncol(large_sets)),
                                 size = fin_cv_size, replace = FALSE)]
    } else {
      res <- NULL
    }
    res
  } else {
    res <- NULL
  }
  res
}, mc.cores = use_cores)
names(cv_lst) <- sapply(seq_len(nrow(cv_params)), FUN = function(i) {
  dat <- cv_params[i, ]
  paste0("Flint=", dat$Flint_Sizes, "_", "Dent=", dat$Dent_Sizes, "_",
         "TS=", dat$TS_Sizes)
})


# Remove all combinations, which could not be set-up.
clean_cv_lst <- Filter(function(x) !is.null(x), x = cv_lst)
# Final check whether all CV-runs within each scenario are unique.
if (any(sapply(clean_cv_lst, FUN = function(x) {
          any(duplicated(x, MARGIN = 2))
        }) == TRUE)) {
  warning("Still duplicated CV-runs in some scenario(s)")
}
saveRDS(clean_cv_lst, 
        file = "./data/processed/possible_cv_scenarios.RDS")


# Plot the distributions of each scenario and each set.
lapply(seq_along(clean_cv_lst), FUN = function(i) {
  dat <- clean_cv_lst[[i]]
  t0 <- colSums(dat == "T0")
  t1 <- colSums(dat == "T1")
  t2 <- colSums(dat == "T2")
  png(file = paste0("./data/derived/cv/",
                    names(clean_cv_lst)[i], ".png"),
      width = 1500, height = 800)
  opar <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
  plot.new()
  text(x = 0.5, y = 0.5,
       labels = names(clean_cv_lst)[i], cex = 2)
  hist(t0, xlim = c(0, max(t0)), main = "VS0", xlab = "")
  abline(v = mean(t0), col = "tomato")
  hist(t1, xlim = c(0, max(t1)), main = "VS1", xlab = "")
  abline(v = mean(t1), col = "tomato")
  hist(t2, xlim = c(0, max(t2)), main = "VS2", xlab = "")
  abline(v = mean(t2), col = "tomato")
  par(opar)
  dev.off()
})


# Convert CV-scenarios into long-format for compatbility with subsequent
# analyses.
long_cv_lst <- vector(mode = "list", length = length(clean_cv_lst))
names(long_cv_lst) <- names(clean_cv_lst)
lapply(seq_along(clean_cv_lst), FUN = function(i) {
  dat <- clean_cv_lst[[i]]
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  colnames(dat) <- seq_len(ncol(dat))
  dat$Sample_ID <- paste0("DF_", base_df$Hybrid)
  m_dat <- melt(dat, 
                id.vars = "Sample_ID",
                variable.name = "Run",
                value.name = "Set")
  m_dat$Run <- as.numeric(as.character(m_dat$Run))
  m_dat[m_dat$Set != "TS", "Set"] <- gsub("T", replacement = "VS", 
                                          m_dat[m_dat$Set != "TS", "Set"])
  m_dat
  saveRDS(m_dat, 
          paste0("./data/processed/cv800_", names(clean_cv_lst[i]), ".RDS"))
})

