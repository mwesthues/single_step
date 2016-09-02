# Single-step prediction of hybrid performance based on Fernando et al.: A
# class of Bayesian methods to combine large numbers of genotyped and
# non-genotyped animals for whole-genome analyses. Genetics Selection Evolution
# 2014 46:50
#
# The approach will be adjusted in a sense that we do have complete genomic
# instead of complete pedigree data (i.e. A -> G) and that we have incomplete
# transcriptomic instead of incomplete genomic data (i.e. G -> T).
#
#
#

if (!require("pacman")) install.packages("pacman")
pacman::p_load("methods", "BGLR", "Matrix", "parallel", "data.table", "dplyr",
               "tidyr", "matrixStats", "testthat", "dtplyr")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load_gh("mwesthues/sspredr")
# Load functions (including Zhang, RadenI and RadenII as part of 'gmat'.
# COLLECT ARGUMENTS -------------------------------------------------------
if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("ITER" = "60000")
  Sys.setenv("MODEL" = "BRR")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
}

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.numeric(Sys.getenv("PRIOR_PI_COUNT"))


## Trait selection
poss_traits <- c("GTM", "GTS", "ADL", "FETT", "RFA", "RPR", "STA", "XZ", "ADF")
if (!all(init_traits %in% poss_traits)) {
  stop("The selected trait does not exist")
}

## Hybrid prediction model
if (!hypred_model %in% c("BRR", "BRR_Kernel", "BayesB", "BayesC")) {
  stop("Choose 'BRR', 'BRR_Kernel', 'BayesB', or 'BayesC' as your model")
}

# Specify the method with which to calculate the variance covariance matrix of
# the features.
if (!g_method %in% c("RadenI", "RadenII", "Zhang", "none")) {
  stop("VCOV method must be either 'RadenI', 'RadenII', 'none' or 'Zhang'")
}


cat(paste0("Trait=", init_traits, "\n"),
    paste0("Iter=", init_iter, "\n"),
    paste0("Model=", hypred_model, "\n"),
    paste0("VCOV=", g_method, "\n"),
    paste0("Pi=", Pi, "\n"),
    paste0("PriorPiCount=", PriorPiCount, "\n"))


########################
# Record the start time of the process.
start_time <- Sys.time()

# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = init_iter)
VerboseModel <- FALSE

for (i in 1:ncol(param_df)) param_df[, i] <- as.character(param_df[, i])
numParam <- nrow(param_df)



## - PHENOTYPIC DATA PREPARATION -------------------------------------------
common_geno <- readRDS("./data/processed/common_genotypes.RDS")
Pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
Pheno <- Pheno[G %in% common_geno$Hybrid, , ]



## - PREDICTOR DATA PREPARATION ---------------------------------------------
snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]

y_mat <- Pheno %>%
  select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  select(-ADL) %>%
  as.data.frame %>%
  tibble::column_to_rownames(var = "G") %>%
  as.matrix
hybrid <- rownames(y_mat)
phenotypes <- colnames(y_mat)
dent <- sapply(strsplit(hybrid, split = "_"), FUN = "[[", 2)
comdent <- unique(dent)
flint <- sapply(strsplit(hybrid, split = "_"), FUN = "[[", 3)
comflint <- unique(flint)
snp <- snp[rownames(snp) %in% c(comdent, comflint), ]
snp_nms <- rownames(snp)
mrna <- mrna[rownames(mrna) %in% c(comdent, comflint), ]

# Names of genotypes for which transcriptomic records exist
nm2 <- rownames(mrna)
# Names of genotypes for which transcriptomic records are missing.
nm1 <- setdiff(snp_nms, nm2)


# Build a kernel from the genomic data
hetgrps <- c("Dent", "Flint")
ETA <- lapply(hetgrps, FUN = function(hetgrp) {
  if (hetgrp == "Dent") {
    geno <- dent
  } else if (hetgrp == "Flint") {
    geno <- flint
  }
  # Unique genotypes in current group without transcriptomic records.
  grp_nm1 <- nm1[nm1 %in% geno]
  # Unique genotypes in current group with transcriptomic records.
  grp_nm2 <- nm2[nm2 %in% geno]
  # Hybrid parents in current group without transcriptomic records.
  geno1 <- geno[geno %in% nm1]
  # Hybrid parents in current group with transcriptomic records.
  geno2 <- geno[geno %in% nm2]
  snp <- snp[match(unique(geno), rownames(snp)), ]

  # Design matrix mapping GCA effects from genotypes without transriptomic
  # records to y1.
  Z1 <- sparse.model.matrix(~-1 + factor(geno1),
                            drop.unused.levels = FALSE)
  colnames(Z1) <- gsub("factor\\(geno1\\)", replacement = "", x = colnames(Z1))
  Z1 <- Z1[, match(grp_nm1, colnames(Z1))]
  expect_identical(colnames(Z1), expected = grp_nm1)
  rownames(Z1) <- geno1
  # Design matrix mapping GCA effects from genotypes with transcriptomic
  # records to y2.
  Z2 <- sparse.model.matrix(~-1 + factor(geno2),
                            drop.unused.levels = FALSE)
  colnames(Z2) <- gsub("factor\\(geno2\\)", replacement = "", x = colnames(Z2))
  Z2 <- Z2[, match(grp_nm2, colnames(Z2))]
  expect_identical(colnames(Z2), expected = grp_nm2)
  rownames(Z2) <- geno2

  # Design matrix mapping the fixed effect ("has transcriptomic records or
  # not") to y.
  geno_fct <- as.factor(as.character(ifelse(geno %in% geno1, yes = 1, no = 0)))
  X <- sparse.model.matrix(~-1 + geno_fct, drop.unused.levels = FALSE)
  rownames(X) <- geno

  mrna <- mrna[grp_nm2, ]
  M2 <- mrna[, colVars(mrna) != 0]
  snp <- snp[, colVars(snp) != 0]
  A <- build_kernel(M = snp, lambda = 0.01, algorithm = "RadenII")
  A11 <- A[grp_nm1, grp_nm1]
  A12 <- A[grp_nm1, grp_nm2]
  A21 <- A[grp_nm2, grp_nm1]
  A22 <- A[grp_nm2, grp_nm2]
  Ainv <- solve(A)
  dimnames(Ainv) <- dimnames(A)
  A_up11 <- Ainv[grp_nm1, grp_nm1]
  A_up12 <- Ainv[grp_nm1, grp_nm2]
  # Eq.21
  M1 <- A12 %*% solve(A22) %*% M2
  expect_identical(rownames(M1), grp_nm1)
  J2 <- matrix(-1, nrow = ncol(A12), ncol = 1)
  # Eq.22
  J1 <- A12 %*% solve(A22) %*% J2
  expect_identical(rownames(J1), grp_nm1)
  # Eq.10
  epsilon <- t(chol(solve(Ainv[grp_nm1, grp_nm1])))
  expect_identical(rownames(epsilon), grp_nm1)
  expect_identical(colnames(epsilon), grp_nm1)
  # Eq.20
  W1 <- Z1 %*% M1
  W2 <- Z2 %*% M2
  W <- as.matrix(rbind(W1, W2))
  W <- W[match(geno, rownames(W)), ]
  expect_identical(rownames(W), geno)
  U1 <- Z1 %*% epsilon
  U2 <- Z2 %*% matrix(0, nrow = ncol(Z2), ncol = ncol(U1))
  expect_true(all(U2 == 0))
  U <- as.matrix(rbind(U1, U2)) 
  U <- U[match(geno, rownames(U)), ]
  expect_identical(rownames(U), geno)
  X_prime <- as.matrix(cbind(X, rbind(Z1 %*% J1, Z2 %*% J2)))
  expect_identical(rownames(X_prime), geno)
  param_lst <- list(X = X_prime, W = W, U = U)

  # For one last time, check whether all matrices are sorted according to the
  # order of the agronomic data.
  if (hetgrp == "Dent") {
    y_nm <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 2)
  } else if (hetgrp == "Flint") {
    y_nm <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 3)
  }
  stopifnot(all(unlist(lapply(param_lst, FUN = function(XWU) {
    identical(rownames(XWU), y_nm)
  }))))
  param_lst
})
# Add the model-type to each ETA-list element.
ETA <- unlist(ETA, recursive = FALSE)
for (i in seq_along(ETA)) {
  if (names(ETA)[i] == "X") {
    ETA[[i]] <- list(X = ETA[[i]], model = "FIXED")
  } else {
    ETA[[i]] <- list(X = ETA[[i]], model = hypred_model)
  }
}
names(ETA) <- NULL


# Run LOOCV.
full_lst <- vector(mode = "list", length = nrow(param_df))
y_length <- nrow(y_mat)
for (i in seq_len(nrow(param_df))) {
  loocv_lst <- mclapply(seq_len(y_length), FUN = function(j) {
    trait <- param_df[i, "Phenotype"]
    res <- data.frame(Phenotype = NA_character_,
                      Hybrid = NA_character_,
                      Dent = NA_character_,
                      Flint = NA_character_,
                      y = NA_complex_,
                      yhat = NA_complex_)

	  # make the training set. Exclude any hybrid that has either parent of the 
    # test sample
    ctrain <- intersect(grep(dent[j], x = hybrid, invert = TRUE),
                        grep(flint[j], x = hybrid, invert = TRUE))
	  y <- y_mat[, trait]
	  y[-ctrain] <- NA

	  # run the model (GBLUP)
    mod_BGLR <- BGLR(y = y, 
                     ETA = ETA,
                     nIter = init_iter,
                     burnIn = init_iter / 2,
                     verbose = VerboseModel)
    # store results
    res$Phenotype <- trait
    res$Hybrid <- hybrid[j]
    res$Dent <- dent[j]
    res$Flint <- flint[j]
    res$y <- y_mat[j, trait]
    res$yhat <- mod_BGLR$yHat[j]
    res
	}, mc.cores = use_cores)
  full_lst[[i]] <- rbindlist(loocv_lst)
  full_lst
}
DT <- rbindlist(full_lst)
out_nm <- paste0("Trait=", init_traits, "_Iter=", init_iter,
                 "_Model=", hypred_model, "_VCOV=", g_method, 
                 "_Pi=", Pi, "_PriorPiCount=", PriorPiCount)
saveRDS(DT, file = paste0("./data/derived/ssBLUP/", out_nm, ".RDS"))
