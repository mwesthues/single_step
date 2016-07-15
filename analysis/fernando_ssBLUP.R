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
pacman::p_load("methods", "BGLR", "Matrix", "parallel", "data.table", "plyr",
               "reshape2", "matrixStats", "testthat")

# Load functions (including Zhang, RadenI and RadenII as part of 'gmat'.
source("./analysis/snp_functions.R")

# COLLECT ARGUMENTS -------------------------------------------------------
if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "2")
  Sys.setenv("CV_RUNS" = "2")
  Sys.setenv("CV_SCHEME" = "custom")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("ITER" = "50000")
  Sys.setenv("MODEL" = "BRR")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
}

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
cv_runs <- as.character(Sys.getenv("CV_RUNS"))
user_cv_scheme <- as.character(Sys.getenv("CV_SCHEME"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.numeric(Sys.getenv("PRIOR_PI_COUNT"))


## CV scheme
# Load the list with test parameters, which was used to estimate the
# computation time.
if (length(user_cv_scheme) != 1 || any(grepl(",", x = user_cv_scheme))) {
  stop("Multiple CV-schemes have to be specified via 'cust_args'!")
}
if (user_cv_scheme == "custom") {
  cust_cv <- readRDS("./data/input/cust_cv.RDS")
  cv_scheme <- cust_cv
  expect_gte(length(cv_scheme), expected = 1)
} else if (user_cv_scheme != "custom") {
  cv_scheme <- user_cv_scheme
} 
cat(cv_scheme, sep = "\n")



## Number of cross-validation runs
# In order to speed up computations, the entire process can be split up into
# multiple subprocesses. Hereto, the possibility of specifying a range for the
# cross-validation runs was enabled. The range must consist of two numbers that
# are separated through "-".
if (isTRUE(unlist(lapply(gregexpr("-", cv_runs), "[[", 1)) != -1)) {
  cv_start <- as.integer(unlist(strsplit(cv_runs, split = "-"))[1])
  cv_end <- as.integer(unlist(strsplit(cv_runs, split = "-"))[2])
  use_runs <- seq(from = cv_start, to = cv_end)
} else {
  use_runs <- seq_len(cv_runs)
}

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


cat(paste0("CV.Runs=", cv_runs, "\n"),
    paste0("Trait=", init_traits, "\n"),
    paste0("Iter=", init_iter, "\n"),
    paste0("CV.Scheme=", cv_scheme, "\n"),
    paste0("Model=", hypred_model, "\n"),
    paste0("VCOV=", g_method, "\n"),
    paste0("Pi=", Pi, "\n"),
    paste0("PriorPiCount=", PriorPiCount, "\n"))


########################
# Record the start time of the process.
start_time <- Sys.time()

# Load the cross-validation scheme
cv_lst <- lapply(seq_along(cv_scheme), FUN = function(i) {
  readRDS(paste0("./data/processed/", cv_scheme[i]))
})
names(cv_lst) <- cv_scheme

# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = init_iter, 
                        CV_Scheme = cv_scheme)
VerboseModel <- FALSE

for (i in 1:ncol(param_df)) param_df[, i] <- as.character(param_df[, i])
numParam <- nrow(param_df)



## - PHENOTYPIC DATA PREPARATION -------------------------------------------
pheno_fls <- list.files("./data/processed/", 
                        pattern = "pheno_BLUE")
y_lst <- lapply(seq_along(pheno_fls), FUN = function(i) {
  Pheno <- read.table(paste0("./data/processed/", pheno_fls[i]),
                      header = TRUE)
  trait_pos <- regexpr("(?<=stage2_)[A-Z]+(?=.txt)", text = pheno_fls[i], 
                       perl = TRUE)
  trait_nm <- substring(pheno_fls[i], first = trait_pos,
                        last = trait_pos + attr(trait_pos, "match.length") - 1)
  Pheno$trait <- trait_nm 
  Pheno <- Pheno[Pheno$check == 0 & Pheno$dent.GTP != 0 &
                 Pheno$flint.GTP != 0, ]
  Pheno <- droplevels(Pheno)
  Pheno <- Pheno[grep("DF_", as.character(Pheno$G)), ]
  Pheno$G <- as.character(Pheno$G)
  Pheno$dent.GTP <- as.character(Pheno$dent.GTP)
  Pheno$flint.GTP <- as.character(Pheno$flint.GTP)
  Pheno
})
Pheno <- rbindlist(y_lst)



## - PREDICTOR DATA PREPARATION ---------------------------------------------
snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna <- mrna[rownames(mrna) %in% rownames(snp), ]
endo_lst <- list(A = snp, M = mrna)

comhybrid <- Pheno[grepl(paste0(rownames(snp), collapse = "|"), 
                         x = Pheno[, G, ]), unique(G), ]
ptlhybrid <- Pheno[dent.GTP %in% rownames(snp) & flint.GTP %in% 
                   rownames(snp), ][, G, ]
comhybrid <- intersect(comhybrid, ptlhybrid)
Pheno <- droplevels(Pheno[Pheno$G %in% comhybrid, , ])
# Store agronomic phenotypes in a matrix.
Y <- dcast.data.table(Pheno, formula = G ~ trait, value.var = "EST")
y_mat <- as.matrix(data.frame(Y[, .(ADF, FETT, GTM, GTS, RPR, STA, XZ), ]))
rownames(y_mat) <- Y[, G, ]
y_mat <- y_mat[complete.cases(y_mat), ]
phenotypes <- colnames(y_mat)
# Extract a vector of parental dent and flint lines, respectively, which have 
# the same length as the number of hybrids. These vectors will be used 
# subsequently for the augmentation of the predictor matrices.
hybrid <- rownames(y_mat)
dent <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 2)
comdent <- unique(dent)
flint <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 3)
comflint <- unique(flint)
snp <- snp[rownames(snp) %in% c(comdent, comflint), ]
snp_nms <- rownames(snp)
mrna <- mrna[rownames(mrna) %in% c(comdent, comflint), ]

# Names of genotypes for which transcriptomic records exist
nm2 <- rownames(mrna)
nm2_dent <- nm2[nm2 %in% comdent]
nm2_flint <- nm2[nm2 %in% comflint]
# Names of genotypes for which transcriptomic records are missing.
nm1 <- setdiff(snp_nms, nm2)
nm1_dent <- nm1[nm1 %in% comdent]
nm1_flint <- nm1[nm1 %in% comflint]



# Build a kernel from the genomic data
hetgrps <- c("Dent", "Flint")
ETA <- lapply(hetgrps, FUN = function(hetgrp) {
  if (hetgrp == "Dent") {
    snp <- snp[match(comdent, rownames(snp)), ]
    nm1 <- nm1_dent
    nm2 <- nm2_dent
    Z1 <- sparse.model.matrix(~-1 + factor(dent[dent %in% nm1]),
                              drop.unused.levels = FALSE)
    Z2 <- sparse.model.matrix(~-1 + factor(dent[dent %in% nm2]),
                              drop.unused.levels = FALSE)
    dent_nm1 <- as.factor(as.character(ifelse(dent %in% nm1, yes = 1, no = 0)))
    X <- sparse.model.matrix(~-1 + dent_nm1, drop.unused.levels = FALSE)
  } else if (hetgrp == "Flint") {
    snp <- snp[match(comflint, rownames(snp)), ]
    nm1 <- nm1_flint
    nm2 <- nm2_flint
    Z1 <- sparse.model.matrix(~-1 + factor(flint[flint%in% nm1]),
                              drop.unused.levels = FALSE)
    Z2 <- sparse.model.matrix(~-1 + factor(flint[flint %in% nm2]),
                              drop.unused.levels = FALSE)
    flint_nm1 <- as.factor(as.character(ifelse(flint %in% nm1,
                                               yes = 1, no = 0)))
    X <- sparse.model.matrix(~-1 + flint_nm1, drop.unused.levels = FALSE)

  }
  M2 <- tcrossprod(mrna) / ncol(mrna)
  M2 <- M2[nm2, nm2]
  snp <- snp[, colVars(snp) != 0]
  A <- gmat[[get("g_method")]](snp, lambda = 0.01)
  A11 <- A[nm1, nm1]
  A12 <- A[nm1, nm2]
  A21 <- A[nm2, nm1]
  A22 <- A[nm2, nm2]
  Ainv <- t(chol(A))
  A_up11 <- Ainv[nm1, nm1]
  A_up12 <- Ainv[nm1, nm2]
  # Eq.21
  M1 <- solve(A_up11, -A_up12 %*% M2)
  J2 <- matrix(-1, nrow = ncol(A12), ncol = 1)
  # Eq.22
  J1 <- solve(A_up11, -A_up12 %*% J2)
  # Eq.10
  epsilon <- A11 - A12 %*% solve(A22) %*% A21
  # Eq.20
  W1 <- Z1 %*% M1
  W2 <- Z2 %*% M2
  W <- as.matrix(rbind(W1, W2))
  U <- as.matrix(rbind(Z1 %*% epsilon, 
                       matrix(0, nrow = nrow(Z2), ncol = ncol(Z1))))
  X_prime <- as.matrix(cbind(X, rbind(Z1 %*% J1, Z2 %*% J2)))
  list(X = X_prime, W = W, U = U)
})
ETA <- unlist(ETA, recursive = FALSE)
for (i in seq_along(ETA)) {
  if (names(ETA)[i] == "X") {
    ETA[[i]] <- list(X = ETA[[i]], model = "FIXED")
  } else {
    ETA[[i]] <- list(X = ETA[[i]], model = hypred_model)
  }
}
names(ETA) <- NULL


full_lst <- vector(mode = "list", length = nrow(param_df))
for (i in seq_len(nrow(param_df))) {
  loocv_lst <- mclapply(seq_len(nrow(y_mat)), FUN = function(j) {
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
out_nm <- paste0("Pred=", predictor, "_Trait=", init_traits,
                 "_Model=", hypred_model, "_VCOV=", g_method, "_Iter=", 
                 niter, "_SnpFilter=", snp_filtr, "_Imputed=", imputed,
                 "_Comparison=", comparison)

