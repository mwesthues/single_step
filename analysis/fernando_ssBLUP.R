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
# Load functions for the creation of BGLR-kernels.
source("./analysis/fernando_ssBLUP_functions.R")

# COLLECT ARGUMENTS -------------------------------------------------------
if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("ITER" = "60000")
  Sys.setenv("MODEL" = "BRR")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  Sys.setenv("IMPUTATION" = "FALSE")
}

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.numeric(Sys.getenv("PRIOR_PI_COUNT"))
imputation <- as.logical(Sys.getenv("IMPUTATION"))


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
    paste0("PriorPiCount=", PriorPiCount, "\n"),
    paste0("Imputation=", imputation, "\n"))


########################
# Record the start time of the process.
start_time <- Sys.time()
print(start_time)

# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = init_iter)
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
if (isTRUE(imputation)) {
  mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
  mrna <- t(mrna)
  mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
  mrna <- mrna[rownames(mrna) %in% rownames(snp), ]
  endo_lst <- list(A = snp, M = mrna)
} else {
  endo_lst <- list(A = snp)
}

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
flint <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 3)
grp_geno <- list(Dent = dent, Flint = flint)


# Build a kernel from the genomic data
hetgrps <- c("Dent", "Flint")
if (isTRUE(imputation)) {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    impute_ETA(snp = snp, mrna = mrna, grp_geno = grp_geno, hetgrp = hetgrp,
               bglr_model = hypred_model)
  })
} else {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    create_snp_ETA(snp = snp, grp_geno = grp_geno, hetgrp = hetgrp, 
                   bglr_model = hypred_model)
  })
}
ETA <- unlist(ETA, recursive = FALSE)



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
                 "_Pi=", Pi, "_PriorPiCount=", PriorPiCount,
                 "_Imputed=", imputation)
saveRDS(DT, file = paste0("./data/derived/ssBLUP/", out_nm, ".RDS"))

# When did the job finish?
print(Sys.time())
