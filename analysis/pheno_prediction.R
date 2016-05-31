# Imputation of mRNA values for samples with genotypic information based on
# samples with both, genotypic and transcriptomic information.
# The imputation is based on the BGLR package.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("BGLR", "parallel", "data.table", "matrixStats", 
               "Matrix", "methods")
source("./analysis/snp_functions.R")

if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "2")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("PREDICTOR" = "snp")
  Sys.setenv("ITER" = "50000")
  Sys.setenv("MODEL" = "BRR_Kernel")
  Sys.setenv("SNP_FILTER" = "FALSE")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("COMPARISON" = "mrna")
  Sys.setenv("IMPUTED" = "both")
}

## ---------------------------------------------------------------------------
# BGLR-parameters
use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))

# phenotypic trait, which will be predicted
init_traits <- as.character(Sys.getenv("TRAIT"))

# the type of predictor, options: snp, mrna, roots
predictor <- as.character(Sys.getenv("PREDICTOR"))
poss_data <- c("mrna", "roots", "snp")
if (!all(predictor %in% poss_data)) {
  stop("Your selection is not part of the set of predictors")
}

# number of MCMC iterations in BGLR
niter <- as.integer(Sys.getenv("ITER"))

# algorithms to use (e.g. "BRR_Kernel" for GBLUP, others: BRR, BayesC, etc.)
hypred_model <- as.character(Sys.getenv("MODEL"))

# Shall the full set of quality-checked SNPs be used for the prediction?
snp_filtr <- as.logical(Sys.getenv("SNP_FILTER"))

# kernel method (only for MODEL=BRR_Kernel), options: Zhang, RadenI, RadenII
g_method <- as.character(Sys.getenv("VCOV"))

# if SNPs are used as predictors, to which predictor shall they be compared?
comparison <- as.character(Sys.getenv("COMPARISON"))
if (comparison == "none" && predictor == "snp") {
  stop("Specify an endophenotype to which SNPs shall be compared")
}

# Concatenate the names of predictors and comparisons. This will be useful if
# the predictors are SNPs since we want to reduce all data (endophenotypes as
# well as phenotypes) to the common set of genotypes to allow for a meaningful
# comparison of prediction results.
if (isTRUE(comparison != "none")) {
  pred_comp <- c(predictor, comparison)
} else {
  pred_comp <- predictor
}

# shall only imputed endophenotypes, only observed endophenotypes or the union 
# of both be used in prediction?
imputed <- as.character(Sys.getenv("IMPUTED"))
imp_opts <- data.frame(IMPUTED = c("no", "yes", "both"),
                       Meaning = c("only observed endophenotypes",
                                   "only imputed endophenotypes",
                                   paste("union of observed and imputed", 
                                         "endohenotpyes")))
print(imp_opts)
burnin <- niter / 2
thinning <- round(niter / 1000)
if (grepl("Bayes", x = hypred_model)) {
  hypred_model <- hypred_model
} else bglr_model <- "BRR"



# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Predictor = predictor,
                        Phenotype = init_traits, 
                        Iter = niter, 
                        stringsAsFactors = FALSE)
VerboseModel <- FALSE
numParam <- nrow(param_df)


## --------------------------------------------------------------------------
# Imputed mRNA data
if (isTRUE(any(c(predictor, comparison) %in% "mrna"))) {
    imp_mrna <- readRDS("./data/processed/imputed_mrna.RDS")
    obs_mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
    obs_mrna <- t(as.matrix(obs_mrna))
  if (isTRUE(imputed == "yes")) {
    mrna <- imp_mrna
  } else if (isTRUE(imputed == "no")) {
    mrna <- obs_mrna
  } else if (isTRUE(imputed == "both"))
    mrna <- rbind(imp_mrna, obs_mrna)
}

if (predictor %in% "snp") {
  snp <- readRDS("./data/processed/snp_mat.RDS")
  if (isTRUE(snp_filtr)) {
    equi_snp_nms <- readRDS("./data/processed/equidistant_snps.RDS")
    snp <- snp[, match(equi_snp_nms, colnames(snp))]
  }
}


# Agronomic data
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
dent_ids <- Pheno[, unique(dent.GTP), ]
flint_ids <- Pheno[, unique(flint.GTP), ]


## --------------------------------------------------------------------------
# Collect all existing endophenotypic feature matrices in one list.
endo_lst <- list()
for (i in seq_along(pred_comp)) {
  pred_nm <- pred_comp[i]
  if (exists(pred_nm)) {
    endo_lst <- c(endo_lst, list(get(pred_nm)))
  }
}
names(endo_lst) <- pred_comp

# Collect all existing endophenotypic feature matrices in one list.
# Get the intersect of genotypes for which all selected endophenotypic data are
# available.
comgeno <- Reduce("intersect", lapply(endo_lst, FUN = rownames))
comgeno <- intersect(comgeno, c(Pheno$dent.GTP, Pheno$flint.GTP))
endo_lst[] <- lapply(seq_along(endo_lst), FUN = function(i) {
  dat <- endo_lst[[i]]
  dat[match(comgeno, rownames(dat)), ]
})

comhybrid <- Pheno[grepl(paste0(comgeno, collapse = "|"), 
                         x = Pheno[, G, ]), unique(G), ]
ptlhybrid <- Pheno[dent.GTP %in% comgeno & flint.GTP %in% comgeno, ][, G, ]
comhybrid <- intersect(comhybrid, ptlhybrid)
Pheno <- droplevels(Pheno[Pheno$G %in% comhybrid, , ])

## ---------------------------------------------------------------------------
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

# After the reduction of endophenotypes and phenotypes to their common core of
# mutual genotypes, remove the 'comparison' endophenotypes from the list of
# endophenotypes. Although not specified in the code, this step will only be
# carried out if a 'comparison' was specified.
pred_lst <- endo_lst[Reduce(`&`,
                            lapply(predictor, 
                                   FUN = grepl,
                                   x = names(endo_lst)))]
# Build kernels from the raw features for GBLUP.
kernels <- lapply(seq_along(pred_lst), FUN = function(i) {
  M_inbred <- pred_lst[[i]]
  pred_nm <- names(pred_lst)[i]
  # Select inbred lines for which phenotypic records are present in the
  # material and augment the matrices so that they match the number of
  # phenotypic records for the hybrids.
  M_flint <- M_inbred[match(comflint, rownames(M_inbred)), ]
  M_dent <- M_inbred[match(comdent, rownames(M_inbred)), ]
  # Purge all monomorphic features.
  M_flint <- M_flint[, colVars(M_flint) != 0]
  M_dent <- M_dent[, colVars(M_dent) != 0]
  
  # Set up kernels, if specified in advance. This will considerably speed up
  # subsequent computations and it equivalent to the results from Bayesian
  # Ridge Regression and GBLUP, respectively, if the number of iterations is
  # sufficiently high.
  if (isTRUE(pred_nm %in% c("mrna", "roots", "snp"))) {
    MG_flint <- gmat[[get("g_method")]](M_flint, lambda = 0.01)
    MG_dent <- gmat[[get("g_method")]](M_dent, lambda = 0.01)
    # generate the design matrices from the kernels (this is (MM')^0.5)
    L_flint <- t(chol(MG_flint))
    L_dent <- t(chol(MG_dent))
  }

  # for the design matrices we need factors with all levels of dent
  # and flint inbred lines that appear as parents
  flint_factor <- factor(flint, levels = comflint)
  dent_factor <- factor(dent, levels = comdent)
  
  # mapping hybrids to inbreds
  Z_flint <- sparse.model.matrix(~-1 + flint_factor, 
                                 drop.unused.levels = FALSE)
  colnames(Z_flint) <- gsub("flint_factor", replacement = "", 
                            x = colnames(Z_flint))
  rownames(Z_flint) <- flint
  Z_dent <- sparse.model.matrix(~-1 + dent_factor,
                                drop.unused.levels = FALSE)
  colnames(Z_dent) <- gsub("dent_factor", replacement = "", 
                           x = colnames(Z_dent))
  rownames(Z_dent) <- dent
  ZL_flint <- as(Z_flint %*% L_flint, "matrix")
  ZL_dent <- as(Z_dent %*% L_dent, "matrix")
  list(Flint = ZL_flint, Dent = ZL_dent)
})
names(kernels) <- names(pred_lst)
# Unlist 'kernels' for compliance with the 'BGLR()' function.
kernels <- unlist(kernels, recursive = FALSE)
names(kernels) <- gsub("[.]", replacement = "_", x = names(kernels))
 
# Add hyperparameters to each list.
ETA <- kernels
for (i in seq_along(kernels)) {
  if (isTRUE(grepl("Bayes", hypred_model))) {
    # If the component is a pedigree element, do not use a BayesB or BayesC,
    # respectively, but use Bayesian Ridge Regression instead.
    if (isTRUE(grepl("ped", names(ETA)[i]))) {
      ETA[[i]] <- list(X = ETA[[i]], model = "BRR")
    } else if (Pi == "unspecified" && PriorPiCount == "unspecified") {
      ETA[[i]] <- list(X = ETA[[i]], model = bglr_model)
    } else {
      ETA[[i]] <- list(X = ETA[[i]], model = bglr_model, probIn = Pi,
                       counts = PriorPiCount)
    } 
  }
  if (!grepl("Bayes", hypred_model)) {
    ETA[[i]] <- list(X = ETA[[i]], model = bglr_model)
  }
}
names(ETA) <- NULL 

## ---------------------------------------------------------------------------
### Run BGLR
full_lst <- vector(mode = "list", length = nrow(param_df))
for (i in seq_len(nrow(param_df))) {
  loocv_lst <- mclapply(seq_len(nrow(y_mat)), FUN = function(j) {
    trait <- param_df[i, "Phenotype"]
    kernel_nms <- param_df[i, "Predictor"]
    res <- data.frame(Predictor = NA_character_,
                      Phenotype = NA_character_,
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
                     nIter = niter,
                     burnIn = burnin,
                     thin = thinning,
                     verbose = VerboseModel)
    # store results
    res$Predictor <- kernel_nms
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
saveRDS(DT,
        file = paste0("./data/derived/pheno_prediction-", out_nm, ".RDS"))

