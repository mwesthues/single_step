# Imputation of mRNA values for samples with genotypic information based on
# samples with both, genotypic and transcriptomic information.
# The imputation is based on the BGLR package.

if (!require("pacman")) install.packages("pacman")
pacman::p_load("BGLR", "parallel", "data.table", "matrixStats", 
               "Matrix", "methods")
source("./analysis/snp_functions.R")

# BGLR-parameters
use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
niter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
snp_filtr <- as.logical(Sys.getenv("SNP_FILTER"))
g_method <- as.character(Sys.getenv("VCOV"))
het_grp <- as.character(Sys.getenv("POOL"))
burnin <- niter / 2
thinning <- round(niter / 1000)
if (grepl("Bayes", x = hypred_model)) {
  hypred_model <- hypred_model
} else bglr_model <- "BRR"


if (isTRUE(interactive())) {
  use_cores <- 3
  niter <- 5e+04
  snp_filtr <- FALSE
  het_grp <- "Dent"
  burnin <- niter / 2
  thinning <- 5
  hypred_model <- "BRR_Kernel"
  bglr_model <- "BRR"
  g_method <- "RadenII"
}


## --------------------------------------------------------------------------
# mRNA data
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(as.matrix(mrna))
 
# SNP data
snp <- readRDS("./data/processed/snp_mat.RDS")
if (isTRUE(snp_filtr)) {
  equi_snp <- readRDS("./data/processed/equidistant_snps.RDS")
  snp <- snp[, match(equi_snp, colnames(snp))]
}
endo_lst <- list(mrna = mrna, snp = snp)

# Collect all existing endophenotypic feature matrices in one list.
# Get the intersect of genotypes for which all selected endophenotypic data are
# available.
comgeno <- Reduce("intersect", lapply(endo_lst, FUN = rownames))
mrna <- mrna[match(comgeno, rownames(mrna)), ]

# For each genotyped inbred line without mRNA data, add a new row to the mRNA
# matrix and fill it with missing values for each mRNA.
# The missing values are supposed to be imputed based on BGLR later on.
mrna_na_nms <- rownames(snp)[!rownames(snp) %in% rownames(mrna)]
mrna_na_mat <- matrix(NA_real_, 
                      nrow = length(mrna_na_nms),
                      ncol = ncol(mrna),
                      dimnames = list(mrna_na_nms,
                                      colnames(mrna)))
mrna_aug <- rbind(mrna, mrna_na_mat)
mrna_aug <- mrna_aug[match(rownames(snp), rownames(mrna_aug)), ]


# Load GTP_IDs of Dent and Flint lines, respectively.
# Select Dent and Flint lines for which both, mRNA and SNP data, are available.
if (het_grp == "Dent") {
  dent <- readRDS("./data/processed/dent_nms.RDS")
  comdent <- dent[dent %in% rownames(snp)]
  inbred <- comdent 
} else if (het_grp == "Flint") {
  flint <- readRDS("./data/processed/flint_nms.RDS")
  comflint <- flint[flint %in% rownames(snp)]
  inbred <- comflint
}
mrna_aug <- mrna_aug[match(inbred, rownames(mrna_aug)), ]
# Build kernels from the raw features for GBLUP.
M_inbred <- snp[match(inbred, rownames(snp)), ]
stopifnot(identical(rownames(mrna_aug), rownames(M_inbred)))
M_inbred <- M_inbred[, colVars(M_inbred) != 0]
MG_inbred <- gmat[[g_method]](M_inbred, lambda = 0.01)
ZL_inbred <- t(chol(MG_inbred))

# Set-up parameters for BGLR
if (isTRUE(grepl("Bayes", hypred_model))) {
  if (Pi == "unspecified" && PriorPiCount == "unspecified") {
    ETA <- list(X = ZL_inbred, model = bglr_model)
  } else {
    ETA <- list(X = ZL_inbred, model = bglr_model, probIn = Pi,
                counts = PriorPiCount)
  } 
} else if (!grepl("Bayes", hypred_model)) {
  ETA <- list(X = ZL_inbred, model = bglr_model)
}
ETA <- list(ETA)


## ---------------------------------------------------------------------------
### Run BGLR
loocv_lst <- mclapply(seq_len(ncol(mrna_aug)), FUN = function(i) {
  mrna_nm <- colnames(mrna_aug)[i]
  res <- list(mRNA = NA_character_,
              y = NA_complex_,
              yhat = NA_complex_)
	# make the training set. Exclude any hybrid that has either parent of the 
  # test sample
	# run the model (GBLUP)
  mod_BGLR <- BGLR(y = mrna_aug[, mrna_nm],
                   ETA = ETA,
                   nIter = niter,
                   burnIn = burnin,
                   thin = thinning,
                   verbose = FALSE)
              
	# store results
  res$mRNA <- mrna_nm
  res$y <- mrna_aug[, mrna_nm]
  res$yhat <- mod_BGLR$yHat
  res
}, mc.cores = use_cores)

out_nm <- paste0("Pool=", het_grp, "_Model=", hypred_model, "_VCOV=", g_method,
                 "_Iter=", niter, "_SnpFilter=", snp_filtr)
saveRDS(loocv_lst,
        file = paste0("./data/derived/mrna_imputation-", out_nm, ".RDS"))

