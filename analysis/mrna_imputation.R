# Imputation of mRNA values for samples with genotypic information based on
# samples with both, genotypic and transcriptomic information.
# The imputation is based on the BGLR package.

if (!require("pacman")) install.packages("pacman")
pacman::p_load("BGLR", "parallel", "data.table", "matrixStats", 
               "Matrix", "methods")
source("./analysis/snp_functions.R")

# BGLR-parameters
if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "2")
  Sys.setenv("ITER" = "50000")
  Sys.setenv("REL_SOURCE" = "snp")
  Sys.setenv("PREDICTOR" = "roots")
  Sys.setenv("MODEL" = "BRR_Kernel")
  Sys.setenv("SNP_FILTER" = "FALSE")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("POOL" = "Dent")
}

# Number of cores to use for computations.
use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
# number of MCMC iterations in BGLR
niter <- as.integer(Sys.getenv("ITER"))
# Relationship source (pedigree, snp)
rel_nm <- as.character(Sys.getenv("REL_SOURCE"))
# which predictor shall be imputed?
predictor <- as.character(Sys.getenv("PREDICTOR"))
# algorithms to use (e.g. "BRR_Kernel" for GBLUP, others: BRR, BayesC, etc.)
hypred_model <- as.character(Sys.getenv("MODEL"))
# Shall the full set of quality-checked SNPs be used for the prediction?
snp_filtr <- as.logical(Sys.getenv("SNP_FILTER"))
# kernel method (only for MODEL=BRR_Kernel), options: Zhang, RadenI, RadenII
g_method <- as.character(Sys.getenv("VCOV"))
# for which heterotic group shall the computation be carried out?
het_grp <- as.character(Sys.getenv("POOL"))
burnin <- niter / 2
thinning <- round(niter / 1000)
if (grepl("Bayes", x = hypred_model)) {
  hypred_model <- hypred_model
} else bglr_model <- "BRR"



## --------------------------------------------------------------------------
#  endophenotypic data
if (predictor == "mrna") {
  endo <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
  endo <- t(as.matrix(endo))
} else if (predictor == "roots") {
  endo <- readRDS("./data/processed/mpi_imputed_root_blues.RDS")
  # Rename genotypes in the root data set from their pretty names to their
  # GTP-IDs.
  smp_conv <- data.frame(fread("./data/processed/uhoh_smp_id_P235-update.txt",
                               sep = "\t"))
  geno <- unique(rownames(endo))
  geno <- geno[grep("x|Check", x = geno, invert = TRUE)]
  geno <- geno[!is.na(smp_conv[match(geno, smp_conv$name.pretty), "id.GTP"])]
  endo <- endo[match(geno, rownames(endo)), ]
  gtp_id <- smp_conv[match(geno, smp_conv$name.pretty), "id.GTP"]
  geno_match <- data.frame(Pretty_Name = geno, GTP_ID = gtp_id)
  rownames(endo) <- geno_match[match(rownames(endo), geno_match$Pretty_Name),
                               "GTP_ID"]
}
 
# relationship information source
if (rel_nm == "snp") {
  rel_src <- readRDS("./data/processed/snp_mat.RDS")
  if (isTRUE(snp_filtr)) {
    equi_snp <- readRDS("./data/processed/equidistant_snps.RDS")
    rel_src <- rel_src[, match(equi_snp, colnames(rel_src))]
  }
} else if (rel_nm == "ped") {
  rel_src <- readRDS("./data/processed/ped-datafull-GTP.RDS") 
  rel_src <- as.matrix(rel_src) * 2
  rel_src[is.na(rel_src)] <- 0
}

endo_lst <- list(endo = endo, rel_src = rel_src)

# Collect all existing endophenotypic feature matrices in one list.
# Get the intersect of genotypes for which all selected endophenotypic data are
# available.
comgeno <- Reduce("intersect", lapply(endo_lst, FUN = rownames))
endo <- endo[match(comgeno, rownames(endo)), ]

# For each genotyped inbred line without mRNA data, add a new row to the mRNA
# matrix and fill it with missing values for each mRNA.
# The missing values are supposed to be imputed based on BGLR later on.
endo_na_nms <- rownames(rel_src)[!rownames(rel_src) %in% rownames(endo)]
endo_na_mat <- matrix(NA_real_, 
                      nrow = length(endo_na_nms),
                      ncol = ncol(endo),
                      dimnames = list(endo_na_nms,
                                      colnames(endo)))
endo_aug <- rbind(endo, endo_na_mat)
endo_aug <- endo_aug[match(rownames(rel_src), rownames(endo_aug)), ]


# Load GTP_IDs of Dent and Flint lines, respectively.
# Select Dent and Flint lines for which both, mRNA and SNP data, are available.
if (het_grp == "Dent") {
  dent <- readRDS("./data/processed/dent_nms.RDS")
  comdent <- dent[dent %in% rownames(rel_src)]
  inbred <- comdent 
} else if (het_grp == "Flint") {
  flint <- readRDS("./data/processed/flint_nms.RDS")
  comflint <- flint[flint %in% rownames(rel_src)]
  inbred <- comflint
}
endo_aug <- endo_aug[match(inbred, rownames(endo_aug)), ]
# Build kernels from the raw features for GBLUP.
M_inbred <- rel_src[match(inbred, rownames(rel_src)), ]
stopifnot(identical(rownames(endo_aug), rownames(M_inbred)))
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
loocv_lst <- mclapply(seq_len(ncol(endo_aug)), FUN = function(i) {
  endo_nm <- colnames(endo_aug)[i]
  res <- list(mRNA = NA_character_,
              y = NA_complex_,
              yhat = NA_complex_)
	# make the training set. Exclude any hybrid that has either parent of the 
  # test sample
	# run the model (GBLUP)
  mod_BGLR <- BGLR(y = endo_aug[, endo_nm],
                   ETA = ETA,
                   nIter = niter,
                   burnIn = burnin,
                   thin = thinning,
                   verbose = FALSE)
              
	# store results
  res$mRNA <- endo_nm
  res$y <- endo_aug[, endo_nm]
  res$yhat <- mod_BGLR$yHat
  res
}, mc.cores = use_cores)

out_nm <- paste0("Pool=", het_grp, "_Model=", hypred_model, "_VCOV=", g_method,
                 "_Predictor=", predictor, "_Iter=", niter, 
                 "_SnpFilter=", snp_filtr, "_RelSource=", rel_nm)
saveRDS(loocv_lst,
        file = paste0("./data/derived/endo_imputation-", out_nm, ".RDS"))

