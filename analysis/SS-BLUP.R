if (!require("pacman")) install.packages("pacman")
pacman::p_load("cpgen", "parallel", "data.table")

if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "2")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("ITER" = "50000")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("CV_Scheme" = "Flint=80_Dent=100_TS=450")
  Sys.setenv("CV_Runs" = "800")
}

## ---------------------------------------------------------------------------
# BGLR-parameters
use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))

# number of cross validation runs
use_runs <- as.numeric(Sys.getenv("CV_Runs"))

# phenotypic trait, which will be predicted
init_traits <- as.character(Sys.getenv("TRAIT"))

# number of MCMC iterations in BGLR
n_iter <- as.integer(Sys.getenv("ITER"))

# kernel method (only for MODEL=BRR_Kernel), options: Zhang, RadenI, RadenII
g_method <- as.character(Sys.getenv("VCOV"))

# Specify the CV-scheme
cv_name <- as.character(Sys.getenv("CV_Scheme"))

# Specify how many iterations to discard as burn in.
n_burnin <- n_iter / 2

# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = n_iter, 
                        stringsAsFactors = FALSE)
VerboseModel <- FALSE
numParam <- nrow(param_df)

pred_info <- readRDS("./data/processed/legarra_inverses.RDS")
y_mat <- pred_info$Agronomic
comhybrid <- rownames(y_mat)
dent <- sapply(strsplit(comhybrid, split = "_"), FUN = "[[", 2)
flint <- sapply(strsplit(comhybrid, split = "_"), FUN = "[[", 3)

# Design matrices
Zdent <- pred_info$DesignDent
Zflint <- pred_info$DesignFlint

# Kernels
GinvDent <- pred_info$GinvDent
GinvFlint <- pred_info$GinvFlint


# Load the cross-validation scheme
cv_lst <- lapply(seq_along(cv_name), FUN = function(i) {
  cv <- readRDS(paste0("./data/processed/cv800_", cv_name[i], ".RDS"))
  cv[cv$Set == "VS0", "Set"] <- "T0"
  cv[cv$Set == "VS1", "Set"] <- "T1"
  cv[cv$Set == "VS2", "Set"] <- "T2"
  cv
})
names(cv_lst) <- cv_name

# Combine all factor levels and store the values in a data frame so that every
# possible computation has its own row.
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = n_iter, 
                        Run = use_runs,
                        CV_Scheme = cv_name)
VerboseModel <- FALSE
numParam <- nrow(param_df)


## THE LOOP -----------------------------------------------------------------
# Get IDs of samples that are supposed to be stored. We want to store the 
# predictions of five cross validation runs (arbitrary number!). Thereby, we 
# can check the predictions for each model without storing a huge amount of
# clutter.
scenario_list <- mclapply(seq_len(numParam), FUN = function(Param) {
  cur_trait <- as.character(param_df[Param, "Phenotype"])
  cur_run <- as.integer(param_df[Param, "Run"])
  cur_iter <- as.integer(param_df[Param, "Iter"])
  cur_cv_name <- as.character(param_df[Param, "CV_Scheme"])
  cur_burnin <- cur_iter / 2
  thinning <- round(cur_iter / 1000)
  ###############################
  # 3. Load CV scheme
  cv <- cv_lst[[cur_cv_name]]
  cv_curr <- cv[cv$Run == as.character(cur_run), ]
  # Sort the CV-scheme so that it matches the phenotypic data.
  cv_curr <- cv_curr[match(comhybrid, cv_curr$Sample_ID), ]

  # Get indices of TS genotypes and set other to NA
  Parents.CV <- as.data.frame(matrix((unlist(strsplit(x = cv_curr$Sample_ID,
                                                      split = "_"))),
                                     ncol = 3, byrow = TRUE))
  colnames(Parents.CV) <- c("DF", "Pd", "Pf")
  for (i in 1:ncol(Parents.CV)) {
    Parents.CV[, i] <- as.character(Parents.CV[, i])
  } 
  stopifnot(all(Parents.CV$Pd %in% dent))
  stopifnot(all(Parents.CV$Pf %in% flint))
  cv_curr <- cbind(cv_curr, Parents.CV[, c("Pd", "Pf")])
  cv_curr <- as.data.frame(cv_curr)

  # Get the indices of all test-set (T0, T1, T2) hybrids, define another vector 
  # of thenotypic records and set the values of the latter to 'NA' if they belong
  # to a genotype that is part of the test set 'tst'.
  tst0 <- cv_curr$Set == "T0"
  tst1 <- cv_curr$Set == "T1"
  tst2 <- cv_curr$Set == "T2"
  tst <- as.logical(tst0 + tst1 + tst2)
  y <- y_mat[, match(cur_trait, colnames(y_mat))]
  yNA <- y
  yNA[tst] <- NA
  
  ##################################
  cv_upd <- cur_cv_name
  cv_upd <- sub(".RDS", replacement = "", x = cv_upd)

  out_name <- paste0(Res.Dir, "Pred=", print_predictor, "_Dom=", 
                     dominance, "_CV.Runs=", cv_runs, "_CV.Scheme=",
                     user_cv_scheme, "_Trait=", init_traits, "_Iter=",
                     init_iter, "_Model=", hypred_model,"_Scale_Pred=",
                     scale_predictor, "_SNP_Filter=", snp_filtr, "_Run=",
                     cur_run, "_VCOV=", g_method, "_FEAT_WEIGHTS=",
                     feat_weight, "_Pi=", Pi, "_PriorPiCount=", PriorPiCount,
                     "_H2=", H2_orig, "_")

  # RUN the GBLUP model.
  mod <- clmm(y = yNA,
              Z = list(Zdent, Zflint),
              ginv = list(GinvDent, GinvFlint),
              niter = cur_iter,
              burnin = cur_burnin,
              verbose = VerboseModel)
  
  # Reassign the name of the predictor to the BGLR model output.
  eta_names <- names(krnl_lst)
  names(mod_BGLR$ETA) <- eta_names
  mod_BGLR 
  # ---> mc.preschedule: task assignment -> TRUE=static, FALSE=dynamic
  # For a discussion of this topic, see pp. 348-350 of Norman Matloff's 
  # "The Art of R Programming".
  }, mc.preschedule = TRUE, mc.cores = use_cores)#end paramaters
names(scenario_list) <- paste0("Predictor=", param_df$Predictor,
                               "_Phenotype=", param_df$Phenotype,
                               "_Iter=", param_df$Iter,
                               "_Run=", param_df$Run,
                               "_CV_Scheme=", as.character(param_df$CV_Scheme),
                               "_VCOV=", g_method,
                               "_FEAT_WEIGHTS=", feat_weight,
                               "_H2_THRESHOLD=", H2_orig)

out_name <- paste0("Pred=", print_predictor, "_Dom=", 
                   dominance, "_CV.Runs=", cv_runs, "_CV.Scheme=",
                   user_cv_scheme, "_Trait=", init_traits, "_Iter=",
                   init_iter, "_Model=", hypred_model, "_Scale_Pred=",
                   scale_predictor,
                   "_SNP_Filter=", snp_filtr, "_VCOV=", g_method,
                   "_FEAT_WEIGHTS=", feat_weight, "_Pi=", Pi,
                   "_PriorPiCount=", PriorPiCount,
                   "_H2_THRESHOLD=", H2_orig)


# Record how much time it took for the process to finish.
end_time <- Sys.time() 
total_time <- as.double(end_time - start_time, units = "hours")
time_df <- data.frame(Start.Time = start_time,
                      End.Time = end_time,
                      Elapsed.Hours = total_time,
                      Predictor = print_predictor,
                      Trait = init_traits,
                      CV.Scheme = user_cv_scheme,
                      CV.Runs = cv_runs,
                      Used.Cores = use_cores,
                      Iter = init_iter,
                      Model = hypred_model)

# Store the BGLR model output.
saveRDS(scenario_list,
        paste0("./81_GBLUP/_results/", out_name, ".RDS"))
if (file.exists("./81_GBLUP/_temp/pred_log.txt")) {
  write.table(time_df, file = "./81_GBLUP/_temp/pred_log.txt",
              sep = "\t", append = TRUE, col.names = FALSE)
} else {
  write.table(time_df, file = "./81_GBLUP/_temp/pred_log.txt",
              sep = "\t")
}

