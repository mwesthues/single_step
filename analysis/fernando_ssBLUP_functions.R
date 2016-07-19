impute_ETA <- function(snp, mrna, grp_geno, hetgrp, bglr_model) {
  ## Purpose: Combine data from genotypes with and without transcriptomic data
  ## for single-step hybrid prediction.
  #
  ## Input arguments
  # snp: matrix with SNP genotypes. Genotype-names in rows, marker names in
  # columns.
  #
  # mrna: matrix with mRNA-BLUEs. Genotype names in rows, oligo names in
  # columns.
  #
  # grp_geno: list with genotype names for Dent and Flint corresponding to 
  # hybrid parents in the matrix of agronomic data.
  #
  # hetgrp: vector with string "Dent" or "Flint"
  #
  # bglr_model: algorithm to be used by BGLR

  # Input tests
  test_that("input is correct", {
    expect_is(snp, class = "matrix")
    expect_is(mrna, class = "matrix")
    expect_is(grp_geno, class = "list")
    expect_named(grp_geno, expected = c("Dent", "Flint"))
    expect_is(hetgrp, class = "character")
    expect_match(hetgrp, regexp = "Dent|Flint")          
    expect_gt(nrow(snp), expected = nrow(mrna))
    expect_lt(nrow(snp), expected = ncol(snp))
    expect_lt(nrow(mrna), expected = ncol(mrna))
})

  # Select genotype names based on the current heterotic group.
  if (hetgrp == "Dent") {
    geno <- grp_geno[["Dent"]]
  } else if (hetgrp == "Flint") {
    geno <- grp_geno[["Flint"]]
  }

  # Reduce input data to the available genotypes.
  comgeno <- unique(geno)
  snp <- snp[rownames(snp) %in% comgeno, ]
  mrna <- mrna[rownames(mrna) %in% comgeno, ]
  # Names of genotypes for which transcriptomic records exist
  nm2 <- rownames(mrna)
  # Names of genotypes for which transcriptomic records are missing.
  nm1 <- setdiff(rownames(snp), nm2)

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
  rownames(Z1) <- geno1
  # Design matrix mapping GCA effects from genotypes with transcriptomic
  # records to y2.
  Z2 <- sparse.model.matrix(~-1 + factor(geno2),
                            drop.unused.levels = FALSE)
  colnames(Z2) <- gsub("factor\\(geno2\\)", replacement = "", x = colnames(Z2))
  Z2 <- Z2[, match(grp_nm2, colnames(Z2))]
  rownames(Z2) <- geno2

  # Design matrix mapping the fixed effect ("has transcriptomic records or
  # not") to y.
  geno_fct <- as.factor(as.character(ifelse(geno %in% geno1, yes = 1, no = 0)))
  X <- sparse.model.matrix(~-1 + geno_fct, drop.unused.levels = FALSE)
  rownames(X) <- geno

  mrna <- mrna[grp_nm2, ]
  M2 <- mrna[, colVars(mrna) != 0]
  snp <- snp[, colVars(snp) != 0]
  A <- gmat[[get("g_method")]](snp, lambda = 0.01)
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
  # Eq.10
  epsilon <- t(chol(solve(Ainv[grp_nm1, grp_nm1])))
  # Eq.20
  W1 <- Z1 %*% M1
  W2 <- Z2 %*% M2
  W <- as.matrix(rbind(W1, W2))
  W <- W[match(geno, rownames(W)), ]
  U1 <- Z1 %*% epsilon
  U2 <- Z2 %*% matrix(0, nrow = ncol(Z2), ncol = ncol(U1))
  U <- as.matrix(rbind(U1, U2)) 
  U <- U[match(geno, rownames(U)), ]
  X_prime <- as.matrix(cbind(X, rbind(Z1 %*% J1, Z2 %*% J2)))
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

  # Output tests
  test_that("BGLR-input matrices are sorted and have correct dimensions", {
    expect_identical(colnames(Z1), expected = grp_nm1)
    expect_identical(colnames(Z2), expected = grp_nm2)
    expect_identical(rownames(J1), grp_nm1)
    expect_identical(rownames(epsilon), grp_nm1)
    expect_identical(colnames(epsilon), grp_nm1)
    expect_identical(rownames(W), geno)
    expect_true(all(U2 == 0))
    expect_identical(rownames(U), geno)
    expect_identical(rownames(X_prime), geno)
  })

  list(list(X = X_prime, model = "FIXED"),
       list(X = W, model = bglr_model),
       list(X = U, model = bglr_model))
}


create_snp_ETA <- function(snp, grp_geno, hetgrp, bglr_model) {
  ## Purpose: Set up the SNP-kernel for BGLR. 
  #
  ## Input arguments
  # snp: matrix with SNP genotypes. Genotype-names in rows, marker names in
  # columns.
  #
  # grp_geno: list with genotype names for Dent and Flint corresponding to 
  # hybrid parents in the matrix of agronomic data.
  #
  # hetgrp: vector with string "Dent" or "Flint"
  #
  # bglr_model: the algorithm to be used by BGLR

  # Input tests
  test_that("input is correct", {
    expect_is(snp, class = "matrix")
    expect_is(grp_geno, class = "list")
    expect_named(grp_geno, expected = c("Dent", "Flint"))
    expect_is(hetgrp, class = "character")
    expect_match(hetgrp, regexp = "Dent|Flint")          
    expect_is(bglr_model, class = "character")
    expect_lt(nrow(snp), expected = ncol(snp))
  })

  # Select genotype names based on the current heterotic group.
  if (hetgrp == "Dent") {
    geno <- grp_geno[["Dent"]]
  } else if (hetgrp == "Flint") {
    geno <- grp_geno[["Flint"]]
  }

  comgeno <- unique(geno)
  snp <- snp[rownames(snp) %in% comgeno, ]
  # Sort 
  snp <- snp[match(unique(geno), rownames(snp)), ]
  snp <- snp[, colVars(snp) != 0]

  # Design matrix to map GCA effects to corresponding parents in y.
  Z <- sparse.model.matrix(~-1 + factor(geno), drop.unused.levels = FALSE)
  colnames(Z) <- gsub("factor\\(geno\\)", replacement = "", x = colnames(Z))
  Z <- Z[, match(comgeno, colnames(Z))]
  rownames(Z) <- geno

  # SNP-kernel
  G <- gmat[[get("g_method")]](snp, lambda = 0.01)
  MG <- as(Z %*% G, "matrix")

  # Output tests
  test_that("order of output matrix corresponds to phenotypic data", {
    expect_identical(colnames(MG), expected = comgeno)
    expect_identical(rownames(MG), expected = geno)
  })
  list(list(X = MG, model = bglr_model))
}



run_bglr_loocv <- function(Pheno, ETA, trait, iter, ncores = 1L, 
                           verbose = FALSE) {
  # Goal: leave-one-out cross-validation implemented in BGLR
  #
  ##### Input
  # Pheno: matrix with agronomic traits to be predicted. genotypes in rows and
  # traits in columns.
  #
  # ETA: list with BGLR inpt matrices
  #
  # trait: agronomic trait, which shall be predicted
  #
  # iter: number of MCMC-iterations
  #
  # ncores: number of cores for runs in parallel
  #
  # verbose: logical. shall BGLR-output be printed to stdout?

  # Packages
  require("BGLR")
  require("methods")
  require("parallel")

  # Operations required for input checks.
  nrow_eta <- unique(vapply(ETA, FUN = function(x) {
    nrow(x[["X"]])
  }, FUN.VALUE = integer(1)))

  # Input testing.
  test_that("input is correct", {
    expect_is(y_mat, "matrix")
    expect_gt(nrow(y_mat), ncol(y_mat))
    expect_is(ETA, "list")
    expect_match(rownames(y_mat), regexp = "DF_[0-9]*_[0-9]*")
    expect_length(nrow_eta, 1)
    expect_identical(nrow_eta, expected = nrow(y_mat))
    expect_type(trait, "character")
    expect_length(trait, 1)
    expect_type(iter, "integer")
  })

  # Hybrid progeny and parent genotypes names for LOOCV-matching.
  hybrid <- rownames(Pheno)
  dent <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", 2, 
                 FUN.VALUE = character(1))
  flint <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", 3, 
                  FUN.VALUE = character(1))

  # BGLR implementation
  loocv_lst <- mclapply(seq_len(nrow(Pheno)), FUN = function(run) {
    res <- data.frame(Phenotype = NA_character_,
                      Hybrid = NA_character_,
                      Dent = NA_character_,
                      Flint = NA_character_,
                      y = NA_complex_,
                      yhat = NA_complex_)

	  # set-up the training set. Exclude any hybrid that has either parent of the 
    # test sample
    ctrain <- intersect(grep(dent[run], x = hybrid, invert = TRUE),
                        grep(flint[run], x = hybrid, invert = TRUE))
	  y <- Pheno[, trait]
	  y[-ctrain] <- NA_character_

	  # run the model (GBLUP)
    mod_BGLR <- BGLR(y = y, 
                     ETA = ETA,
                     nIter = iter,
                     burnIn = iter / 2,
                     verbose = verbose)

    # store results
    res$Phenotype <- trait
    res$Hybrid <- hybrid[run]
    res$Dent <- dent[run]
    res$Flint <- flint[run]
    res$y <- y_mat[run, trait]
    res$yhat <- mod_BGLR$yHat[run]

    # Output testing.
    test_that("prediction worked", {
      expect_is(res, "data.frame")
      expect_false(anyNA(res))
      expect_identical(nrow(res), nrow(Pheno))
    })

    # Return results.
    res
  }, mc.cores = ncores)
  DT <- rbindlist(loocv_lst)
  DT[, CV_Scheme := "LOOCV", ]
  DT
}




run_bglr_cv800 <- function(Pheno, ETA, trait, iter, burnin, cv, 
                           verbose = FALSE, run = 1L, out) {
  # Goal: CV800 cross-validation implemented in BGLR
  #
  ##### Input
  # Pheno: matrix with agronomic traits to be predicted. genotypes in rows and
  # traits in columns.
  #
  # ETA: list with BGLR inpt matrices
  #
  # trait: agronomic trait to be predicted
  #
  # iter: number of MCMC iterations
  #
  # burnin: number of MCMC iterations to discard as burn-in
  # 
  # cv: cross-validation scheme order
  #
  # verbose: logical. shall BGLR-output be printed to stdout?
  #
  # run: current CV run
  # 
  # out: location where unwanted BGLR-output will be stored

  # Packages
  require("BGLR")
  require("methods")
  require("parallel")
  require("data.table")

  # Hybrid progeny and parent genotypes names for LOOCV-matching.
  hybrid <- rownames(Pheno)
  dent <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", 2, 
                 FUN.VALUE = character(1))
  flint <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", 3, 
                  FUN.VALUE = character(1))

  ###############################
  # 3. Load CV scheme
  # BGLR implementation
  cv_curr <- cv[cv$Run == as.character(run), ]
  cv_curr <- cv_curr[match(hybrid, cv_curr$Sample_ID), ]
  
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
 
  # Get the indices of all test-set (T0, T1, T2) hybrids, define another vector 
  # of thenotypic records and set the values of the latter to 'NA' if they belong
  # to a genotype that is part of the test set 'tst'.
  tst0 <- cv_curr$Set == "VS0"
  tst1 <- cv_curr$Set == "VS1"
  tst2 <- cv_curr$Set == "VS2"
  tst <- as.logical(tst0 + tst1 + tst2)
  y <- yNA <- Pheno[, trait]
  yNA[tst] <- NA_real_
   
  # RUN BGLR     
  mod_BGLR <- BGLR(y = yNA, 
                   ETA = ETA,
                   nIter = iter,
                   burnIn = burnin,
                   saveAt = out,
                   verbose = verbose)
  pred_ability <- cor(y[tst], mod_BGLR$yHat[tst])
  data.table(Run = run,
             Pred_Ability = pred_ability)
}
