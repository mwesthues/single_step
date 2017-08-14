# Goal: Evaluate whether reducing ETA objects to the number of genotypes that
# will be used as TRN- or TST-hybrids compared to setting the complement of
# genotypes as 'NA' (in the vector of phenotypic values) makes a difference in
# predictive ability.
# Knowing this is crucial for determining how many ETA objects need to build,
# in particular for the maize hybrid data set, which might require more than
# 100,000 ETA objects.
pacman::p_load("tidyverse", "BGLR")

data(wheat)
n_iter = 40000
burnin = 20000

rownames(wheat.X) <- paste0("Line_", seq_len(nrow(wheat.X)))
rownames(wheat.Y) <- rownames(wheat.X)

# Declare genotypes that will neither be used as training set nor as test set
# genotypes as missing.
missing <- paste0(
  "Line_",
  c(12, 90, 95, 135, 215, 290, 310, 335, 400, 410, 420, 534, 555)
)

# Predict these genotypes
set.seed(3490)
tst <- sample(setdiff(rownames(wheat.X), missing), size = 80, replace = FALSE)


## -- KEEP MISSINGS IN ETA OBJECT ---------------------------------------------
# In this scenario we will retain all genotypes, even the ones that we do not
# need as training or test set hybrids, in the ETA object.
# Instead, we will set the phenotypic values for these genotypes, together with
# actual test set hybrids, as missing.
X <- scale(wheat.X) / sqrt(ncol(wheat.X))
y <- wheat.Y[, 1]
y_na <- y
y_na[match(c(tst, missing), names(y_na))] <- NA_real_

fm1 <- BGLR::BGLR(
  y = y_na,
  ETA = list(mrk = list(X = X, model = 'BRR')),
  nIter = n_iter,
  burnIn = burnin,
  saveAt = 'brr_',
  verbose = FALSE
)
yhat <- fm1$yHat
fm1_cor <- cor(
  yhat[match(tst, names(yhat))],
  y[match(tst, names(y))]
)


## -- REMOVE MISSINGS IN ETA OBJECT ---------------------------------------------
# In this scenario we will set all genotypes, that are neither elements of the
# training set nor of the test set, as missing.
X_red <- wheat.X[!rownames(wheat.X) %in% missing, ]
X_red <- scale(X_red) / sqrt(ncol(X_red))
y_red <- wheat.Y[!rownames(wheat.Y) %in% missing, 1]
stopifnot(identical(rownames(X_red), names(y_red)))
y_na_red <- y_red
y_na_red[match(tst, names(y_na_red))] <- NA_real_

fm2 <- BGLR::BGLR(
  y = y_na_red,
  ETA = list(mrk = list(X = X_red, model = 'BRR')),
  nIter = n_iter,
  burnIn = burnin,
  saveAt = 'brr_',
  verbose = FALSE
)
yhat_red <- fm2$yHat
fm2_cor <- cor(
  yhat_red[match(tst, names(yhat_red))],
  y_red[match(tst, names(y_red))]
)

all.equal(fm1_cor, fm2_cor)


## -- CONCLUSION --------------------------------------------------------------
# Apparently, the difference between the two procedures is miniscule, so a
# single ETA object will sufficie for all genotypes (except for core-set inbred
# lines).
