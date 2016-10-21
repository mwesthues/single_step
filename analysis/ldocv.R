# Leave-d-out cross-validation
# Gianola et al. (2016), doi:10.1534/g3.116.033381
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table", "sspredr", "BGLR", "Matrix")


## -- DATA SELECTION -----------------------------------------------------
genos <- readRDS("./data/processed/common_genotypes.RDS")
# hybrids
hybrid <- genos %>%
  filter(Pool == "Hybrid") %>%
  .$G

# snp
snp <- readRDS("./data/processed/snp_mat.RDS")
snp_nms <- genos %>%
  filter(Data_Type == "snp") %>% 
  .$G
snp <- snp %>%
  .[match(snp_nms, rownames(.)), ]


## -- PREDICTOR AND AGRONOMIC DATA PREPARATION ----------------------------
# Agronomic data
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  as_data_frame() %>%
  filter(G %in% hybrid) %>%
  dplyr::select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  dplyr::select(-ADL) %>%
  as.data.frame %>%
  remove_rownames() %>%
  column_to_rownames(var = "G") %>%
  as.matrix


# Extract the names of the dent and flint parents, respectively.
dent_flint_lst <- pheno %>%
  rownames() %>%
  as.list() %>%
  map(~strsplit(., split = "_")) %>%
  at_depth(.depth = 2, .f = function(x) x[c(2, 3)]) %>%
  flatten() %>%
  transpose()
names(dent_flint_lst) <- c("Dent", "Flint")




# Function for the setup of the feature matrix. 
build_ZL <- function(M, Z, use_kernel = TRUE, bglr_model = "BRR") {
  # Input tests
  if (class(M) != "matrix") stop("'M' is not a matrix")
  stopifnot(bglr_model %in% c("FIXED", "BRR", "BayesA", "BayesB", "BayesC",
                              "RKHS"))
  comgeno <- unique(rownames(Z))
  M <- M[rownames(M) %in% comgeno, ]
  # Sort
  M <- M[match(comgeno, rownames(M)), ]
  M <- M[, matrixStats::colVars(M) != 0]

  if (isTRUE(use_kernel)) {
    # scale with the standard deviation of original features
    W <- scale(M, center = TRUE, scale = TRUE)
    # G = WW' / m, where m is the number of features
    G <- tcrossprod(W) / ncol(W)
    diag(G) <- diag(G) + 0.01
    L <- G %>% t() %>% chol()
  } else {
    L <- M
  }
   # x-kernel
  ZL <- as.matrix(Z %*% L)
  # Output tests
  stopifnot(identical(rownames(ZL), rownames(Z)))
  list(X = ZL, model = bglr_model)
}
create_Z_matrix <- function(x, geno = NULL) {
   # Design matrix to map GCA effects to corresponding parents in y.
  comx <- unique(x)
  if (is.null(geno)) {
    geno <- x
  }
  Z <- Matrix::sparse.model.matrix(~-1 + factor(geno),
                                   drop.unused.levels = FALSE)
  rownames(Z) <- geno
  colnames(Z) <- gsub("factor\\(geno\\)", replacement = "", x = colnames(Z))
  Z <- Z[, match(comx, colnames(Z))]
  Z
}

ZL_lst <- snp %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  mutate(Pool = ifelse(G %in% dent_flint_lst$Dent,
                       yes = "Dent", no = "Flint")) %>%
  split(.$Pool) %>%
  map(~remove_rownames(.)) %>%
  map(~column_to_rownames(., var = "G")) %>%
  map(function(x) {
    x$Pool <- NULL
    x
  }) %>%
  map(as.matrix) %>%
  map2(.y = dent_flint_lst %>% map(create_Z_matrix),
       .f = build_ZL, use_kernel = TRUE)


# Application to multiple traits.
bglr_lst <- pheno %>%
  as.data.frame() %>%
  as.list() %>%
  map(~BGLR(y = .,
            ETA = ZL_lst,
            nIter = 15000,
            burnIn = 5000,
            verbose = TRUE,
            saveAt = "./tmp/"))

# Compute lambda
varA_lst <- bglr_lst %>%
  map("ETA") %>%
  at_depth(.depth = 2, .f = function(x) x$varB)
varE_lst <- bglr_lst %>%
  map(function(x) x$varE)
lambda_lst <- varA_lst %>%
  map(transpose) %>%
  flatten() %>% 
  map2(.y = varE_lst,
       .f = function(vA, vE) vE / vA) %>%
  transpose() %>%
  map(function(x) {
    names(x) <- names(varA_lst)
    x
  })


G_lst <- ZL_lst %>%
  map(function(x) x[["X"]]) %>%
  map(function(x) {
    x[match(unique(rownames(x)), rownames(x)), ]
  })

Z_lst <- dent_flint_lst %>%
  map(function(x) {
    Matrix::sparse.model.matrix(~-1 + factor(x), 
                                drop.unused.levels = FALSE)
  }) %>%
  map(function(x) {
    colnames(x) <- gsub("factor\\(x\\)", x = colnames(x), replacement = "")
    x 
  }) %>%
#  map2(.y = dent_flint_lst,
#       .f = function(x, y) {
#    x[, match(y, colnames(x))]
#  }) %>%
  map2(.y = dent_flint_lst,
       .f = function(x, y) {
    rownames(x) <- y
    x
  })

# Eq. 109
C_mat <- lapply(seq_len(ncol(pheno)), FUN = function(i) {

  Cinv <- lambda_lst %>%
    map(i) %>%
    map2(.y = G_lst,
         .f = function(lambda, G) {
      diag(1, nrow = nrow(G), ncol = ncol(G)) + solve(G) * lambda 
    }) %>%
    map(solve)

    map2(.x = Z_lst, .y = Cinv, .f = function(x, y) x %*% y %*% t(x))
})
names(C_mat) <- colnames(pheno)
    

lapply(seq_len(pheno), FUN = function(i) {
  y <- bglr_lst %>%
    map("y") %>%
    .[[i]] %>%
    as.matrix()

  C_mat %>%
    .[[i]] %>%
    map(function(x) x %*% y) %>%
    map(as.matrix) %>%
    map(c) %>%
    map(summary)
})

