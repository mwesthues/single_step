# Genotype-Wise mRNA-Imputation
pacman::p_load("sspredr", "tidyverse", "forcats", "purrr", "methods", 
               "parallel")
use_cores <- 3L

# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
dent <- genos %>% 
  filter(Pool == "Dent", Data_Type == "mrna") %>%
  .$G
flint <- genos %>%
  filter(Pool == "Flint", Data_Type == "mrna") %>%
  .$G

# Remove the 'Group' variable after splitting Dent and Flint data.
rm_grp <- function(x) {
  x$Group <- NULL
  x
}

# mRNA
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna_lst <- mrna %>%
  as_data_frame() %>%
  mutate(G = rownames(mrna)) %>%
  mutate(Group = ifelse(G %in% dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(~as.matrix(.))

# SNPs
snp <- readRDS("./data/processed/snp_mat.RDS")
snp_lst <- snp %>%
  as_data_frame() %>%
  mutate(G = rownames(snp)) %>%
  filter(G %in% c(dent, flint)) %>%
  mutate(Group = ifelse(G %in% dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(~as.matrix(.))

# Pedigree
ped <- readRDS("./data/processed/ped-datafull-GTP.RDS")
ped_lst <- ped %>%
  mutate(G = rownames(ped)) %>%
  filter(G %in% c(dent, flint)) %>%
  mutate(Group = ifelse(G %in% dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(as.matrix) %>%
  map(function(x) x[, match(rownames(x), colnames(x))]) %>%
  map(function(x) x * 2) %>%
  map(function(x) {
    x[is.na(x)] <- 0
    x
  })


## Function for the imputation of a single mRNA.
impute_mrna_genotype <- function(x, y) {
  geno <- rownames(x)
  # Names of genotypes for which transcriptomic records exist
  nm2 <- rownames(y)
  # Names of genotypes for which transcriptomic records are missing.
  nm1 <- setdiff(rownames(x), nm2)
  # Hybrid parents not in y.
  geno1 <- geno[geno %in% nm1]
  # Hybrid parents in y.
  geno2 <- geno[geno %in% nm2]
  x <- x[match(unique(geno), rownames(x)), ]
  # Design matrix mapping GCA effects from genotypes with transcriptomic
  # records to y2.
  Z2 <- Matrix::sparse.model.matrix(~-1 + factor(geno2),
                                    drop.unused.levels = FALSE)
  colnames(Z2) <- gsub("factor\\(geno2\\)", replacement = "", x = colnames(Z2))
  Z2 <- Z2[, match(nm2, colnames(Z2))]
  rownames(Z2) <- geno2
  # Hybrid parents not in y.
  geno1 <- geno[geno %in% nm1]
  # Hybrid parents in y.
  geno2 <- geno[geno %in% nm2]
  x <- x[match(unique(geno), rownames(x)), ]
  y <- y[nm2, ]
  M2 <- y[, matrixStats::colVars(y) != 0]
  x <- x[, matrixStats::colVars(x) != 0]
  A <- build_kernel(M = x, lambda = 0.01, algorithm = "RadenII")
  A11 <- A[nm1, nm1, drop = FALSE]
  A12 <- A[nm1, nm2, drop = FALSE]
  A21 <- A[nm2, nm1, drop = FALSE]
  A22 <- A[nm2, nm2, drop = FALSE]
  Ainv <- solve(A)
  dimnames(Ainv) <- dimnames(A)
  A_up11 <- Ainv[nm1, nm1]
  A_up12 <- Ainv[nm1, nm2]
  # Eq.21
  M1 <- A12 %*% solve(A22) %*% M2
  J2 <- matrix(-1, nrow = ncol(A12), ncol = 1)
  # Eq.22
  J1 <- A12 %*% solve(A22) %*% J2
  # Eq.10
  epsilon <- t(chol(solve(Ainv[nm1, nm1])))
  # Eq.20
  W2 <- Z2 %*% M2
  w1 <- A12 %*% solve(A22) %*% Z2 %*% M2
  as.matrix(w1)
}

# Shuffled version of the imputation function above.
shuffled_imputation <- function(x, y) {
  cor_lst <- lapply(seq(from = 1, to = nrow(y) - 2), FUN = function(i) {
    idx <- sample(nrow(y), size = i, replace = FALSE)
    z <- y[-idx, ]
    imp_mrna <- impute_mrna_genotype(x = x, y = z)
    orig_mrna <- y[match(rownames(imp_mrna), rownames(y)), , drop = FALSE]
    stopifnot(identical(rownames(imp_mrna), rownames(orig_mrna)))
    cor(c(imp_mrna), c(orig_mrna))
  })
  cor_lst %>%
    flatten_dbl() %>%
    tbl_df() %>%
    rownames_to_column(var = "Number_NA_Geno") %>%
    rename(Value = value)
}


## Correlations between vectors of original and imputed mRNAs.
### General idea:
#  Impute transcriptomic data for different numbers of genotypes (from 1 to
#  number of genotypes with transcriptomic data).
#  For each number of genotypes without transcriptomic data, randomly sample 
#  100 sets of genotypes, which will be imputed.
#  This will allow us to aggregate and summarize the correlation between 
#  observed and imputed mRNAs with high confidence.
repl <- 1000
donors <- c("snp_lst", "ped_lst")
imp_lst <- lapply(donors, FUN = function(type) {
  par_lst <- vector(mode = "list", length = repl)
  names(par_lst) <- paste0("Replication_", seq_len(repl))
  par_lst[] <- mclapply(par_lst, FUN = function(iter, ...) {
    try(map2(.x = get(type), .y = mrna_lst, .f = shuffled_imputation))
  }, mc.cores = use_cores)
  par_lst
})
names(imp_lst) <- donors

saveRDS(imp_lst, 
        file = "./data/derived/impute_single_mrna.RDS")

