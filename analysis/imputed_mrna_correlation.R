# Genotype-Wise mRNA-Imputation
pacman::p_load("sspredr", "tidyverse", "forcats", "methods", "parallel")

# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
mrna_inbreds <- genos %>%
  filter(Data_Type == "mrna") %>%
  split(.$Pool) %>%
  map("G")

# Remove the 'Group' variable after splitting Dent and Flint data.
rm_grp <- function(x) {
  x$Group <- NULL
  x
}

# mRNA
mrna <- "./data/derived/uhoh/mrna.RDS" %>%
  readRDS() %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  mutate(Group = ifelse(G %in% mrna_inbreds$Dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(~as.matrix(.))

# SNPs 
snp <- "./data/derived/uhoh/snp_matrix.RDS" %>%
  readRDS() %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  filter(G %in% (mrna_inbreds %>% flatten_chr())) %>%
  mutate(Group = ifelse(G %in% mrna_inbreds$Dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(~as.matrix(.)) %>%
  map(., ~sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.025, 
    any_missing = FALSE, remove_duplicated = TRUE
  ))

# Pedigree
ped <- "./data/derived/uhoh/pedigree_matrix.RDS" %>%
  readRDS() %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  filter(G %in% (mrna_inbreds %>% flatten_chr())) %>%
  mutate(Group = ifelse(G %in% mrna_inbreds$Dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  column_to_rownames(var = "G") %>%
  split(.$Group) %>%
  map(rm_grp) %>%
  map(as.matrix) %>%
  map(function(x) x[, match(rownames(x), colnames(x))])



## Function for the imputation of a single mRNA
impute_mrna <- function(x, y, is_pedigree = TRUE) {
  # Input tests
  stopifnot(class(x) == "matrix")
  stopifnot(class(y) == "matrix")
  stopifnot(all(rownames(y) %in% rownames(x)))
  
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

  y <- y[nm2, ]
  M2 <- y[, matrixStats::colVars(y) != 0]
  x <- x[match(c(nm1, nm2), rownames(x)),
         matrixStats::colVars(x) != 0]

  if (isTRUE(is_pedigree)) {
    A <- x[match(unique(geno), rownames(x)),
           match(unique(geno), colnames(x))]
    diag(A) <- diag(A) + 0.01
  } else {
    A <- tcrossprod(scale(x)) / ncol(x)
    diag(A) <- diag(A) + 0.01
  }
  A11 <- A[match(nm1, rownames(A)), match(nm1, colnames(A))]
  A12 <- A[match(nm1, rownames(A)), match(nm2, colnames(A))]
  A21 <- A[match(nm2, rownames(A)), match(nm1, colnames(A))]
  A22 <- A[match(nm2, rownames(A)), match(nm2, colnames(A))]
  Ainv <- solve(A)
  dimnames(Ainv) <- dimnames(A)
  A_up11 <- Ainv[match(nm1, rownames(Ainv)), match(nm1, colnames(Ainv))]
  # Eq.10
  epsilon <- t(chol(solve(A_up11)))

  # Eq.20
  W2 <- Z2 %*% M2
  w1 <- as.matrix(A12 %*% solve(A22) %*% Z2 %*% M2)
  rownames(w1) <- nm1
  w1
}


# Shuffled version of the imputation function above.
shuffle_then_impute <- function(x, y, is_pedigree) {
  cor_lst <- lapply(seq(from = 1, to = nrow(y) - 2), FUN = function(i) {
    idx <- sample(nrow(y), size = i, replace = FALSE)
    z <- y[-idx, ]
    imp_mrna <- impute_mrna(x = x, y = z, is_pedigree = is_pedigree)
    orig_mrna <- y[match(rownames(imp_mrna), rownames(y)), , drop = FALSE]
    stopifnot(identical(rownames(imp_mrna), rownames(orig_mrna)))
    cor(c(imp_mrna), c(orig_mrna))
  })
  cor_lst %>%
    flatten_dbl() %>%
    as_data_frame() %>%
    rownames_to_column(var = "Number_NA_Geno") %>%
    rename(Value = value)
}


## Correlations between vectors of original and imputed mRNAs.
### General idea:
# *    Impute transcriptomic data for different numbers of genotypes (from 1 to
#      number of genotypes with transcriptomic data).
# *    For each number of genotypes without transcriptomic data, randomly sample 
#      1000 sets of genotypes, which will be imputed.
#      This will allow us to aggregate and summarize the correlation between 
#      observed and imputed mRNAs with high confidence.
repeatedly_impute_mrna <- function(x, y, is_pedigree, iter, use_cores = 1L) {
  par_lst <- vector(mode = "list", length = iter)
  names(par_lst) <- paste0("Replication_", seq_len(iter))
  par_lst[] <- mclapply(par_lst, FUN = function(i) {
    try(shuffle_then_impute(x, y, is_pedigree = is_pedigree))
  }, mc.cores = use_cores)
  par_lst
}


ped_imputed_mrna <- map2(
  ped, mrna, .f = repeatedly_impute_mrna, is_pedigree = TRUE, 
  iter = 1000, use_cores = 4
  ) %>%
  map(~keep(., .p = is.data.frame)) %>%
  map(bind_rows, .id = "Replication") %>%
  bind_rows(.id = "Group")

snp_imputed_mrna <- map2(
  snp, mrna, .f = repeatedly_impute_mrna, is_pedigree = FALSE, 
  iter = 1000, use_cores = 4
  ) %>%
  map(~keep(., .p = is.data.frame)) %>%
  map(bind_rows, .id = "Replication") %>%
  bind_rows(.id = "Group")

snp_imputed_mrna %>%
  group_by(Group, Number_NA_Geno) %>%
  summarize(r = mean(Value))

# Pearson correlation coefficients between original and imputed mRNA genotypes 
# for two heterotic groups and two sources of information (i.e. genomic and 
# pedigree), depending on the number of genotypes without mRNA records. 
# Genotypes without mRNA records were declared as such at random.
par_lst <- readRDS("./data/derived/impute_single_mrna.RDS")
par_lst %>%
  map(~keep(., is_list)) %>%
  at_depth(.depth = 2, ~bind_rows(., .id = "Group")) %>%
  map(~bind_rows(., .id = "Replication")) %>%
  bind_rows(.id = "Data_Type") %>%
  mutate(Data_Type = fct_recode(Data_Type,
                                Genomic = "snp_lst",
                                Pedigree = "ped_lst")) %>%
  mutate(Number_NA_Geno = as.factor(Number_NA_Geno)) %>%
  mutate(Number_NA_Geno = fct_inorder(Number_NA_Geno)) %>%
  ggplot(aes(x = Number_NA_Geno, y = Value)) +
  geom_boxplot() +
  facet_wrap(Data_Type ~ Group, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  xlab("Number of Genotypes Without mRNA-Data") +
  ylab("Correlation Between Original and Imputed mRNA BLUEs")
