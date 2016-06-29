# Goal: Compute the relationship matrix when SNP information are complete
# whereas mRNA records are not.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("cpgen", "foreach", "doMC", "data.table", "matrixStats", 
               "ggplot2")
source("./analysis/snp_functions.R")


predictor <- "mrna"
g_method <- "RadenII"
datall <- readRDS(paste0("~/bwunicluster/silomais_pred2015/91_SI_pedigree/", 
                         "_input/datall.RDS"))

snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(as.matrix(mrna))
endo_lst <- list(mrna = mrna, snp = snp)
comgeno <- Reduce("intersect", lapply(endo_lst, FUN = rownames))
# Remove genotypes for which SNP data are not available.
mrna <- mrna[match(comgeno, rownames(mrna)), ]
# Update the list with endophenotypes.
endo_lst <- list(mrna = mrna, snp = snp)
# Genotypes for which mRNA records are available
n2_nms <- rownames(mrna)
# Genotypes for which mRNA records are unavailable
n1_nms <- setdiff(rownames(snp), n2_nms)

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
comhybrid <- Pheno[grepl(paste0(c(n1_nms, n2_nms), collapse = "|"), 
                         x = Pheno[, G, ]), unique(G), ]
ptlhybrid <- Pheno[dent.GTP %in% c(n1_nms, n2_nms) & flint.GTP %in% 
                   c(n1_nms, n2_nms), ][, G, ]
comhybrid <- intersect(comhybrid, ptlhybrid)
Pheno <- droplevels(Pheno[Pheno$G %in% comhybrid, , ])

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


# Build relationship matrices from the raw features for GBLUP.
rel_lst <- lapply(seq_along(endo_lst), FUN = function(i) {
  M_inbred <- endo_lst[[i]]
  pred_nm <- names(endo_lst)[i]
  # Select inbred lines for which phenotypic records are present in the
  # material and augment the matrices so that they match the number of
  # phenotypic records for the hybrids.
  M_flint <- M_inbred[rownames(M_inbred) %in% comflint, ]
  M_dent <- M_inbred[rownames(M_inbred) %in% comdent, ]
  # Purge all monomorphic features.
  M_flint <- M_flint[, colVars(M_flint) != 0]
  M_dent <- M_dent[, colVars(M_dent) != 0]
  MG_flint <- gmat[[get("g_method")]](M_flint, lambda = 0.001)
  MG_dent <- gmat[[get("g_method")]](M_dent, lambda = 0.001)
  list(Flint = MG_flint, Dent = MG_dent)
})
names(rel_lst) <- names(endo_lst)
# Unlist 'kernels' for compliance with the 'BGLR()' function.
rel_lst <- unlist(rel_lst, recursive = FALSE)
names(rel_lst) <- gsub("[.]", replacement = "_", x = names(rel_lst))


# Generate Design matrices and the inverse of the combined mRNA and SNP
# relationship matrices according to equation (4) from Legarra et al. (2009). 
pred_info <- list(GinvDent = NULL,
                  GinvFlint = NULL,
                  DesignDent = NULL,
                  DesignFlint = NULL,
                  NamesDent = NULL,
                  NamesFlint = NULL,
                  Agronomic = NULL)

# ---------------------------------------------------------------------------
invert_H <- function(x, y, equation = c("Legarra_Eq4", "Legarra_Eq6", 
                                        "Christ_Eq8")) {
  # x: matrix with records on all observations
  # y: matrix with records on a subset of observations
  # equation: which method shall be used for inversion of H
  stopifnot(is.matrix(x))
  stopifnot(is.matrix(y))
  x_nms <- rownames(x)
  nms_2 <- rownames(y)
  nms_1 <- setdiff(x_nms, nms_2)
  n1 <- length(nms_1)
  n2 <- length(nms_2)
  x <- x[match(c(nms_1, nms_2), rownames(x)),
         match(c(nms_1, nms_2), rownames(x))]
  x_11 <- x[nms_1, nms_1]
  x_12 <- x[nms_1, nms_2]
  x_21 <- x[nms_2, nms_1]
  x_22 <- x[nms_2, nms_2]
  stopifnot(identical(dimnames(y), dimnames(x_22)))

  if (isTRUE(equation == "Legarra_Eq4")) {
    H11 <- x_11 + x_12 %*% solve(x_22) %*% (y - x_22) %*% solve(x_22) %*% x_21
    H12 <- x_12 %*% solve(x_22) %*% solve(y)
    H21 <- y %*% solve(x_22) %*% x_21
    H22 <- y
    H <- cbind(rbind(H11, H21), rbind(H12, H22))
    Hinv <- solve(H)
  }
  if (isTRUE(equation == "Legarra_Eq6")) {
    block11 <- x_12 %*% solve(x_22) %*% (y - x_22) %*% solve(x_22) %*% x_21
    block12 <- x_12 %*% solve(x_22) %*% (y - x_22)
    block21 <- (y - x_22) %*% solve(x_22) %*% x_21
    block22 <- y - x_22
    H <- x + cbind(rbind(block11, block21), rbind(block12, block22))
    Hinv <- solve(H)
  }
  if (isTRUE(equation == "Christ_Eq8")) {
    block11 <- matrix(0, nrow = n1, ncol = n1)
    block12 <- matrix(0, nrow = n1, ncol = n2)
    block21 <- matrix(0, nrow = n2, ncol = n1)
    block22 <- solve(y) - solve(x_22)
    Hinv <- solve(x) + cbind(rbind(block11, block21), rbind(block12, block22))
  }
  Hinv
}
# ---------------------------------------------------------------------------
#### Check distribution of coefficients in H^{-1} for different algorithms
# Extract the SNP kernel
snp_krnl <- rel_lst[["snp_Dent"]]
# Extract the mRNA kernel
mrna_krnl <- rel_lst[["mrna_Dent"]]
Hinv_a <- invert_H(x = snp_krnl, y = mrna_krnl, equation = "Legarra_Eq4")
Hinv_b <- invert_H(x = snp_krnl, y = mrna_krnl, equation = "Legarra_Eq6")
Hinv_c <- invert_H(x = snp_krnl, y = mrna_krnl, equation = "Christ_Eq8")
tri_hinv_a <- Hinv_a[upper.tri(Hinv_a, diag = FALSE)]
tri_hinv_b <- Hinv_b[upper.tri(Hinv_b, diag = FALSE)]
tri_hinv_c <- Hinv_c[upper.tri(Hinv_c, diag = FALSE)]
hinv_df <- data.frame(Values = c(tri_hinv_a, tri_hinv_b, tri_hinv_c),
                      Type = c(rep("Legarra_Eq4", times = length(tri_hinv_a)),
                               rep("Legarra_Eq6", times = length(tri_hinv_b)),
                               rep("Christ_Eq8", times = length(tri_hinv_c))))
hinv_df$Type <- as.factor(hinv_df$Type)
ggplot(hinv_df, aes(x = Values, fill = Type)) +
  geom_density(alpha = 0.4) +
  theme(legend.position = "top")


groups <- c("Dent", "Flint")
for (grp in groups) {
  # Extract the SNP kernel
  snp_krnl <- rel_lst[[grep(paste0("snp_", grp), names(rel_lst))]]
  # Extract the mRNA kernel
  mrna_krnl <- rel_lst[[grep(paste0("mrna_", grp), names(rel_lst))]]
  Hinv <- invert_H(x = snp_krnl, y = mrna_krnl, equation = "Legarra_Eq4")
  # Create design matrix for the considered heterotic group
  nm_col <- ifelse(grp == "Dent", yes = 2, no = 3)
  pheno_nms <- sapply(strsplit(rownames(y_mat), split = "_"), 
                      FUN = "[[", nm_col)
  krnl_nms <- rownames(Hinv)
  nm_factor <- factor(pheno_nms, levels = krnl_nms)
  # mapping hybrids to inbreds
  Z <- sparse.model.matrix(~-1 + nm_factor,
                           drop.unused.levels = FALSE)
  colnames(Z) <- gsub("nm_factor", replacement = "", 
                      x = colnames(Z))
  rownames(Z) <- pheno_nms
  pred_info[[grep(paste0("Ginv", grp), x = names(pred_info))]] <- Hinv
  pred_info[[grep(paste0("Design", grp), x = names(pred_info))]] <- Z 
  pred_info[[grep(paste0("Names", grp), x = names(pred_info))]] <- krnl_nms
}
hybs <- rownames(y_mat)
hybs <- gsub("DF_", replacement = "", x = hybs)
pred_info$NamesHybrid <- hybs
pred_info$Agronomic <- y_mat
saveRDS(pred_info, "./data/processed/legarra_inverses.RDS")

