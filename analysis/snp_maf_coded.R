pacman::p_load("matrixcalc", "corpcor")
source("./analysis/snp_functions.R")
load("./data/processed/snp_qc-output.RData")
load("./data/processed/imputed_snps.RData")


# KINSHIP MATRICES --------------------------------------------------------
# Dent preparation
dent_polynm <- maf.fun(x = combi_mat[match(par_dent, rownames(combi_mat)), ],
                       mafThresh = maf, output = "markerNames")
d_alph <- combi_mat[match(par_dent, rownames(combi_mat)),
                    match(dent_polynm, colnames(combi_mat))]
d_major <- ref_alleles[["major_allele"]]
d_major <- d_major[match(colnames(d_alph), names(d_major))]
d_minor <- ref_alleles[["minor_allele"]]
d_minor <- d_minor[match(colnames(d_alph), names(d_minor))]
d_num <- recode.fun(x = d_alph, 
                    major = d_major, 
                    minor = d_minor, 
                    major_coding = 1,
                    minor_coding = -1, 
                    het_coding = 999,
                    na_coding = 999)
d_num <- matrix(as.numeric(d_num), 
                nrow = nrow(d_num), ncol = ncol(d_num),
                dimnames = dimnames(d_num))

# Flint preparation.
flint_polynm <- maf.fun(x = combi_mat[match(par_flint, rownames(combi_mat)), ],
                        mafThresh = maf, output = "markerNames")
f_alph <- combi_mat[match(par_flint, rownames(combi_mat)),
                    match(flint_polynm, colnames(combi_mat))]
f_major <- ref_alleles[["major_allele"]]
f_major <- f_major[match(colnames(f_alph), names(f_major))]
f_minor <- ref_alleles[["minor_allele"]]
f_minor <- f_minor[match(colnames(f_alph), names(f_minor))]
f_num <- recode.fun(x = f_alph, 
                    major = f_major, 
                    minor = f_minor, 
                    major_coding = 1,
                    minor_coding = -1, 
                    het_coding = 999,
                    na_coding = 999)
f_num <- matrix(as.numeric(f_num), 
                nrow = nrow(f_num), ncol = ncol(f_num),
                dimnames = dimnames(f_num))

# Recode comb_mat to "2" (major) and "0" (minor), respectively.
snp_mat <- rbind(d_num[, intersect(colnames(d_num), colnames(f_num))],
                 f_num[, intersect(colnames(d_num), colnames(f_num))])
snp_mat[snp_mat == 1] <- 2 
snp_mat[snp_mat == -1] <- 0

saveRDS(snp_mat, file = "./data/processed/snp_mat.RDS")
