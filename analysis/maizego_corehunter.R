# Goal: Sample core sets of genotypes for the Yan-lab data and scenario two of
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
options(java.parameters = "-Xmx25G")
pacman::p_load("tidyverse", "data.table", "dtplyr", "corehunter")
pacman::p_load_gh("mwesthues/sspredr")

# Load the names of genotypes, which are covered by all data types (i.e.
# phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/common_snp-mrna-pheno_genotypes.RDS"
)


## --- Scenario 2 SNP preparation.
# Explore the kinship matrix.
snp_path <- "./data/input/maizego/513allchr-imputation.txt"
snp <- snp_path %>%
  fread(header = TRUE) %>%
  as.data.table()

# Collect meta information on the genomic data such as the marker name, the 
# chromosome and the position of the marker on the chromosome.
snp_meta_vars <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                   "center", "protLSID", "assayLSID", "paneLSID", "QCcode")
snp_selection_vars <- c(common_genotypes, snp_meta_vars)
snp_meta_info <- snp %>%
  select(one_of(snp_meta_vars))

# Keep 
snp_mat <- snp %>%
  select(-one_of(snp_meta_vars)) %>%
  select(one_of(common_genotypes)) %>%
  as.matrix() %>%
  t()
colnames(snp_mat) <- snp_meta_info %>%
  select(`rs#`) %>%
  flatten_chr()
snp_geno_nms <- rownames(snp_mat)
rm(snp)

# Remove all marker loci with more than 5% missing values.
high_cf_snp_nms <- snp_mat %>%
  sspredr::compute_cf(., output = "markerNames", missing = "NN", 
                      callThresh = 0.95)
low_na_snp <- snp_mat %>%
  .[, high_cf_snp_nms]
rm(snp_mat)

# Impute missing values in the genotype marker data by replacing each missing 
# marker allele with the most frequent one.
# Select SNPs with a minor allele frequency (MAF) >= 0.05.
low_maf_marker_names <- compute_maf(low_na_snp,
                                    output = "markerNames",
                                    missing = "NN",
                                    mafThresh = 0.05)
low_maf_snp <- low_na_snp[, colnames(low_na_snp) %in% low_maf_marker_names]
major_allele_per_locus <- low_maf_snp %>%
  sspredr::compute_maf(output = "genoList", missing = "NN", 
                       mafThresh = 0) %>%
  .["major_allele"] %>%
  flatten_chr()
stopifnot(identical(length(major_allele_per_locus),
                    ncol(low_maf_snp)))

no_na_snp <- lapply(seq_len(ncol(low_maf_snp)), FUN = function(i) {
  x <- low_maf_snp[, i]
  y <- major_allele_per_locus[i]
  if (any(x == "NN")) {
    x[x == "NN"] <- y
  }
  unique_genotypes <- x %>%
    unique() %>%
    length()
  stopifnot(unique_genotypes == 2)
  x
}) %>%
  combine() %>%
  matrix(., ncol = ncol(low_maf_snp), byrow = FALSE)
rownames(no_na_snp) <- rownames(low_maf_snp)
colnames(no_na_snp) <- colnames(low_maf_snp)
rm(low_na_snp, low_maf_snp)


# Recode marker loci as 0 (no major allele), 0.5 (heterozygous) and 
# 1 (both major alleles).
geno_list <- no_na_snp %>%
  sspredr::compute_maf(output = "genoList", mafThresh = 0)

num_snp <- no_na_snp %>%
  recode_snps(major = geno_list$major_allele, 
              minor = geno_list$minor_allele,
              major_coding = 2,
              minor_coding = 0,
              het_coding = 1,
              na_coding = NA_real_) %>%
  unique(MARGIN = 2)
rm(no_na_snp)

geno <- num_snp %>%
  genotypes(format = "biparental") %>%
  coreHunterData(genotypes = .)
saveRDS(geno, "./data/derived/maizego/scenario2_biparental_genotypes.RDS")
geno <- readRDS("./data/derived/maizego/scenario2_biparental_genotypes.RDS")

core_seq <- seq(from = 0.1, to = 0.9, by = 0.1)
core_lst <- core_seq %>%
  map(., ~sampleCore(data = geno,
           obj = objective(type = "EN", measure = "MR", weight = 1),
           size = .)
  )
names(core_lst) <- paste0("selected_fraction_", core_seq)
saveRDS(core_lst, "./data/derived/maizego/scenario2_snp_core_list.RDS")
