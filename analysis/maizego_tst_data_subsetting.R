if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "viridis", "stringr",
               "corehunter")
pacman::p_load_gh("mwesthues/sspredr")


# Get labels for the different subpopulations of the genotypes.
pop_struc <- fread("./data/input/maizego/genotype_annotation.txt")

# Keep only genotypes from the tropical/subtropical (TST) subset to avoid 
# hassles with extreme population structure.
tst_genotypes <- pop_struc %>%
  filter(Subpopulations == "TST") %>%
  select(Lines) %>%
  flatten_chr()

# Load the phenotypic BLUPs.
pheno_path <- "./data/input/maizego/blup_traits_final.csv"
pheno <- pheno_path %>%
  read_csv() %>%
  rename(Genotype = `<Trait>`) %>%
  gather(key = Trait, value = Value, -Genotype) %>%
  filter(Genotype %in% tst_genotypes)


# Load the gene expression data.
mrna_path <-  "./data/input/maizego/Expression_kernel_finalNormalized_28850.txt"
mrna <- mrna_path %>%
  fread(header = TRUE) %>%
  gather(key = Genotype, value = Expression, -Gene_ID) %>%
  as.data.table() %>%
  filter(Genotype %in% tst_genotypes)



# Explore the kinship matrix.
snp_path <- "./data/input/maizego/513allchr-imputation.txt"
snp <- snp_path %>%
  fread(header = TRUE) %>%
  as.data.table()


# Collect meta information on the genomic data such as the marker name, the 
# chromosome and the position of the marker on the chromosome.
snp_meta_vars <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                   "center", "protLSID", "assayLSID", "paneLSID", "QCcode")
snp_selection_vars <- c(tst_genotypes, snp_meta_vars)
snp_meta_info <- snp %>%
  select(one_of(snp_meta_vars))

# Keep 
snp_mat <- snp %>%
  select(-one_of(snp_meta_vars)) %>%
  select(one_of(tst_genotypes)) %>%
  as.matrix() %>%
  t()
colnames(snp_mat) <- snp_meta_info %>%
  select(`rs#`) %>%
  flatten_chr()
snp_geno_nms <- rownames(snp_mat)

# Find the common set of genotypes for all data (phenotypic, genotypic, 
# transcriptomic).
# This will be useful for scenario 2, where all predictor sets initially cover
# the same set of genotypes but one predictor will be reduced in increments of
# 10-percentage points by employing core-sampling with uniform coverage of the 
# genetic space.
geno_lst <- list(mrna = mrna, pheno = pheno) %>%
  map(., ~select(., Genotype)) %>%
  map(flatten_chr) %>%
  map(unique)
geno_lst$snp <- snp_geno_nms
common_genotypes <- geno_lst %>%
  reduce(intersect)
saveRDS(common_genotypes,
        "./data/derived/maizego/common_snp-mrna-pheno_genotypes.RDS",
        compress = FALSE)

# Determine the number of overlapping elements between each data type. Then, 
# keep the largest intersect between 'pheno' (agronomic data) and either 'mrna'
# or 'snp' data. This will be useful for scenarios 1 and 3 where one predictor
# covers the same set of genotypes as the agronomic data and the other 
# predictor covers less genotypes.
wide_coverage_pred <- geno_lst %>%
  stack() %>%
  table() %>%
  crossprod() %>%
  .[, "pheno"] %>%
  discard(names(.) == "pheno") %>%
  which.max() %>%
  names()
wide_coverage_names <- geno_lst %>%
  keep(names(.) %in% c(wide_coverage_pred, "pheno")) %>%
  reduce(intersect)

# Save the corresponding data.
pheno %>%
  filter(Genotype %in% wide_coverage_names) %>%
  saveRDS(., "./data/derived/maizego/tst_pheno_tibble.RDS")

snp_mat %>%
  .[rownames(.) %in% wide_coverage_names, ] %>%
  saveRDS(., "./data/derived/maizego/tst_raw_snp_matrix.RDS")

mrna %>%
  saveRDS("./data/derived/maizego/tst_raw_mrna_datatable.RDS")
