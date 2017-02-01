if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "data.table")

# Load the phenotypic BLUPs.
pheno_path <- "./data/input/maizego/blup_traits_final.csv"
pheno <- pheno_path %>%
  read_csv() %>%
  rename(Genotype = `<Trait>`) %>%
  gather(key = Trait, value = Value, -Genotype)

# Load the gene expression data.
mrna_path <-  "./data/input/maizego/Expression_kernel_finalNormalized_28850.txt"
mrna <- mrna_path %>%
  fread(header = TRUE) %>%
  gather(key = Genotype, value = Expression, -Gene_ID) %>%
  as.data.table()

# Check whether genotypes between phenotypic and gene expression data overlap.
common_genotypes <- list(mrna = mrna, pheno = pheno) %>%
  map(., ~select(., Genotype)) %>%
  map(flatten_chr) %>%
  map(unique) %>%
  reduce(intersect)

# Number of overlapping genotypes.
length(common_genotypes)

# Number of genotypes that are part of one source of information additional to
# the common set of genotypes shared between different sources.
list(mrna = mrna, pheno = pheno) %>%
  map(., ~select(., Genotype)) %>%
  map(flatten_chr) %>%
  map(unique) %>%
  map(~setdiff(., common_genotypes)) %>%
  map(length)


# Explore the kinship matrix.
kinship_path <- "./data/input/maizego/513lines_27229snps_kinship_110608.txt"
kinship_path %>%
  fread(header = FALSE)
count.fields(kinship_path)


# Get labels for the different subpopulations of the genotypes.


# Determine the structure of the data using genotypic and gene expression data,
# respectively.


# Based on the results of a PCA, decide whether it is even necessary to treat
# subpopulations as such or whether all inbred lines can be treated as a single
# group.


# Build a genomic relationship matrix.


# Evaluate the distribution of gene expression values.


# Apply centering, scaling, near-zero-variance removal and a Box-Cox
# transformation to the gene expression data.
