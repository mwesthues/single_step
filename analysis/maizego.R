if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "viridis", "stringr",
               "corehunter")
pacman::p_load_gh("mwesthues/sspredr")
#source("https://bioconductor.org/biocLite.R")
#biocLite("LEA")
library("LEA")


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
snp_path <- "./data/input/maizego/513allchr-imputation.txt"
snp <- snp_path %>%
  fread(header = TRUE) %>%
  as.data.table()


# Collect meta information on the genomic data such as the marker name, the 
# chromosome and the position of the marker on the chromosome.
snp_meta_vars <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                   "center", "protLSID", "assayLSID", "paneLSID", "QCcode")
snp_meta_info <- snp %>%
  select(one_of(snp_meta_vars))

# Keep 
snp_mat <- snp %>%
  select(-`rs#`, -alleles, -chrom, -pos, -strand, -`assembly#`, -center, 
         -protLSID, -assayLSID, -paneLSID, -QCcode) %>%
  as.matrix() %>%
  t()
colnames(snp_mat) <- snp_meta_info %>%
  select(`rs#`) %>%
  flatten_chr()
geno_nms <- rownames(snp_mat)
rm(snp)


# Remove all marker loci with more than 10% missing values.
call_frequency <- 0.1
many_na_markers <- snp_mat %>%
  ncol() %>%
  seq_len() %>%
  map(function(i) {
    mean(snp_mat[, i] == "NN") >= call_frequency
  }) %>%
  flatten_lgl()
mean(many_na_markers)
low_na_snp <- snp_mat %>%
  .[, !many_na_markers] %>%
  as.data.table()
rm(snp_mat)


# Impute missing values in the genotype marker data by replacing each missing 
# marker allele with the most frequent one.
replace_missing_with_major_allele <- function(x) {
  if (any(x == "NN")) {
    major_genotype <- x %>%
      table() %>%
      sort(decreasing = TRUE) %>%
      names() %>%
      .[1]
    x[x == "NN"] <- major_genotype   
  }
  x
}
# !!! Running the following code snippet takes about 30 minutes and requires 
# roughly 20GB of RAM!!!
cols <- colnames(low_na_snp)
low_na_snp[,
           (cols) := lapply(.SD, replace_missing_with_major_allele), 
           .SDcols = cols
           ]
any(low_na_snp == "NN")
saveRDS(low_na_snp, "./data/derived/maizego/imputed_snp_DT.RDS")


# Select SNPs with a minor allele frequency (MAF) >= 0.05.
low_maf_nms <- low_na_snp %>%
  as.matrix() %>%
  sspredr::compute_maf(output = "markerNames", mafThresh = 0.05)
low_maf_snp <- low_na_snp %>%
  as.matrix() %>%
  .[, colnames(.) %in% low_maf_nms]
rm(low_na_snp)


# Recode marker loci as 0 (no major allele), 0.5 (heterozygous) and 
# 1 (both major alleles).
geno_list <- low_maf_snp %>%
  sspredr::compute_maf(output = "genoList", mafThresh = 0)

num_snp <- low_maf_snp %>%
  recode_snps(major = geno_list$major_allele, 
              minor = geno_list$minor_allele,
              major_coding = 1,
              minor_coding = 0,
              het_coding = 0.5,
              na_coding = NA_real_)
rownames(num_snp) <- geno_nms
saveRDS(num_snp, "./data/derived/maizego/numeric_snp_matrix.RDS")


# Build a genomic relationship matrix.
num_snp <- readRDS("./data/derived/maizego/numeric_snp_matrix.RDS")
G <- build_kernel(M = num_snp)
saveRDS(G, "./data/derived/maizego/genomic_relationship.RDS")
write.geno(num_snp, "./data/derived/maizego/snp.geno")
rm(low_maf_snp, num_snp, G)


# Get labels for the different subpopulations of the genotypes.
pop_struc <- fread("./data/input/maizego/genotype_annotation.txt")


# Determine the structure of the data using genotypic and gene expression data,
# respectively.
snp_pc <- LEA::pca(input.file = "./data/derived/maizego/snp.geno",
                   scale = TRUE)
# Perform Tracy-Widom tests on all eigenvalues to determine the optimal number
# of principal components.
snp_tw <- tracy.widom(snp_pc)
tw_K <- snp_tw %>%
  .["pvalues"] %>%
  flatten_dbl() %>%
  keep(~ .x <= 0.01) %>%
  length()
tw_K <- min(c(tw_K, 10))

snp_tw %>%
  .["percentage"] %>%
  flatten_dbl() %>%
  plot()

# Qualitatively assign the genotypes to the previously determined 
# subpopulations.
project <- snmf(input.file = "./data/derived/maizego/snp.geno",
                K = seq_len(tw_K),
                project = "new",
                repetitions = 10,
                CPU = 3,
                entropy = TRUE)

# Based on the results of a PCA, decide whether it is even necessary to treat
# subpopulations as such or whether all inbred lines can be treated as a single
# group.
project <- load.snmfProject("./data/derived/maizego/snp.snmfProject")

pc_mat <- snp_pc$projections
rownames(pc_mat) <- geno_nms
pc_with_assignment <- pc_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "Lines") %>%
  filter(Lines %in% pop_struc$Lines) %>%
  gather(key = PC, value = Score, -Lines) %>%
  left_join(y = pop_struc, by = "Lines") %>%
  mutate(PC = str_replace_all(PC, pattern = "V", replacement = "PC"))

pc_with_assignment %>%
  filter(PC %in% c("PC1", "PC2")) %>%
  spread(key = PC, value = Score) %>%
  mutate(Lines = as.factor(Lines)) %>%
  ggplot(., aes(x = PC1, y = PC2, color = Subpopulations)) +
  geom_point() 

pc_with_assignment %>%
  filter(PC %in% c("PC1", "PC3")) %>%
  spread(key = PC, value = Score) %>%
  mutate(Lines = as.factor(Lines)) %>%
  ggplot(., aes(x = PC1, y = PC3, color = Subpopulations)) +
  geom_point() 


# Evaluate the population structure when considering the first K principal 
# components.
K <- 10
ce <- cross.entropy(project, K = K)
best <- which.min(ce)
barplot(t(Q(project, K = K, run = best)), col = viridis(n = K))



# Assemble ten core sets with high genetic diversity.


# Evaluate the distribution of gene expression values.


# Apply centering, scaling, near-zero-variance removal and a Box-Cox
# transformation to the gene expression data.
