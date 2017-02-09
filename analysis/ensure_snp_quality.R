if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("tidyverse", "data.table", "dtplyr", "viridis", "stringr",
               "corehunter")
pacman::p_load_gh("mwesthues/sspredr")
#source("https://bioconductor.org/biocLite.R")
#biocLite("LEA")
library("LEA")


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
rm(snp)

# Find the common set of genotypes for all data (phenotypic, genotypic, 
# transcriptomic).
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






quality_snps <- ensure_snp_quality(snp = snp_mat,
                                   callfreq_check = TRUE,
                                   callfreq_threshold = 0.95,
                                   maf_check = TRUE,
                                   maf_threshold = 0.05,
                                   any_missing = TRUE,
                                   missing_value = "NN",
                                   remove_duplicated = TRUE)

# Recode marker loci as 0 (no major allele), 0.5 (heterozygous) and 
# 1 (both major alleles).
geno_list <- quality_snps %>%
  sspredr::compute_maf(output = "genoList", mafThresh = 0)

num_snp <- quality_snps %>%
  recode_snps(major = geno_list$major_allele, 
              minor = geno_list$minor_allele,
              major_coding = 1,
              minor_coding = 0,
              het_coding = 0.5,
              na_coding = NA_real_)
saveRDS(num_snp, "./data/derived/maizego/numeric_snp_matrix.RDS")


# Build a genomic relationship matrix.
num_snp <- readRDS("./data/derived/maizego/numeric_snp_matrix.RDS")
G <- build_kernel(M = num_snp)
saveRDS(G, "./data/derived/maizego/genomic_relationship.RDS")
write.lfmm(num_snp, "./data/derived/maizego/snp.lfmm")
rm(no_na_snp, num_snp, G)



# Determine the structure of the data using genotypic and gene expression data,
# respectively.
snp_pc <- LEA::pca(input.file = "./data/derived/maizego/snp.lfmm",
                   scale = TRUE)
# Perform Tracy-Widom tests on all eigenvalues to determine the optimal number
# of principal components.
snp_tw <- tracy.widom(snp_pc)
tw_K <- snp_tw %>%
  .["pvalues"] %>%
  flatten_dbl() %>%
  keep(~ .x <= 0.0001) %>%
  length()
tw_K <- min(c(tw_K, 10))

snp_tw %>%
  .["percentage"] %>%
  flatten_dbl() %>%
  plot()

# Qualitatively assign the genotypes to the previously determined 
# subpopulations.
lfmm2geno("./data/derived/maizego/snp.lfmm", 
          output.file = "./data/derived/maizego/snp.geno")
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
