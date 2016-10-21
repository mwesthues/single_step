# Goal: Generate nine subsets of predictors for later comparisons between their
# predictive abilities:
# a. 100% SNP
# b. 100% pedigree
# c. 42% mRNA
# d. 42% pedigree
# e. 42% SNP
# f. 100% pedigree & 77% SNP
# g. 100% pedigree & 42% mRNA
# h. 100% SNP & 42% mRNA
# i. 100% pedigree & 77% SNP & 42% mRNA
#
# Make comparisons between the following scenarios:
# a, b, f, g, h, i
# c, d, e
#
# Throughout this script, digits behind the variable name of each predictor
# (ped: pedigree, snp: genomic, mrna: transcriptomic) indicate the size of 
# the subset relative to the full set of genotypes (n = 245). For example,
# we have transcriptomic data for 103 genotypes, which are roughly 42% of 
# 245 genotypes, hence the name 'mrna42' and so forth.

pacman::p_load("tidyverse", "data.table")
devtools::install_github("mwesthues/sspredr")



# -- Genotype information.
# Information on the set of genotypes.
genos <- readRDS("./data/processed/common_genotypes.RDS")
genos %>%
  filter(Pool %in% c("Dent", "Flint")) %>%
  group_by(Pool, Data_Type) %>%
  count()

# Get the names of the Dent Hybrid-parents.
geno <- genos %>%
  filter(Pool == "Hybrid") %>%
  .$G %>%
  map(~stringr::str_split(., pattern = "_")) %>%
  at_depth(.depth = 2, .f = 2) %>%
  flatten() %>%
  as_vector()

# -- Load predictor data.
# pedigree
ped100 <- readRDS("./data/processed/ped-datafull-GTP.RDS")
ped_nms <- genos %>%
  filter(Data_Type == "ped") %>%
  .$G
ped100 <- as.matrix(ped100) * 2
ped100 <- ped100 %>%
  .[rownames(.) %in% ped_nms, colnames(.) %in% ped_nms]
ped100[is.na(ped100)] <- 0

# snp
snp100 <- readRDS("./data/processed/snp_mat.RDS")
snp_nms <- genos %>%
  filter(Data_Type == "snp") %>%
  .$G
snp100 <- snp100 %>%
  .[rownames(.) %in% snp_nms, ]

# mrna
mrna42 <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna_nms <- genos %>%
  filter(Data_Type == "mrna") %>%
  .$G
mrna42 <- mrna42 %>%
  as.matrix() %>%
  t() %>%
  .[rownames(.) %in% mrna_nms, ]

# -- Resample predictor data 
# Remove some genomic records so that we can impute them via pedigree
# information.
set.seed(314)
snp35 <- snp100 %>%
  as_tibble() %>%
  mutate(G = rownames(snp100)) %>%
  filter(!G %in% rownames(mrna42)) %>%
  left_join(genos %>% filter(Data_Type == "snp") %>% select(-Data_Type),
            by = "G") %>%
  group_by(Pool) %>%
  sample_frac(size = 0.6) %>%
  ungroup() %>%
  select(-Pool) %>%
  as.data.frame() %>%
  column_to_rownames(var = "G") %>%
  as.matrix()
snp77 <- snp100 %>%
  .[rownames(.) %in% rownames(mrna42), ] %>%
  rbind(snp35)
stopifnot(all(rownames(mrna42) %in% rownames(snp77)))

ped42 <- ped100 %>%
  .[match(rownames(mrna42), rownames(.)), 
    match(rownames(mrna42), colnames(.))]
snp42 <- snp100 %>%
  .[match(rownames(mrna42), rownames(.)), ]

# Store the new objects
storage_dir <- "./data/derived/predictor_subsets/"
ls() %>%
  .[stringr::str_detect(., pattern = "[:digit:]")] %>%
  as.list() %>%
  map_if(.p = ~get(.) %>% is.matrix,
         .f = ~saveRDS(get(.), file = paste0(storage_dir, ., ".RDS"),
                       compress = TRUE))

