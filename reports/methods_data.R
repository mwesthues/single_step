# Goal: Gather all data on the number of genotypes and features that were used
# throughout this study in a comprehensive manner so that it will we easy to 
# write and update the manuscript.



if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "dtplyr", "data.table")



# UHOH --------------------------------------------------------------------
# *** ./analysis/common_genotypes.R ***
# Number of Dent and Flint parental inbred lines as well as hybrids per 
# predictor.
common_genotypes <- "./data/processed/common_genotypes.RDS" %>%
  readRDS()
common_genotypes %>%
  group_by(Pool, Data_Type) %>%
  count()

# *** ./analysis/uhoh_data_preparation.R ***
# Number of transcripts."
uhoh_mrna <- "./data/derived/uhoh/mrna.RDS" %>%
  readRDS() 

uhoh_mrna %>%
  ncol()

# Names of parental inbred lines in the UHOH maize hybrid data set that are 
# covered by transcriptomic data.
uhoh_mrna_nms <- rownames(uhoh_mrna)

# Number of hybrids whose parents are covered by transcriptomic data.
"./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  split(.$Trait) %>% 
  map(., ~select(., Genotype)) %>% 
  map(flatten_chr) %>% 
  reduce(intersect) %>% 
  as_data_frame() %>% 
  separate(value, into = c("DF", "Dent", "Flint")) %>% 
  filter(Dent %in% uhoh_mrna_nms,
         Flint %in% uhoh_mrna_nms) %>% 
  dim()

# Number of pedigree records.
"./data/derived/uhoh/pedigree_matrix.RDS" %>%
  readRDS() %>%
  ncol()

# SNP data
# *** ./analysis/snp_preparation.R ***
## Number of SNP marker before applying quality checks.
"./data/raw/marker_sample_matrix_FWD_P235-update.txt" %>%
  data.table::fread() %>%
  nrow()

# *** ./analysis/uhoh_data_preparation.R ***
## MAF: 5%, Call rate: 95%, Heterozygosity rate: 5%
uhoh_imp_snp <- "./data/derived/uhoh/snp_matrix.RDS" %>%
  readRDS()

# Number of SNPs for all genotypes.
uhoh_imp_snp %>%
  ncol()

# *** ./analysis/prediction.R ***
# 1) Number of SNPs for parental inbred lines of all 1,521 hybrids.
# 2) Number of SNPs for parental inbred lines covered by mRNA data.
common_genotypes %>%
  filter(Data_Type %in% c("snp", "mrna")) %>%
  split(.$Pool, .$Data_Type) %>%
  map(., ~split(., .$Data_Type)) %>%
  at_depth(.depth = 2, "G") %>%
  at_depth(.depth = 2, function(x) {
    uhoh_imp_snp[rownames(uhoh_imp_snp) %in% x, ]
  }) %>%
  at_depth(.depth = 2, .f = ~sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )) %>%
  at_depth(.depth = 2, ncol)



# MAIZEGO INBRED LINES ----------------------------------------------------
# *** ./analysis/maizego_tst_data_subsetting.R ***
# Keep only genotypes from the tropical/subtropical (TST) subset to avoid 
# hassles with extreme population structure.
pop_struc <- fread("./data/input/maizego/genotype_annotation.txt")
tst_genotypes <- pop_struc %>%
  filter(Subpopulations == "TST") %>%
  select(Lines) %>%
  flatten_chr()

# For each predictor, load the names of genotypes that it covers.
maizego_genotypes <- "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS" %>%
  readRDS()

# Number of lines with agronomic information.
"./data/derived/maizego/tst_pheno_tibble.RDS" %>%
  readRDS() %>%
  filter(Genotype %in% tst_genotypes,
         Trait %in% c("100grainweight", "cobweight", "Plantheight", 
                      "Silkingtime", "Kernelwidth", "Eardiameter")) %>%
  group_by(Trait) %>%
  count() %>%
  select(n) %>%
  flatten_int() %>%
  reduce(intersect)


# Number of SNPs
"./data/derived/maizego/tst_raw_snp_matrix.RDS" %>%
  readRDS() %>%
  ncol()

# MAF: 5%, Heterozygosity rate: 5%, Call frequency: 95%
maizego_imp_snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS()
  
maizego_imp_snp %>%
  ncol()

# *** ./analysis/prediction.R ***
# 1) Number of SNPs for inbred lines covered by mRNA data.
# 2) Number of SNPs for all inbred lines.
maizego_genotypes %>%
  discard(names(.) == "pheno") %>%
  map(function(x) {
    maizego_imp_snp[rownames(maizego_imp_snp) %in% x, ]
  }) %>%
  map(., .f = ~sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )) %>%
  map(ncol)



# Number of transcripts.
"./data/derived/maizego/tst_raw_mrna_datatable.RDS" %>%
  readRDS() %>%
  group_by(Genotype) %>%
  count() %>%
  select(n) %>%
  flatten_int() %>%
  reduce(intersect)
