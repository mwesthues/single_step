# Goal: Sample core sets of genotypes for the Yan-lab data and scenario two of
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
options(java.parameters = "-Xmx25G")
pacman::p_load("tidyverse", "data.table", "dtplyr", "corehunter", "LEA", 
               "ggthemes")
pacman::p_load_gh("mwesthues/sspredr")

# For the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- common_genotypes %>%
  reduce(intersect)

## --- Scenario 2 SNP preparation.
# Explore the kinship matrix.
snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  .[rownames(.) %in% common_genotypes, ] %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  )

# Rogers distance.
build_rogers_dist_matrix <- function(x) {
  mat <- x %>%
    poppr::rogers.dist() %>%
    matrix(., nrow = nrow(x), ncol = nrow(x))
  diag(mat) <- 0
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  dimnames(mat) <- rerun(.n = 2, rownames(x))
  stopifnot(isSymmetric(mat))
  corehunter::distances(mat)
}
rogers_mat <- snp %>%
  build_rogers_dist_matrix()

# Genotype data.
geno <- snp %>%
  genotypes(format = "biparental")

# Corehunter input data.
geno_dist_data <- coreHunterData(genotypes = geno,
                                 distances = rogers_mat)
saveRDS(geno_dist_data, 
        file = "./data/derived/geno_dist.RDS")

core_seq <- seq(from = 0.1, to = 0.9, by = 0.1)
core_lst <- core_seq %>%
  map(., ~sampleCore(data = geno,
           obj = objective(type = "EN", measure = "MR", weight = 1),
           size = .)
  )
names(core_lst) <- paste0("selected_fraction_", core_seq)

# Add the full set of genotypes that are covered by mRNA data to the list for 
# the equivalent of a "1.0" scenario to exist.
# This is of interest because it simplifies the prediction script.
core_lst[["selected_fraction_1.0"]][["sel"]] <- common_genotypes
core_lst[["selected_fraction_1.0"]][["EN"]][["MR"]] <- NA_real_
saveRDS(core_lst, "./data/derived/maizego/scenario2_snp_core_list.RDS")



# EVALUATE CORES ----------------------------------------------------------
# Get the names of the core members for each core size.
core_members <- "./data/derived/maizego/scenario2_snp_core_list.RDS" %>%
  readRDS() %>%
  map("sel")

# Get the names of the non-core members for each core size.
core_complement <- core_members %>%
  map(., .f = setdiff, x = common_genotypes)

# Evaluate the modified rogers distance for each core size, separately for 
# included and excluded genotypes.
core_statistics <- list(Members = core_members,
     Complement = core_complement) %>%
  at_depth(.depth = 2, .f = corehunter::evaluateCore,
           data = geno_dist_data,
           objective = objective(type = "EN", measure = "MR")) %>%
  map(bind_rows, .id = "Selected_Fraction") %>%
  bind_rows(.id = "Core_Group") %>%
  gather(key = Fraction, value = MR, -Core_Group) %>%
  filter(MR != 1) %>%
  spread(key = Core_Group, value = MR)


# PCA
write.lfmm(snp, "./data/derived/maizego/common_tst_snp.lfmm")
snp_pc <- LEA::pca("./data/derived/maizego/common_tst_snp.lfmm", scale = TRUE)

group_membership <- list(Members = core_members,
     Complement = core_complement) %>%
  map(stack) %>%
  bind_rows(.id = "Core_Group") %>%
  as_data_frame() %>%
  rename(G = values,
         Fraction = ind) %>%
  mutate(Fraction = gsub(pattern = "selected_fraction_", 
                         replacement = "", 
                         x = Fraction))

pc_mat <- snp_pc$projections
rownames(pc_mat) <- rownames(snp)
colnames(pc_mat) <- paste0("PC_", seq_len(ncol(pc_mat)))

pc_df <- pc_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "G") %>%
  gather(key = PC, value = Score, -G) %>%
  as_data_frame() %>%
  filter(PC %in% paste0("PC_", seq_len(2))) %>%
  left_join(y = group_membership, by = "G") %>%
  mutate(
    PC = gsub("PC_", replacement = "PC", x = PC)
  ) %>%
  spread(key = PC, value = Score)
saveRDS(pc_df, "./data/derived/maizego/core_sampling_pca.RDS")

