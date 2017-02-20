pacman::p_load("sspredr", "tidyverse", "caret", "e1071", "viridis")

# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
hybrids <- genos %>%
  filter(Pool == "Hybrid") %>%
  select(G) %>%
  flatten_chr()



## -- AGONOMIC DATA ----------------------------------------------------
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
# Keep only genotypes that are progeny of parental inbred lines for which at
# least one predictor has records.
pheno %>%
  select(G, EST, Trait) %>%
  filter(G %in% hybrids,
         Trait != "ADF") %>%
  rename(Genotype = G,
         Value = EST) %>%
  as_tibble() %>%
  saveRDS(., "./data/derived/uhoh/agro_tibble.RDS")


## -- TRANSCRIPTOMIC DATA ------------------------------------------------
mrna_inbreds <- genos %>%
  filter(Pool != "Hybrid", Data_Type == "mrna") %>%
  split(.$Pool) %>%
  map("G")

# mRNA
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna_df <- mrna %>%
  as_data_frame() %>%
  mutate(G = rownames(mrna)) %>%
  mutate(Group = ifelse(G %in% mrna_inbreds$Dent, yes = "Dent", no = "Flint"),
         Group = as.factor(Group)) %>%
  select(-G)

# Separately for each heterotic group, run a principal component analysis (PCA) 
# on the "raw" mRNA-BLUEs and return the importance of the first three PCs.
rm_grp <- function(x) {
  x$Group <- NULL
  x
}

## Data pre-processing
# Pre-process the data in the following order:
# 
# 1.    Near-zero-variance filtering.
# 2.    Box-Cox tranformation
# 3.    Centering
# 4.    Scaling
mrna_lst <- mrna_df %>% 
  split(.$Group) %>%
  map(rm_grp) %>%
  map(~as.matrix(.))


# Data transformation
trans_mat <- mrna_lst %>%
  map(preProcess,
      method = c("BoxCox", "center", "scale", "nzv")) %>%
  map2(.y = mrna_lst, .f = predict) %>%
  reduce(rbind)
rownames(trans_mat) <- rownames(mrna)
saveRDS(trans_mat, "./data/derived/uhoh/transformed_mrna.RDS")


## -- PEDIGREE DATA -------------------------------------------------------
ped <- readRDS("./data/processed/ped-datafull-GTP.RDS")
ped_nms <- genos %>%
  filter(Data_Type == "ped") %>%
  .$G
ped <- as.matrix(ped) * 2
ped <- ped %>%
  .[rownames(.) %in% ped_nms, colnames(.) %in% ped_nms]
ped[is.na(ped)] <- 0
if (anyNA(ped)) stop("Missing values in pedigree data not allowed")
saveRDS(ped, "./data/derived/uhoh/pedigree_matrix.RDS")



## -- SNP DATA ------------------------------------------------------------
snp <- readRDS("./data/processed/uhoh/imputed_snp_mat.RDS")
rownames(snp) <- gsub(rownames(snp), pattern = "X", replacement = "")
saveRDS(snp, "./data/derived/uhoh/snp_matrix.RDS")
