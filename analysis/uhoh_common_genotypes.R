# Goal: Determine for which genotypes endophenotypic data are available and
# ensure that for the hybrid progeny of these parent lines, agronomic data are
# available.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "tidyverse")

smp_conv <- fread("./data/input/tGTP.txt", sep = "\t")
id_gtp <- smp_conv[, .(id_GTP, pool_2), 
                   ][pool_2 %in% c("D", "F"), , ]
id_gtp[, id_GTP := as.character(id_GTP), ]

#--- Predictor data
# Genomic data
snp <- fread("./data/raw/uhoh_smp_id_P235-update.txt")
snp_geno <- id_gtp %>%
  filter(id_GTP %in% snp$id.GTP) %>%
  mutate(Data_Type = "snp")

# There is only one genotype for which we have pedigree records that does not
# also have genomic records. Hence, there is almost no difference between 
# the extent of both datasets and we will reduce the pedigree data set to its
# intersection with the genomic data set to alleviate further coding.
ped <- readRDS("./data/processed/ped-datafull-GTP.RDS") 
ped_geno <- id_gtp %>%
  filter(id_GTP %in% rownames(ped)) %>%
  mutate(Data_Type = "ped")
split_pred <- snp_geno %>%
  rbind(ped_geno) %>%
  split(.$Data_Type) 
common_snp_ped <- split_pred %>% map("id_GTP") %>%
  reduce(intersect)
snp_ped <- split_pred %>%
  map(function(x) x[x$id_GTP %in% common_snp_ped, ]) %>%
  bind_rows() %>%
  as_tibble()


# Transcriptomic data
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna_geno <- id_gtp %>%
  filter(id_GTP %in% colnames(mrna)) %>%
  mutate(Data_Type = "mrna")

# Concatenate predictor data 
pred_geno <- rbindlist(list(snp_ped, mrna_geno))
pred_geno <- pred_geno %>%
  rename(G = id_GTP,
         Pool = pool_2) %>%
  mutate(Pool = ifelse(Pool == "D", yes = "Dent", no = "Flint")) %>%
  as_tibble()

#--- Agronomic traits
traits <- c("GTM", "GTS", "ADF", "FETT", "RPR", "STA", "XZ")
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
# Keep only factorial hybrids and unique records for each trait.
agro_hybrid <- pheno %>%
  as_tibble() %>%
  filter(Trait %in% traits, 
         check == 0,
         dent.GTP %in% (pred_geno %>% filter(Pool == "Dent") %>% .$G),
         flint.GTP %in% (pred_geno %>% filter(Pool == "Flint") %>% .$G)) %>%
  select(G, Trait) %>%
  split(.$Trait) %>%
  map("G") %>%
  reduce(intersect) %>%
  stringr::str_split(pattern = "_") %>%
  map(~matrix(., nrow = 1)) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(V2, V3) %>%
  rename(Dent = V2, Flint = V3) %>%
  mutate(Data_Type = "agro",
         Pool = "Hybrid",
         G = paste0("DF_", Dent, "_", Flint))

agro_inbred <- agro_hybrid %>%
  select(Dent, Flint, Data_Type) %>%
  gather(key = "Pool", value = "G", Dent, Flint) %>%
  distinct()
agro_geno <- rbind(agro_inbred, agro_hybrid %>% select(Data_Type, Pool, G))


#--- Combined data
cmb_data <- pred_geno %>%
  semi_join(y = agro_geno %>% filter(Pool %in% c("Dent", "Flint")),
            by = c("Pool", "G")) %>%
  rbind(agro_geno)

cmb_data %>%
  group_by(Pool, Data_Type) %>%
  count()

saveRDS(cmb_data,
        file = "./data/processed/common_genotypes.RDS",
        compress = FALSE)

