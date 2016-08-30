# Goal: Determine for which genotypes endophenotypic data are available and
# ensure that for the hybrid progeny of these parent lines, agronomic data are
# available.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table")

smp_conv <- fread("./data/input/tGTP.txt", sep = "\t")
id_gtp <- smp_conv[, .(id_GTP, pool_2), 
                   ][pool_2 %in% c("D", "F"), , ]
id_gtp[, id_GTP := as.character(id_GTP), ]

#--- Genomic data
snp <- fread("./data/processed/uhoh_smp_id_P235-update.txt")
snp_dent <- id_gtp[match(snp$id.GTP, id_gtp$id_GTP), , 
                   ][pool_2 == "D", id_GTP, ]
snp_flint <- id_gtp[match(snp$id.GTP, id_gtp$id_GTP), , 
                    ][pool_2 == "F", id_GTP, ]
snp_geno <- c(snp_dent, snp_flint)


#--- Agronomic traits
traits <- c("GTM", "GTS", "ADF", "FETT", "RPR", "STA", "XZ")
Pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
Pheno <- Pheno[Trait %in% traits, , ]
Pheno <- Pheno[check == 0 & dent.GTP != 0 & flint.GTP != 0, ]
Pheno[, G := as.character(G), ]
hybrid_dt <- Pheno[, .(G = unique(G)), by = .(Trait)]
hybrid_dt[, Trait := as.factor(Trait), ]
hybrid_lst <- lapply(split(hybrid_dt, f = hybrid_dt$Trait),
                     FUN = function(x) {
  x$G
})
hybrids <- Reduce(intersect, hybrid_lst)
Pheno <- Pheno[G %in% hybrids, , ]
agro_dent <- dcast.data.table(Pheno, G ~ Trait, value.var = "dent.GTP")
agro_dent[, G := NULL, ]
agro_dent <- as.list(agro_dent)
agro_dent <- Reduce("intersect", agro_dent)
agro_flint <- dcast.data.table(Pheno, G ~ Trait, value.var = "flint.GTP")
agro_flint[, G := NULL, ]
agro_flint <- as.list(agro_flint)
agro_flint <- Reduce("intersect", agro_flint)


#--- Transcriptomic data
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna_geno <- rownames(mrna)[grep("Exp|Sigma", x = rownames(mrna), 
                                 invert = TRUE)]
mrna_geno <- mrna_geno[mrna_geno %in% id_gtp$id_GTP]
names(mrna_geno) <- id_gtp[match(mrna_geno, id_gtp$id_GTP), pool_2, ]
mrna_dent <- mrna_geno[names(mrna_geno) == "D"]
mrna_flint <- mrna_geno[names(mrna_geno) == "F"]
(n_mrna_dent <- length(mrna_dent))
(n_mrna_flint <- length(mrna_flint))

dent_lst <- list(snp = snp_dent,
                 mrna = mrna_dent,
                 agro = agro_dent)
dent_lst[] <- lapply(dent_lst, FUN = as.integer)
dent_lst[c("snp", "mrna")] <- lapply(dent_lst[c("snp", "mrna")],
                                     FUN = function(x) {
  x[x %in% dent_lst[["agro"]]]
})
dent_lst[["agro"]] <- intersect(dent_lst[["agro"]], dent_lst[["snp"]])

flint_lst <- list(snp = snp_flint,
                  mrna = mrna_flint,
                  agro = agro_flint)
flint_lst[] <- lapply(flint_lst, FUN = as.integer)
flint_lst[c("snp", "mrna")] <- lapply(flint_lst[c("snp", "mrna")],
                                      FUN = function(x) {
  x[x %in% flint_lst[["agro"]]]
})
flint_lst[["agro"]] <- intersect(flint_lst[["agro"]], flint_lst[["snp"]])
hybrids <- Pheno[dent.GTP %in% dent_lst[["agro"]] &
                 flint.GTP %in% flint_lst[["agro"]], unique(G), ]
genos <- list(Dent = dent_lst, Flint = flint_lst, Hybrid = hybrids)

stopifnot(all(Reduce(union, genos$Dent) %in%
              as.integer(sapply(strsplit(genos$Hybrid, split = "_"),
                                FUN = "[[", 2))))
stopifnot(all(Reduce(union, genos$Flint) %in%
              as.integer(sapply(strsplit(genos$Hybrid, split = "_"),
                                FUN = "[[", 3))))
saveRDS(genos, file = "./data/processed/common_genotypes.RDS", compress = FALSE)

