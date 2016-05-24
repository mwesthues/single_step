## Reduce the number of markers by selecting roughly equidistant SNPs in 
## increments of 0.1 mega base pairs, which is intended to lessen the 
## computational burden for all Bayesian models.
##
## This approach is based on the following post by Matt Dowle, the creater of
## the data table package: 'http://stackoverflow.com/a/15724500/2323832'
pacman::p_load("data.table")

# Load the marker map with physical positions of loci and adjust the names to 
# match other marker data.
map <- fread("./data/raw/marker_map.txt", sep = "\t")
map <- map[complete.cases(map), ]
map <- map[map$chromosome != 0, ]
map_marker_names <- as.character(map[["markername"]])
map_marker_names <- gsub(pattern = "[.]", replacement = "_",
                         x = map_marker_names)
map_marker_names <- gsub(pattern = "[-]", replacement = "_", 
                         x = map_marker_names)
map[, markername := map_marker_names]

# Keep only markers, which passed quality checks.
snp <- readRDS("./data/processed/snp_mat.RDS")
stopifnot(isTRUE(all(colnames(snp) %in% map$markername)))

# Split the chromosome into intervals of 0.1Mbp width, each. This corresponds 
# to the 10 markers per mega base pair used by Technow et al. (2014, Genetics).
map[, `:=`(Bin = findInterval(position,
                              seq(from = 0, to = max(position), 
                                  by = 1e+05))), by = .(chromosome)]
setkey(map, chromosome, Bin, position)
map <- copy(map[markername %in% colnames(snp), ])
# Add the central, physical position of each bin, which will be required to 
# select roughly equally spaced SNP loci later on.
unq_center <- copy(unique(map[, `:=`(Center = Bin * 1e+05 - 5e+04)][
  , .(chromosome, Bin, Center), ]))
setkey(unq_center, chromosome, Bin, Center)
# Binary search and "roll" to the nearest neighbor.
fin_map <- copy(map[unq_center, roll = "nearest"])
equi_snp <- fin_map[["markername"]]
saveRDS(equi_snp, compress = FALSE, 
        file = "./data/processed/equidistant_snps.RDS")


