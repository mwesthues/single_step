# !!! Set the working directory to "/silomais_pred2015" !!!
source("./analysis/snp_functions.R")
pacman::p_load("reshape2", "Matrix", "data.table", "matrixcalc", "MASS")


# ----------------------------------------------------------------------------
## Phenotyic data
# Aggregate phenptypic stage 2 BLUEs for all traits.
pheno_fl_nms <- list.files("./data/processed/", pattern = "pheno_BLUE_stage2")
pheno_lst <- lapply(seq_along(pheno_fl_nms), FUN = function(i) {
  fl_nm <- pheno_fl_nms[i]
  trait_pos <- regexpr("(?<=stage2_)[A-Z]+(?=.txt)", text = fl_nm, perl = TRUE)
  trait_nm <- substring(fl_nm, first = trait_pos, 
                        last = trait_pos + attr(trait_pos, "match.length") - 1)
  dat <- fread(paste0("./data/processed/", fl_nm), sep = "\t")
  dat[, Trait := trait_nm, ]
  dat
})
pheno_dt <- rbindlist(pheno_lst)
# Keep only experimental hybrids.
pheno_dt <- pheno_dt[dent.GTP != 0 & flint.GTP != 0, ]
# Names of parental inbred lines.
dent <- pheno_dt[, unique(dent.GTP), ]
saveRDS(dent, compress = FALSE, "./data/processed/dent_nms.RDS")
flint <- pheno_dt[, unique(flint.GTP), ]
saveRDS(flint, compress = FALSE, "./data/processed/flint_nms.RDS")


# MARKER MATRIX PREPARATION -----------------------------------------------
# Load data frame with meta information on marker names and their positions in 
# the genome.
map <- fread("./data/raw/marker_map.txt", sep = "\t")
map <- map[complete.cases(map), ]
map <- map[map$chromosome != 0, ]
map_marker_names <- as.character(map[["markername"]])
map_marker_names <- gsub(pattern = "[.]", replacement = "_",
                         x = map_marker_names)
map_marker_names <- gsub(pattern = "[-]", replacement = "_", 
                         x = map_marker_names)
map[, markername := map_marker_names]


# Load the marker data.
x_marker <- fread("./data/raw/marker_sample_matrix_FWD_P235-update.txt",
                  sep = "\t")

# Load the meta data on samples from the UHOH database.
x_sample <- fread("./data/raw/collect_samples.txt", sep = "\t")

x_sample$id.GTP <- as.factor(x_sample$id.GTP)
x_sample$GTP.DXM.id <- as.factor(x_sample$GTP.DXM.id)
setnames(x_sample, old = "id.GTP", new = "GTP_ID")

# Names for sample ID conversion.
smp_conv <- fread("./data/raw/uhoh_smp_id_P235-update.txt", 
                  sep = "\t")
setnames(x_marker,
         old = colnames(x_marker)[2:ncol(x_marker)], 
         new = as.character(smp_conv[i = match(colnames(x_marker)[
             2:ncol(x_marker)], smp_conv$GTP.DXM.id)][["id.GTP"]]))

# Keep only marker data for inbred lines, which were used throughout the field
# trials.
x_sample[, GTP_ID := as.integer(as.character(GTP_ID)), ]
x_sample <- copy(x_sample[i = GTP_ID %in% c(dent, flint)])
present_genos <- colnames(x_marker)[!colnames(x_marker) %in% "markername"]
dent_nogeno <- dent[!dent %in% present_genos]
flint_nogeno <- flint[!flint %in% present_genos]
pheno_dt <- copy(pheno_dt[i = c(!dent.GTP %in% dent_nogeno &
                                !flint.GTP %in% flint_nogeno)])
par_dent <- as.character(unique(pheno_dt[["dent.GTP"]]))
par_flint <- as.character(unique(pheno_dt[["flint.GTP"]]))
hybrids <- unique(paste0("DF_", pheno_dt[["dent.GTP"]], "_",
                         pheno_dt[["flint.GTP"]]))

# Stop if not all inbred lines are genotyped.
stopifnot(sum(!c(par_dent, par_flint) %in% present_genos) == 0)


#################### SPECIFYING MARKER SELECTION CRITERIA #####################
marker_origin_excl <- NULL
callfreq <- 0.95 # remove if callfreq <= value
rating <- 3
smp_callrate <- 0
maf <- 0.025 # remove if maf <= value
snp_het <- 0.05 # remove if snp_het > 0.05


############################## SNP QUALITY CHECKS ############################
## Keep all samples equal to or smaller than the rating threshold.
rating_pass <- as.character(x_sample[i = x_sample$rating <= 
                                         rating][["GTP_ID"]])

# All marker origins in the UHOH database.
marker_origin <- c('PZE', 'SYN', 'PTU', 'PZA', 'ZM0', 'PZD', 'PZB', 'PHM',
                   'MISC')

if (is.null(marker_origin_excl))
{
    # All markers are included.
    marker_origin_incl <- marker_origin
    
} else
{
    ### Remove markers from sets, which shall not be included.
    unwanted <- unique(grep(pattern = paste(marker_origin_excl,
                                            collapse = "|"),
                            x = x_marker[, "markername"], value = TRUE))
    
    # Markers from origins, which shall be kept.
    x_marker <- x_marker[!as.character(x_marker$markername) %in% unwanted, ]
    
    # Which marker origins are included in the output?
    marker_origin_incl <- marker_origin[!marker_origin %in% marker_origin_excl]
}


# Convert DF with markers to matrix.
init_markmat <- as.matrix(x_marker[, j = colnames(x_marker) %in% rating_pass,
                                   with = FALSE])

# Extract marker names and convert special character to '_'.
mnames <- as.character(x_marker[["markername"]])
mnames <- gsub(pattern = "[.]", replacement = "_", x = mnames)
mnames <- gsub(pattern = "[-]", replacement = "_", x = mnames)

# Add marker names
rownames(init_markmat) <- mnames

# Change all NAs to '?'.
init_markmat[is.na(init_markmat)] <- "??"
init_markmat[init_markmat == "XX"] <- "??"
init_markmat <- t(init_markmat)

# Keep only markers for which there are information in the map object.
init_markmat <- init_markmat[, intersect(colnames(init_markmat), 
                                         map[["markername"]])]
init_markmat <- init_markmat[, match(map$markername, colnames(init_markmat))]
stopifnot(identical(colnames(init_markmat), map$markername))

# Mutual
call_names <- callfreq.fun(x = init_markmat, output = "markerNames",
                           callThresh = callfreq)
het_names <- heterozygosity.fun(x = init_markmat, output = "markerNames",
                                hetThresh = snp_het)
call_het <- intersect(call_names, het_names)
mut_markmat <- init_markmat[, match(call_het, colnames(init_markmat))]

# Get mutual major and minor genotypes for each locus in the combined Dent
# and Flint data.
poly_names <- maf.fun(x = mut_markmat, output = "markerNames")
poly_mat <- mut_markmat[, match(poly_names, colnames(mut_markmat))]
geno_list <- maf.fun(x = poly_mat, output = "genoList")

# Recode dent genotypes.
combi_marker <- recode.fun(x = poly_mat, 
                           major = geno_list[["major_allele"]],
                           minor = geno_list[["minor_allele"]],
                           major_coding = "XX", minor_coding = "YY", 
                           het_coding = "XY", na_coding = "??")

# Check for duplicated markers prior to imputation.
poly_dent <- combi_marker[rownames(combi_marker) %in% par_dent, ]
ncol(poly_dent) - sum(duplicated(x = poly_dent, MARGIN = 2))
poly_flint <- combi_marker[rownames(combi_marker) %in% par_flint, ]
ncol(poly_flint) - sum(duplicated(x = poly_flint, MARGIN = 2))
geno_map <- copy(map[i = match(colnames(combi_marker), map$markername)])

### 2) Preparation of marker data
# Initial marker matrix for heterotic groups.
# Sort matrices and check congruency.
dta_raw <- t(combi_marker)
save(list = c("geno_map", "dta_raw", "maf", "par_dent", "par_flint",
              "hybrids"),
     file = "./data/processed/snp_qc-output.RData")
