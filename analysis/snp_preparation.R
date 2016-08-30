if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("data.table", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# MARKER MATRIX PREPARATION -----------------------------------------------
# Load data that connect GTP_IDs to pool identifiers (i.e. Dent or Flint)
smp_conv <- fread("./data/input/tGTP.txt", sep = "\t")
id_gtp <- smp_conv[, .(id_GTP, pool_2), 
                   ][pool_2 %in% c("D", "F"), , ]
id_gtp[, id_GTP := as.character(id_GTP), ]

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

# Identify genotypes as 'Dent' or 'Flint'.
snp_dent <- id_gtp[match(smp_conv$id.GTP, id_gtp$id_GTP), , 
                   ][pool_2 == "D", id_GTP, ]
snp_flint <- id_gtp[match(smp_conv$id.GTP, id_gtp$id_GTP), , 
                    ][pool_2 == "F", id_GTP, ]

setnames(x_marker,
         old = colnames(x_marker)[2:ncol(x_marker)], 
         new = as.character(smp_conv[i = match(colnames(x_marker)[
             2:ncol(x_marker)], smp_conv$GTP.DXM.id)][["id.GTP"]]))

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
call_nms <- compute_cf(x = init_markmat, output = "markerNames",
                       callThresh = callfreq)
het_nms <- compute_het(x = init_markmat, output = "markerNames",
                       hetThresh = snp_het)
init_markmat <- init_markmat[, match(intersect(call_nms, het_nms), 
                                     colnames(init_markmat))]
poly_nms <- compute_maf(x = init_markmat, output = "markerNames", 
                        mafThresh = maf)
poly_mat <- init_markmat[, match(poly_nms, colnames(init_markmat))]
# Get mutual major and minor genotypes for each locus in the combined Dent
# and Flint data.
geno_lst <- compute_maf(x = poly_mat, output = "genoList", mafThresh = maf)
# Recode genotypes.
combi_marker <- recode_snps(x = poly_mat, 
                            major = geno_lst[["major_allele"]],
                            minor = geno_lst[["minor_allele"]],
                            major_coding = "XX", minor_coding = "YY",
                            het_coding = "XY", na_coding = "??")

geno_map <- copy(map[i = match(colnames(combi_marker), map$markername)])

### 2) Preparation of marker data
# Initial marker matrix for heterotic groups.
# Sort matrices and check congruency.
dta_raw <- t(combi_marker)
save(list = c("geno_map", "dta_raw", "maf"),
     file = "./51_SNP_qc/_results/snp_qc_output.RData")


# ---
# Store SNP information
smry <- data.frame(Statistic = c("raw_length",
                                 "raw_call_length", 
                                 "raw_het_length",
                                 "raw_callhet_length",
                                 "raw_maf_length",
                                 "imp_maf_dent_length",
                                 "imp_maf_flint_length",
                                 "imp_maf_length",
                                 "imp_equi_length"),
                   Length = c(ncol(init_markmat),
                              length(call_names),
                              length(het_names),
                              length(call_het),
                              length(poly_names),
                              0, 0, 0, 0))
saveRDS(smry, "./54_SNP_Stats/_results/SNP_Summary.RDS", compress = FALSE)


################################ BEAGLE INPUT FILES ###########################
chromosome <- seq_len(10)
setkey(geno_map, chromosome, position)
dta_raw <- dta_raw[geno_map[["markername"]], ]    
setnames(geno_map, old = c("chromosome", "position"), 
         new = c("chr", "pos"))

for (i in chromosome) {
    
    # Marker file of chromosome j
    dta_chr <- dta_raw[geno_map$chr == i, ]    
    # Map file of chromosome j
    map_chr <- geno_map[i = geno_map$chr == i, ]        
    # Data translation for each chromosome
    y <- c()
    
    for (l in 1:ncol(dta_chr)) {
        
        x <- na.omit(dta_chr[, l])
        x1 <- unlist(lapply(X = strsplit(x, ""), FUN = "[[", 1))
        x2 <- unlist(lapply(X = strsplit(x, ""), FUN = "[[", 2))
        y <- cbind(y, x1, x2)
        stopifnot(!any(is.na(y)))
    }
    
    # Generate the first row of the matrix required for Beagle.
    line_names <- rep(colnames(dta_chr), each = 2)
    ifelse(identical(line_names[(1:length(line_names)) %% 2 == 0], 
                     colnames(dta_chr)),
           yes = y1 <- rbind(line_names, y),
           no = warning("Lack of congruency"))
    # Generate the first column required for Beagle.
    I <- c("I", rep("M", nrow(dta_chr)))
    # Generate the second columnn required for Beagle.
    mk_names <- c("id", rownames(dta_chr))
    # Produce the final Beagle input file.
    y2 <- cbind(I, mk_names, y1)        
    upd_map <- copy(map_chr)
    map_chr[, c("alleleA", "alleleB", "chr", "markerindex", "markername") := 
                list(rep("X", times = nrow(map_chr)),
                     rep("Y", times = nrow(map_chr)),
                     NULL,
                     NULL,
                     NULL)]
    map_chr[, pos := as.character(pos)]
    stopifnot(!is.unsorted(as.numeric(map_chr[["pos"]])))
    map_chr_mat <- as.matrix(map_chr)
    rownames(map_chr_mat) <- upd_map[["markername"]]
    colnames(map_chr_mat) <- NULL
    # Output files (== Beagle input files)
    if (identical(rownames(map_chr_mat), as.vector(y2[-1, "mk_names"])))
    {        
        # Marker file
        write.table(y2, paste0("./data/derived/snp_qc/imp_input/", "chr", 
                               get("i"), ".txt"),
                    sep = "\t", quote = FALSE, col.names = FALSE, 
                    row.names = FALSE)
        
        # Map file
        write.table(map_chr_mat, 
                    paste0("./data/derived/snp_qc/imp_input/map_chr",
                           get("i"), ".txt"), sep = "\t", quote = FALSE,
                    col.names = FALSE, row.names = TRUE)
    } else{stop("Mismatched SNP names")}
    message(paste(i, paste0("chr_", i)))
}

############################## BEAGLE OUTPUT FILES ############################
# Extract the names of all Beagle input files.
all_files <- list.files(path = "./data/derived/snp_qc/imp_input")
chr <- paste0("chr", seq_len(10))

# Make sure to adjust the 'Xmx...m' number, depending on available memory.
# Create Beagle output.
for (i in chr) {
    
    chr_files <- all_files[!grepl(pattern = "map", x = all_files)]
    map_files <- all_files[grepl(pattern = "map", x = all_files)]
    
    # Browse the list of geno-files and map-files with matching names,
    # then select the matching geno and map files and impute missing
    # values.
    if (i == "chr1") {
        chr_files <- chr_files[!grepl(pattern = "chr10", x = chr_files)]
        map_files <- map_files[!grepl(pattern = "chr10", x = map_files)]
    }
    chr_matches <- sapply(X = get("i"),
                          FUN = grep, fixed = TRUE, chr_files)
    map_matches <- sapply(X = get("i"),
                          FUN = grep, map_files)
    system(paste("java -Xmx15000m -jar", paste0(getwd(),
                 "/software/beagle.jar"), 
                 paste0("unphased=", getwd(),
                        "/data/derived/snp_qc/imp_input/",
                        chr_files[chr_matches]),
                 paste0("markers=", getwd(),
                        "/data/derived/snp_qc/imp_input/",
                        map_files[map_matches]),
                 "missing=?", 
                 paste0("out=", getwd(),
                        "/data/derived/snp_qc/imp_output/out"),
                 "niterations=25 nsamples=20"))
}


# HOMOZYGOUS BEAGLE OUTPUT ------------------------------------------------
gprobs_gz <- list.files(path = "./data/derived/snp_qc/imp_output/")
gprobs_gz <- gprobs_gz[grepl(pattern = "gprobs", x = gprobs_gz)]

# Unzip all files with genotype probabilities.
for (i in gprobs_gz) {
    
    system(paste("gzip -d -f", 
                 paste0("./data/derived/snp_qc/imp_output/", i)))
}

# Select all unzipped files containing genotype probabilities.
gprobs_unz <- list.files(path = "./data/derived/snp_qc/imp_output/")
gprobs_unz <- gprobs_unz[grepl(pattern = "gprobs", x = gprobs_unz)]
gprobs_unz <- gprobs_unz[!grepl(pattern = "gz", x = gprobs_unz)]

for (i in gprobs_unz) {

    # Specify 'check.names= FALSE' to allow duplicate genotype IDs.
    gprob <- read.table(paste0("./data/derived/snp_qc/imp_output/", i),
                        check.names = FALSE, header = TRUE, row.names = 1)
    
    ### Allele dosage computation.
    # Vector with unique genotypes names.
    geno_names <- unique(names(gprob)[!names(gprob) %in% 
                                          c('alleleA', 'alleleB')])
    # Number of unique genotypes.
    geno_length <- length(unique(names(gprob)[!names(gprob) %in% 
                                                  c('alleleA', 'alleleB')]))
    # Storage matrix for allele 'A' dosage.
    dos_mat <- matrix(rep(NA, times = nrow(gprob) * geno_length),
                      ncol = geno_length, nrow = nrow(gprob),
                      dimnames = list(rownames(gprob), geno_names))
    # Remove 'alleleA' and 'alleleB'.
    X <- gprob[, -c(1, 2)]
    # Split the data into sequences of one genotype with three possible 
    # genotype configurations (i.e. 'XX', 'XY', 'YY'), each.
    seqX <- seq(1, ncol(X), 3)
    # Extract genotype 'j' from the data...
    for (j in 1:length(seqX)) {
        # Retain all markers for this genotypes...
        Y <- X[, seq(seqX[j], seqX[j] + 2)]
        # Compute allele-dosage as ados= 2*p(XX) + p(XY)
        dos_mat[, j] <- 2 * Y[, 1] + Y[, 2]
    }
    # Set genotypes with an allele dosage of 'X' smaller than 1 to '0' and
    # genotypes with an allele dosage of 'X' greater or equal to '1' to '1'.
    dos_mat[dos_mat < 1] <- 0
    dos_mat[dos_mat >= 1] <- 1
    outname <- get("i")
    outname <- gsub(pattern = ".txt.gprobs", replacement = "", x = outname)
    outname <- gsub(pattern = "out.", replacement = "", x = outname)
    outname <- paste0("homoz.", outname)
    write.table(x = dos_mat, 
                file = paste0("./data/derived/snp_qc/imp_output/",
                              "homozygous/", outname))
}


# Concatenate chromosomes from Beagle output.
chromosome <- seq_len(10)
combi_list <- list()
for (i in chromosome) {
    imp_df <- read.table(file = paste0("./data/derived/snp_qc/", 
                                       "imp_output/homozygous/homoz.chr", i),
                         stringsAsFactors = FALSE, header = TRUE)
    combi_list[[i]] <- imp_df
}
combi_mat <- as.matrix(do.call(rbind, combi_list))
colnames(combi_mat) <- gsub(pattern = "X", replacement = "", 
                            x = colnames(combi_mat))
stopifnot(all(combi_mat != "??"))
combi_mat <- t(combi_mat)
# Preparation for 'maf.fun' function.
storage.mode(combi_mat) <- "character"
combi_mat[combi_mat == "1"] <- "AA"
combi_mat[combi_mat == "0"] <- "BB"
stopifnot(all(unique(as.vector(combi_mat)) %in% c("AA", "BB")))


# Get the mutual reference allele for each marker across both populations.
ref_alleles <- compute_maf(x = combi_mat, output = "genoList")
ref_names <- compute_maf(x = combi_mat, output = "markerNames")
ref_alleles <- lapply(X = ref_alleles, FUN = function(x){
    names(x) <- ref_names
    x
})

save(list = c("combi_mat", "ref_alleles"), 
     file = "./data/processed/imputed_snps.RData")


############################## MAF #####################################
# Common genotypes
genos <- readRDS("./data/processed/common_genotypes.RDS")
dent <- genos$Dent$snp
flint <- genos$Flint$snp

dent_polynm <- compute_maf(x = combi_mat[match(dent, rownames(combi_mat)), ],
                           mafThresh = maf, output = "markerNames")
d_alph <- combi_mat[match(dent, rownames(combi_mat)),
                    match(dent_polynm, colnames(combi_mat))]
d_major <- ref_alleles[["major_allele"]]
d_major <- d_major[match(colnames(d_alph), names(d_major))]
d_minor <- ref_alleles[["minor_allele"]]
d_minor <- d_minor[match(colnames(d_alph), names(d_minor))]
d_num <- recode_snps(x = d_alph, 
                     major = d_major, 
                     minor = d_minor, 
                     major_coding = 2,
                     minor_coding = 0, 
                     het_coding = 1,
                     na_coding = 999)

# Flint preparation.
flint_polynm <- compute_maf(x = combi_mat[match(flint, rownames(combi_mat)), ],
                            mafThresh = maf, output = "markerNames")
f_alph <- combi_mat[match(flint, rownames(combi_mat)),
                    match(flint_polynm, colnames(combi_mat))]
f_major <- ref_alleles[["major_allele"]]
f_major <- f_major[match(colnames(f_alph), names(f_major))]
f_minor <- ref_alleles[["minor_allele"]]
f_minor <- f_minor[match(colnames(f_alph), names(f_minor))]
f_num <- recode_snps(x = f_alph, 
                     major = f_major, 
                     minor = f_minor, 
                     major_coding = 2,
                     minor_coding = 0, 
                     het_coding = 1,
                     na_coding = 999)

# Recode comb_mat to "2" (major) and "0" (minor), respectively.
snp_mat <- rbind(d_num[, intersect(colnames(d_num), colnames(f_num))],
                 f_num[, intersect(colnames(d_num), colnames(f_num))])
saveRDS(snp_mat, file = "./data/processed/snp_mat.RDS")
