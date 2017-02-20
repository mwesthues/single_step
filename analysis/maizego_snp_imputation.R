if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("data.table", "parallel")
pacman::p_load_gh("mwesthues/sspredr")

# MARKER MATRIX PREPARATION -----------------------------------------------
# Load data frame with meta information on marker names and their positions in 
# the genome.
geno_map <- readRDS("./data/derived/maizego/snp_meta_info_datatable.RDS")
setnames(geno_map,
         old = c("rs#", "chrom"),
         new = c("markername", "chr"))
geno_map <- geno_map[, c("markername", "chr", "pos"), with = FALSE]

# Load the marker data.
init_markmat <- readRDS("./data/derived/maizego/tst_raw_snp_matrix.RDS")


#################### SPECIFYING MARKER SELECTION CRITERIA #####################
callfreq <- 0.95 # remove if callfreq <= value
maf <- 0.05 # remove if maf <= value
snp_het <- 0.05 # remove if snp_het > 0.05


############################## SNP QUALITY CHECKS ############################
# Mutual
call_nms <- compute_cf(x = init_markmat, output = "marker_names",
                       call_threshold = callfreq, missing_value = "NN")
het_nms <- compute_het(x = init_markmat, output = "marker_names",
                       het_threshold = snp_het, na_coding = "NN")
init_markmat <- init_markmat[, match(intersect(call_nms, het_nms), 
                                     colnames(init_markmat))]
poly_nms <- compute_maf(x = init_markmat, output = "marker_names", 
                        maf_threshold = maf, missing_value = "NN")
poly_mat <- init_markmat[, match(poly_nms, colnames(init_markmat))]
# Get mutual major and minor genotypes for each locus in the combined Dent
# and Flint data.
geno_lst <- compute_maf(x = poly_mat, output = "geno_list", 
                        missing_value = "NN", maf_threshold = maf)
# Recode genotypes.
combi_marker <- recode_snps(x = poly_mat, 
                            major = geno_lst[["major_genotype"]],
                            minor = geno_lst[["minor_genotype"]],
                            major_coding = "XX", minor_coding = "YY",
                            het_coding = "XY", missing_value = "NN", 
                            na_coding = "NN")

geno_map <- copy(geno_map[i = match(colnames(combi_marker),
                                    geno_map$markername)])
setkey(geno_map, chr, pos)

################################ BEAGLE INPUT FILES ###########################
chromosome <- seq_len(10)
dta_raw <- t(combi_marker)
dta_raw <- dta_raw[geno_map[["markername"]], ]    

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
    map_chr[, c("alleleA", "alleleB", "chr", "markername") := 
                list(rep("X", times = nrow(map_chr)),
                     rep("Y", times = nrow(map_chr)),
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
        write.table(y2, paste0("./data/derived/maizego/snp_qc/imp_input/", "chr", 
                               get("i"), ".txt"),
                    sep = "\t", quote = FALSE, col.names = FALSE, 
                    row.names = FALSE)
        
        # Map file
        write.table(map_chr_mat, 
                    paste0("./data/derived/maizego/snp_qc/imp_input/map_chr",
                           get("i"), ".txt"), sep = "\t", quote = FALSE,
                    col.names = FALSE, row.names = TRUE)
    } else{stop("Mismatched SNP names")}
    message(paste(i, paste0("chr_", i)))
}

############################## BEAGLE OUTPUT FILES ############################
# Extract the names of all Beagle input files.
all_files <- list.files(path = "./data/derived/maizego/snp_qc/imp_input")
chr_length <- 10
use_cores <- 4L

# Make sure to adjust the 'Xmx...m' number, depending on available memory.
# Create Beagle output.
mclapply(seq_len(chr_length), FUN = function(iter) {
    i <- paste0("chr", iter)
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
    system(paste("java -Xmx5000m -jar", paste0(getwd(),
                 "/software/beagle.jar"), 
                 paste0("unphased=", getwd(),
                        "/data/derived/maizego/snp_qc/imp_input/",
                        chr_files[chr_matches]),
                 paste0("markers=", getwd(),
                        "/data/derived/maizego/snp_qc/imp_input/",
                        map_files[map_matches]),
                 "missing=N", 
                 paste0("out=", getwd(),
                        "/data/derived/maizego/snp_qc/imp_output/out"),
                 "niterations=25 nsamples=20"))
}, mc.preschedule = FALSE, mc.cores = use_cores)


# HOMOZYGOUS BEAGLE OUTPUT ------------------------------------------------
gprobs_gz <- list.files(path = "./data/derived/maizego/snp_qc/imp_output/")
gprobs_gz <- gprobs_gz[grepl(pattern = "gprobs", x = gprobs_gz)]

# Unzip all files with genotype probabilities.
for (i in gprobs_gz) {
    
    system(paste("gzip -d -f", 
                 paste0("./data/derived/maizego/snp_qc/imp_output/", i)))
}

# Select all unzipped files containing genotype probabilities.
gprobs_unz <- list.files(path = "./data/derived/maizego/snp_qc/imp_output/")
gprobs_unz <- gprobs_unz[grepl(pattern = "gprobs", x = gprobs_unz)]
gprobs_unz <- gprobs_unz[!grepl(pattern = "gz", x = gprobs_unz)]

mclapply(seq_along(gprobs_unz), FUN = function(iter) {
    i <- gprobs_unz[iter]
    # Specify 'check.names= FALSE' to allow duplicate genotype IDs.
    gprob <- read.table(paste0("./data/derived/maizego/snp_qc/imp_output/", i),
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
                file = paste0("./data/derived/maizego/snp_qc/imp_output/",
                              "homozygous/", outname))
}, mc.preschedule = FALSE, mc.cores = use_cores)


# Concatenate chromosomes from Beagle output.
chromosome <- seq_len(10)
combi_list <- list()
for (i in chromosome) {
    imp_df <- read.table(file = paste0("./data/derived/maizego/snp_qc/", 
                                       "imp_output/homozygous/homoz.chr", i),
                         stringsAsFactors = FALSE, header = TRUE)
    combi_list[[i]] <- imp_df
}
combi_mat <- as.matrix(do.call(rbind, combi_list))
stopifnot(all(combi_mat != "??"))
combi_mat <- t(combi_mat)
# Remove all marker loci that are in perfect linkage disequilibrium. Having
# "copies" of any marker locus will only add noise to our prediction models.
combi_mat <- unique(combi_mat, MARGIN = 2)

# Remove all markers, which might violate the minor allele frequency threshold 
# after imputing missing values.
imp_poly_nms <- compute_maf(x = combi_mat, output = "marker_names", 
                            maf_threshold = maf, missing_value = NA_real_)
snp <- combi_mat[, match(imp_poly_nms, colnames(combi_mat))]
saveRDS(snp, file = "./data/processed/maizego/imputed_snp_mat.RDS")
