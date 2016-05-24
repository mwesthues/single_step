#!/usr/bin/Rscript --vanilla --slave

# !!! Set the working directory to "/silomais_pred2015" !!!

load("./data/processed/snp_qc-output.RData")
source("./analysis/snp_functions.R")
pacman::p_load("data.table")

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
combi_mat <- t(combi_mat[, match(c(par_dent, par_flint), colnames(combi_mat))])
# Preparation for 'maf.fun' function.
combi_mat <- matrix(as.character(combi_mat), 
                    nrow = nrow(combi_mat), ncol = ncol(combi_mat),
                    dimnames = dimnames(combi_mat))
combi_mat[combi_mat == "1"] <- "AA"
combi_mat[combi_mat == "0"] <- "BB"
stopifnot(all(unique(as.vector(combi_mat)) %in% c("AA", "BB")))


# Get the mutual reference allele for each marker across both populations.
ref_alleles <- maf.fun(x = combi_mat, output = "genoList")
ref_names <- maf.fun(x = combi_mat, output = "markerNames")
ref_alleles <- lapply(X = ref_alleles, FUN = function(x){
    names(x) <- ref_names
    x
})

save(list = c("combi_mat", "ref_alleles"), 
     file = "./data/processed/imputed_snps.RData")
