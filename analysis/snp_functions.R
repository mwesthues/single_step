gmat <- list(
  Zhang = function(M, lambda = 0.01) {
    if (!typeof(M) %in% c("integer", "double")) {
      stop("The SNP data are neither of type 'integer' nor 'double'")
    }
    if (anyNA(M)) stop("NAs in M not allowed")
    # center
    M_scaled <- scale(M, scale = FALSE)
    # get variances
    vars <- apply(M_scaled, MARGIN = 2, FUN = var) 
    # compute crossproduct and add lambda
    G <- ((1 - lambda) * tcrossprod(M_scaled)) / sum(vars)
    for (i in 1:nrow(G)) G[i, i] <- G[i, i] + lambda
    G
  },
  RadenI = function(M, lambda = 0.01) {
    if (!typeof(M) %in% c("integer", "double")) {
      stop("The SNP data are neither of type 'integer' nor 'double'")
    }
    if (anyNA(M)) stop("NAs in M not allowed")
    Z <- tcrossprod(scale(M, center = TRUE, scale = FALSE))
    denominator <- sum(apply(M, MARGIN = 2, FUN = sd))
    G <- Z / denominator
    for (i in 1:nrow(G)) G[i, i] <- G[i, i] + lambda
    G
  },
  RadenII = function(M, lambda = 0.01, feat_weights = NULL) {
    if (!typeof(M) %in% c("integer", "double")) {
      stop("The SNP data are neither of type 'integer' nor 'double'")
    }
    if (anyNA(M)) stop("NAs in M not allowed")
    # center features and scale ...
    # with the square root of the feature heritability multiplied by the
    # standard deviation of the original feature
    if (!is.null(feat_weights)) {
      weighted <- sapply(seq_len(ncol(M)), FUN = function(i) {
        sqrt(feat_weights[i] * var(M[, i]))
      })
      W <- scale(M, center = TRUE, scale = weighted)
    } else {
      # scale with the standard deviation of original features
      W <- scale(M, center = TRUE, scale = TRUE)
    }     
    # G = WW' / m, where m is the number of features
    G <- tcrossprod(W) / ncol(W)
    # Add a small number to the diagonal elements to avoid singularity in G.
    for (i in 1:nrow(G)) G[i, i] <- G[i, i] + lambda
    G
  }
)



recode.fun <- function(x, major, minor, major_coding, minor_coding, 
                       het_coding, na_coding){
    # This function changes the genotype encoding of matrix elements.
    # x: matrix with genotypes encoded as nucleotides, the content must be
    #   exclusively polymorphic markers (usually output of the 'poly.geno()'
    #   function).
    # major: vector with major genotype
    # minor: vector with minor genotypes
    # major_coding: new encoding of the major genotype
    # minor_coding: new encoding of the minor genotypes
    # het_coding: new encoding of heterozygous genotypes
    # na_coding: new encoding of missing elements
    majorMat <- matrix(major, nrow= 1)
    minorMat <- matrix(minor, nrow= 1)
    x[x == majorMat[rep(1, nrow(x)), ]] <- major_coding
    x[x == minorMat[rep(1, nrow(x)), ]] <- minor_coding
    x[x == "??"] <- na_coding
    x[!x %in% c(get("major_coding"), 
                get("minor_coding"), get("na_coding"))] <- het_coding
    x
}




### Find markers with low call call frequencies.
callfreq.fun <- function(x, output = c("markerNames", "markerCallfreq"),
                         callThresh = NULL){
    # class(x) = matrix, rows = genotypes, columns = marker loci
    # Names of all genotypes present in the marker data matrix.
    # output:   a) markerName = vector with markers who passed the criterion
    #           b) markerCallfreq = named vector with call frequencies for 
    #               each vector
    # callThresh: call frequency threshold    
    output <- match.arg(output)
    allgeno <- as.character(na.exclude(unique(c(x))))
    if(!all(nchar(allgeno) == 2)){
        stop(cat("Alleles must be encoded by two letters."))
    }
    stopifnot(all(nchar(as.character(na.exclude(unique(c(x))))) == 2))
    if(all(allgeno != "??")) warning(paste0("Matrix contains now missing", 
                                            "values or they are not encoded",
                                            "as '??'"))
    # Compute the call frequency for each marker.
    callfreq <- 1 - colMeans(x == "??")
    switch(EXPR = output,
           markerNames = colnames(x)[callfreq >= callThresh],
           markerCallfreq = callfreq)
}


### Remove markers with high heterozygosity levels.
heterozygosity.fun <- function(x, output = c("markerNames",
                                             "markerHeterozygosity"),
                               hetThresh = NULL){
    # x: class(matrix), rows = genotypes, columns = marker
    # output:   markerNames = a named vector with markers which pass the
    #           selection criteria
    #           markerHeterozygosity = named vector with heterozygosity levels
    #           for each marker.
    # hetThresh: heterozygosity threshold
    if(is.null(colnames(x))) stop("Assign colnames to x")
    output <- match.arg(output)
    allgeno <- as.character(na.exclude(unique(c(x))))
    if(!all(nchar(allgeno) == 2)){
        stop(cat("Alleles must be encoded by two letters."))
    }
    # Split genotypes into their alleles...
    namesplit <- strsplit(allgeno, split = "")
    # ... and extract heterozygous genotypes.
    het <- allgeno[unlist(lapply(namesplit, function(x) x[1] != x[2]))]
    # Code heterozygous genotypes as '1' and all others as '0'.
    x[x %in% het] <- 1
    x[x != 1] <- 0
    snpHet <- colSums(x == 1) / nrow(x)
    switch(EXPR = output,
           markerHeterozygosity = colSums(x == 1) / nrow(x),
           markerNames = colnames(x[, snpHet <= hetThresh]))
}


maf.fun <- function(x, output = c("markerNames", "markerMAF", "genoList"),
                    missing = "??", mafThresh = 0){
    # x: class = matrix, rows = genotypes, columns = marker
    # output:   a) markerNames = marker, which passed the MAF threshold
    #           b) markerMAF = minor allele frequency for each marker
    #           c) genoList = list with two vectors: (i) major genotype,
    #                           (ii) minor genotype
    # mafThresh: minor allele frequency threshold
    output <- match.arg(output)
    # Names of all genotypes present in the marker data matrix.
    allgeno <- as.character(na.exclude(unique(c(x))))
    if(!all(nchar(allgeno) == 2)){
        stop(cat("Alleles must be encoded by two letters."))
    }
    if(any(allgeno %in% c("XX", "XY", "YY"))){
        stop("Encoding reserved for output")
    }
    # Separate heterozygous genotypes.
    namesplit <- strsplit(allgeno, split = "")
    heteroz <- allgeno[unlist(lapply(namesplit, function(x) x[1] != x[2]))]
    # Separate homozygous genotypes.
    homoz <- allgeno[!allgeno %in% c(missing, heteroz)]
    # Matrix for storage of genotype frequencies.
    genomatrix <- matrix(data= NA, nrow= length(allgeno), 
                         ncol= ncol(x))
    rownames(genomatrix) <- allgeno
    colnames(genomatrix) <- colnames(x)
    # Frequency matrix of genotypes by marker.
    for(i in allgeno){genomatrix[i, ] <- colSums(x == i)}
    # Genotypes excluding missings ("??").
    excna <- genomatrix[rownames(genomatrix) != missing, ]
    # Remove all markers with merely missing values.
    nona <- excna[, colSums(excna) != 0]
    # Remove all monomorphic markers.
    poly <- nona[, colSums(nona != 0) > 1]
    # Reduce the dimensions of the original data frame so that only polymorphic
    # markers are included.
    polyX <- x[, colnames(x) %in% colnames(poly)]
    # Reduce the genotype matrix to homozygous genotypes.
    homozmat <- poly[homoz, ]
    if(any(colSums(homozmat != 0) > 2)) stop(paste("markers with more than", 
                                                   "two homozygous genotypes",
                                                   "present"))
    # For each locus, select the major and the minor genotype.    
    n <- nrow(homozmat)
    major <- vapply(X = seq_len(ncol(homozmat)), FUN = function(i){
        x <- homozmat[, i]
        majorGeno <- names(x[x == sort(x, partial = n)[n]])
        if(length(majorGeno) != 1){
            majorGeno <- majorGeno[1]
        }
        return(majorGeno)
    }, FUN.VALUE = "character")
    
    minor <- vapply(X = seq_len(ncol(homozmat)), FUN = function(i){
        x <- homozmat[,i]
        minorGeno <- names(x[x == sort(x, partial = n - 1)[n - 1]])
        if(length(minorGeno) != 1){
            minorGeno <- minorGeno[2]
        }
        return(minorGeno)
    }, FUN.VALUE = "character")
    # For which loci are there no true major and minor genotypes?
    stopifnot(all(major != minor))
    
    majorMat <- matrix(major, nrow = 1)[rep(1, nrow(polyX)), ]
    minorMat <- matrix(minor, nrow = 1)[rep(1, nrow(polyX)), ]
    
    # Augment the matrix with reference/minor genotypes to have the same number 
    # of rows as the input (Dent or Flint). Replace genotypes by new 
    # identifiers for major, minor and heterozygous genotypes, respectively.
    polyX[polyX == majorMat] <- "XX"
    polyX[polyX == minorMat] <- "YY"
    polyX[!polyX %in% c("XX", "YY", missing)] <- "XY"
    n_maj <- colSums(polyX == "XX")
    n_min <- colSums(polyX == "YY")
    n_het <- colSums(polyX == "XY")
    # Number of genotypes, which are not missing.
    n_geno <- n_maj + n_min + n_het
    # Frequency of the major allele at any locus.
    p_maj <- ((2 * n_maj) + n_het) / (2 * n_geno)
    # Frequency of the minor allele at any locus.
    p_min <- 1 - p_maj
    stopifnot(max(p_min) <= 0.5)
    
    if(mafThresh == 0){
        markernames <- colnames(polyX)
    } else{
        markernames <- colnames(polyX)[p_min >= mafThresh]
    }
    
    switch(output,
           markerNames = markernames,
           markerMAF = p_min,
           genoList = list(major_allele = major, minor_allele = minor))   
}


equi.dist.marker <- function(markerMapList, distance, rounds){
    # markerMapList: class = list, length = number of chromosomes, each list
    #               element comprises a data frame with the following 
    #               variables:
    #               markername = character vector with names of markers
    #               position = numeric vector with loci position
    #               chromosome = character vector with chromosome numbers
    # distance: distance (in base pairs) between two neighboring loci
    # rounds: number of times to run this function
    library("data.table")
    roundList <- lapply(X = names(markerMapList), FUN = function(x){
        obj <- markerMapList[[x]]
        snpPos <- obj[, c("markername", "position")]
        snpPos <- snpPos[order(snpPos$position), ]
        # Vector of starting positions to ensure the selection of different, 
        # equi-distant loci in each round.
        startPos <- seq(from = distance/rounds, to = distance, 
                        by = distance/rounds)
        equiList <- lapply(X = seq_len(rounds), FUN = function(i){
            # Chromosome position at which a marker should be picked.
            specPos <- seq(from = startPos[i], to = max(snpPos[, 2]), 
                           by = distance)
            # Actual position of available markers on the chromosome.
            actPos <- snpPos$position
            # Pick the markers, which are closest to the specified positions.
            pos <- data.table(actPos, val = actPos)
            setattr(x = pos, name = "sorted", value = "actPos")
            dtval <- pos[J(specPos), roll = "nearest"]
            equiDistPos <- dtval$val
            equiDistMarker <- as.character(snpPos[snpPos$position %in%
                                                      equiDistPos,
                                                  "markername"])
            names(equiDistMarker) <- rep(i, times = length(equiDistMarker))
            return(equiDistMarker)
        })
        roundsDf <- stack(unlist(equiList))
        names(roundsDf)[names(roundsDf) %in%
                            c("values", "ind")] <- c("SNP", "Round")
        roundsDf$CHR <- as.factor(x)
        return(roundsDf)  
    })
    return(do.call("rbind", roundList))
}


rounds.kinship <- function(roundsMarkerDf, genoMat){
    # roundsMarkerDf:   class = data.frame
    #   variables:  SNP (character) = marker nanme, 
    #               Round (factor) = identifier for set of selected markers,
    #               CHR (factor) = chromosome
    # genoMat:  matrix with genotypes in rows and marker loci in columns
    
    # Split the data frame with markers by round so that the kinship correction
    # of Rincent et al. (2014) can be applied to each of them separately.
    library("Matrix")
    roundList <- split(x = roundsMarkerDf, f = list(roundsMarkerDf$Round))
    kinList <- lapply(X = seq_len(length(roundList)), FUN = function(i){
        if(is.null(colnames(genoMat))) stop("No marker names in genoMat")
        markerMap <- roundList[[i]]
        chrKin <- lapply(X = seq_len(10), FUN = function(j){
            # Rincent et al. (2014)
            marker <- markerMap[markerMap$CHR != j, "SNP"]
            geno <- genoMat[, marker]
            # Yang et al. (2010)
            G <- tcrossprod(scale(geno)) / ncol(geno)
            Gvec <- as.vector(G)
            names(Gvec) <- rep(j, times = length(Gvec))
            return(Gvec)
        })
        # Get the number of records for each list element.
        matLength <- sapply(X = chrKin, FUN = length)
        chrDf <- as.data.frame(do.call("c", chrKin))
        colnames(chrDf) <- "Value"
        chrDf$CHR <- as.factor(as.character(rep(x = seq_len(10), 
                                                times = matLength)))
        chrDf$ID <- unlist(lapply(X = matLength, FUN = seq_len))
        chrDf$Round <- as.factor(as.character(i))
        return(chrDf)
    })
    fullKinDf <- do.call("rbind", kinList)
    return(fullKinDf)
}



rincent14.fun <- function(markerMap, genoMat){
    # markerMap:    class = data.frame, 
    #               variables:  SNP (chr) = vector with marker names
    #                           CHR (factor) = vector with chromosome numbers
    # genoMat:  class = matrix
    #           variables:  rows = genotypes, columns = marker loci
    #
    # Output:   class = (nested) list 
    #           [["KinMatrix"]] = kinship matrix (Rincent et al., 2014)
    #           [["Marker"]] = names of used markers (for kinship matrix)
    library("Matrix")
    lapply(X = seq_len(length(levels(markerMap$CHR))), FUN = function(j){
        # Rincent et al. (2014)
        marker <- markerMap[markerMap$CHR != j, "SNP"]
        geno <- genoMat[, marker]
        # Yang et al. (2010)
        G <- tcrossprod(scale(geno)) / ncol(geno)
        list(KinMatrix = G, Marker = marker)
    })
}



rincent14.sca <- function(markMap, dent, flint, trait_matrix, dentMat,
                          flintMat){
    # markMap:    class = data.frame, 
    #               variables:  SNP (chr) = vector with marker names
    #                           CHR (factor) = vector with chromosome numbers
    #
    # dent: class = vector, contains corresponding parental dent lines of 
    #               hybrids
    #
    # flint: class = vector, contains corresponding parental flint lines of 
    #               hybrids
    # 
    # trait_matrix: class = matrix
    #               variables: rows: hybrid genotypes, columns: trait vectors
    #               row names: hybrid genotypes
    #
    # dentMat:  class = matrix
    #           variables:  rows = dent genotypes, columns = marker loci
    #
    # flintMat: class = matrix
    #           variables: rows = flint genotypes, columns = marker loci
    #
    # Output:   class = (nested) list 
    #           [["dent_inv"]] = Moore-Penrose inverted dent kinship matrix
    #           [["flint_inv"]] = Moore-Penrose inverted flint kinship matrix
    #           [["hyb_inv"]] = Moore-Penrose inverted hybrid kinship matrix
    library("matrixcalc")
    library("Matrix")
    rin <- lapply(X = seq_len(length(levels(markMap$CHR))), FUN = function(i){
        
        marker <- markMap[markMap$CHR != i, "SNP"]
        # Dent
        sel_dent <- dentMat[, marker]
        g_dent <- tcrossprod(scale(sel_dent)) / ncol(sel_dent)
        pd_dent <- as(object = nearPD(x = g_dent, corr = FALSE, 
                                      keepDiag = TRUE,
                                      do2eigen = TRUE,
                                      ensureSymmetry = TRUE,
                                      doSym = TRUE, posd.tol = 1e-3)$mat,
                      Class = "matrix")
        stopifnot(is.positive.definite(pd_dent))
        aug_dent_names <- match(dent, rownames(g_dent))
        aug_g_dent <- pd_dent[aug_dent_names, aug_dent_names]
        # Flint
        sel_flint <- flintMat[, marker]
        g_flint <- tcrossprod(scale(sel_flint)) / ncol(sel_flint)
        pd_flint <- as(object = nearPD(x = g_flint, corr = FALSE, 
                                       keepDiag = TRUE,
                                       do2eigen = TRUE,
                                       ensureSymmetry = TRUE,
                                       doSym = TRUE, posd.tol = 1e-3)$mat,
                       Class = "matrix")
        stopifnot(is.positive.definite(pd_flint))
        aug_flint_names <- match(flint, rownames(g_flint))
        aug_g_flint <- pd_flint[aug_flint_names, aug_flint_names]
        # Hybrids
        g_hybrid <- hadamard.prod(aug_g_dent, aug_g_flint)
        pd_hybrid <- as(object = nearPD(x = g_hybrid, corr = FALSE, 
                                        keepDiag = TRUE,
                                        do2eigen = TRUE,
                                        ensureSymmetry = TRUE,
                                        doSym = TRUE, posd.tol = 1e-3)$mat,
                        Class = "matrix")
        stopifnot(is.positive.definite(pd_hybrid))
        dimnames(pd_hybrid) <- list(rownames(trait_matrix),
                                    rownames(trait_matrix))
        # Invert the genomic relationship matrices for asreml-r.
        dent_inv <- solve(pd_dent)
        dimnames(dent_inv) <- dimnames(pd_dent)
        flint_inv <- solve(pd_flint)
        dimnames(flint_inv) <- dimnames(pd_flint)
        hyb_inv <- solve(pd_hybrid)
        dimnames(hyb_inv) <- dimnames(pd_hybrid)
        list(dent_inv = dent_inv, flint_inv = flint_inv,
             hyb_inv = hyb_inv)  
    })
    names(rin) <- paste0("Chr_", seq_len(length(levels(markMap$CHR))))
    rin
}



rincent14.diag <- function(markMap, dent, flint, trait_matrix, dentMat,
                           flintMat){
    # markMap:    class = data.frame, 
    #               variables:  SNP (chr) = vector with marker names
    #                           CHR (factor) = vector with chromosome numbers
    #
    # dent: class = vector, contains corresponding parental dent lines of 
    #               hybrids
    #
    # flint: class = vector, contains corresponding parental flint lines of 
    #               hybrids
    # 
    # trait_matrix: class = matrix
    #               variables: rows: hybrid genotypes, columns: trait vectors
    #               row names: hybrid genotypes
    #
    # dentMat:  class = matrix
    #           variables:  rows = dent genotypes, columns = marker loci
    #
    # flintMat: class = matrix
    #           variables: rows = flint genotypes, columns = marker loci
    #
    # Output:   class = (nested) list 
    #           [["dent_inv"]] = Moore-Penrose inverted dent kinship matrix
    #           [["flint_inv"]] = Moore-Penrose inverted flint kinship matrix
    #           [["hyb_inv"]] = Moore-Penrose inverted hybrid kinship matrix
    library("matrixcalc")
    library("Matrix")
    rin <- lapply(X = seq_len(length(levels(markMap$CHR))), FUN = function(i){
        
        marker <- markMap[markMap$CHR != i, "SNP"]
        # Dent
        sel_dent <- dentMat[, marker]
        g_dent <- tcrossprod(scale(sel_dent)) / ncol(sel_dent)
        ifelse(is.positive.definite(g_dent),
               yes = g_dent,
               no = diag(g_dent) <- diag(g_dent) + 1e-6)
        stopifnot(is.positive.definite(g_dent))
        aug_dent_names <- match(dent, rownames(g_dent))
        aug_g_dent <- g_dent[aug_dent_names, aug_dent_names]
        # Flint
        sel_flint <- flintMat[, marker]
        g_flint <- tcrossprod(scale(sel_flint)) / ncol(sel_flint)
        ifelse(is.positive.definite(g_flint),
               yes = g_flint,
               no = diag(g_flint) <- diag(g_flint) + 1e-6)
        stopifnot(is.positive.definite(g_flint))
        aug_flint_names <- match(flint, rownames(g_flint))
        aug_g_flint <- g_flint[aug_flint_names, aug_flint_names]
        # Hybrids
        g_hybrid <- hadamard.prod(aug_g_dent, aug_g_flint)
        ifelse(is.positive.definite(g_hybrid),
               yes = g_hybrid,
               no = diag(g_hybrid) <- diag(g_hybrid) + 1e-6)
        stopifnot(is.positive.definite(g_hybrid))
        dimnames(g_hybrid) <- list(rownames(trait_matrix),
                                   rownames(trait_matrix))
        # Invert the genomic relationship matrices for asreml-r.
        dent_inv <- solve(g_dent)
        dimnames(dent_inv) <- dimnames(g_dent)
        flint_inv <- solve(g_flint)
        dimnames(flint_inv) <- dimnames(g_flint)
        hyb_inv <- solve(g_hybrid)
        dimnames(hyb_inv) <- dimnames(g_hybrid)
        list(dent_inv = dent_inv, flint_inv = flint_inv,
             hyb_inv = hyb_inv)  
    })
    names(rin) <- paste0("Chr_", seq_len(length(levels(markMap$CHR))))
    rin
}



technow14_kinship.fun <- function(x, y, ref_allele){
    # Rationale: Compute kinship matrix according to Technow et al. (2014),
    #           Genetics
    # x: class = matrix, rows: genotypes, columns: marker loci,
    #   coding: "AA", "AB", "BB"
    # y: class = matrix, rows: genotypes, columns: marker loci,
    #   coding: 0 (no ref-alleles), 1 (one ref-allele), 2 (two ref-alleles)
    # ref_allele: List of reference alleles [["major_allele"]] and the cor-
    #               responding minor alleles [["minor_allele"]], respectively.
    #
    # Output: list()
    #       1) [["G"]]: Genomic relationship matrix.
    #       2) [["Ginverse"]]: Inverted genomic relationship matrix.
    library("matrixcalc")
    
    stopifnot(identical(colnames(x), colnames(y)))
    # Keep only reference alleles for markers, which are part of the two 
    # marker matrices.
    ref_vec <- ref_allele[["major_allele"]]
    ref_vec <- ref_vec[match(colnames(x), names(ref_vec))]
    stopifnot(identical(names(ref_vec), colnames(x)))
    
    # Compute the frequency of the reference allele at each locus.
    stopifnot(all(as.character(unique(c(x))) %in% c("AA", "AB", "BB")))
    major_freq <- lapply(X = seq_len(ncol(x)), FUN = function(i){
        (sum(x[, i] == ref_allele[["major_allele"]][i]) +
             sum(x[, i] == "AB") / 2) / nrow(x)
    })
    major_freq <- unlist(major_freq)
    
    # Scale the Dent marker matrix according to Technow et al. (2014)
    w_list <- lapply(X = seq_len(ncol(x)), FUN = function(i){
        xuv <- y[, i]
        pv <- major_freq[i]
        wuv <- (xuv - 2 * pv) / sqrt(4 * pv * (1 - pv))
        wuv
    })
    w_mat <- do.call("cbind", w_list)
    
    # Compute the Dent K matrix according to VanRaden (2008)
    g_mat <- tcrossprod(w_mat) / ncol(w_mat)
    stopifnot(is.positive.definite(g_mat))
    ginv_mat <- solve(g_mat)
    list(G = g_mat, Ginverse = ginv_mat)
}



# Function, which generates design matrices without intercepts from factors.
createZmatrixForFactor <- function(f.vector, factor.name = "",
                                   ordered.factor = FALSE){
    # takes a factor as its input and returns a design matrix without 
    # an intercept
    if (sum(is.na(f.vector)) > 0){
        stop("Missing values in factor not allowed.")
    }
    if (class(f.vector) == "data.frame"){
        f.vector <- f.vector[, ncol(f.vector)]
    }
    if (!is.factor(f.vector)){
        f.levels <- unique(as.character(sort(f.vector)))
        f.vector <- factor(f.vector, levels = f.levels, 
                           ordered = ordered.factor)
    }
    temp.data <- data.frame(y = 1:length(f.vector), new.factor = f.vector)
    Zmatrix <- model.matrix(as.formula(y ~ new.factor), data = temp.data)
    Zmatrix[, 1] <- 0
    Zmatrix[f.vector == levels(f.vector)[1], 1] <- 1
    colnames(Zmatrix) <- paste(factor.name, levels(f.vector), sep = "")
    return(Zmatrix)
}




dupfreqm.fun <- function(x, output = c("names", "copynumber")){
    # x: class = matrix, rows = genotypes, columns = marker loci
    #
    # The output is a named vector with the frequency of duplicated columns
    # of the input object 'x'.
    output <- match.arg(output)
    if(is.null(colnames(x))) stop("Name columns for identification")
    # Matrix of unique markers, which have copies in the full matrix.
    dupMat <- unique(x[, duplicated(x, MARGIN = 2)], MARGIN = 2)
    # Frequency
    freq.fun <- function(dup){
        dupFreq <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(i){
            # Matrix full of copies of the i-th duplicated marker. 
            dupI <- matrix(rep(dup[, i], times = ncol(x)),
                           nrow = nrow(x), ncol = ncol(x))
            # Amount of copies of duplicated marker in full matrix.
            sum(colSums(dupI == x) == nrow(x))
        }))
        names(dupFreq) <- colnames(dup)
        return(dupFreq)
    }
    # Names
    name.fun <- function(dup){
        dupMarker <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(i){
            # Matrix full of copies of the i-th duplicated marker. 
            dupI <- matrix(rep(dup[, i], times = ncol(x)),
                           nrow = nrow(x), ncol = ncol(x))
            # Names of marker with identical copies in full matrix.
            colnames(x)[colSums(dupI == x) == nrow(x)]
        }))
        return(dupMarker)
    }
    
    switch(EXPR = output,
           names = name.fun(dupMat),
           copynumber = freq.fun(dupMat))
}



freq.fun <- function(x){
    # x: class = matrix, rownames = genotypes, colnames = markers
    
    if(is.null(colnames(x))) stop("No names for loci were assigned")
    
    # Get genotypes of all duplicated markers.
    dup <- unique(x[, duplicated(x, MARGIN = 2)], MARGIN = 2)
    
    # Return the names of all duplicated markers.
    dupMark <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(i){
        dupI <- matrix(rep(dup[, i], times = ncol(x)),
                       nrow = nrow(x), ncol = ncol(x))
        return(which(colSums(dupI == x) == nrow(x)))
    }))
    
    # Assign a frequency of 1 to all non-duplicated markers.
    unqMark <- colnames(x)[!colnames(x) %in% names(dupMark)]
    unqFreq <- rep(1, times = length(unqMark))
    names(unqFreq) <- unqMark
    
    # Return the frequency of all duplicated marker genotypes.
    dupFreq <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(j){
        dupI <- matrix(rep(dup[, j], times = ncol(x)),
                       nrow = nrow(x), ncol = ncol(x))
        return(sum(colSums(dupI == x) == nrow(x)))
    }))
    names(dupFreq) <- colnames(dup)
    # Return the frequency of all evaluated markers.
    return(c(unqFreq, dupFreq))
}


dup.fun <- function(x, map, copy, mdist){
    # x:    class = matrix, rows = genotypes, columns = snp loci
    # map:  class = data.frame, 
    #       $markername = names of snp loci, 
    #       $position = position of marker on chromosome
    # copy: class = numeric, number of allowed duplicates of a marker genotype
    # mdist:    maximum distance (in Mbp) between all copies of marker genotype
    # 
    # returns names of markers
    dup <- unique(x[, duplicated(x, MARGIN = 2)], MARGIN = 2)    
    copyMarker <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(i){
        dupI <- matrix(rep(dup[, i], times = ncol(x)),
                       nrow = nrow(x), ncol = ncol(x))
        mnames <- colnames(x)[colSums(dupI == x) == nrow(x)]
        dupPos <- marker_map[marker_map$markername %in% mnames, "position"]
        ifelse(max(dupPos) - min(dupPos) < mdist,
               out <- mnames,
               out <- "empty")
        return(out)
    }))
    copyMarker <- copyMarker[copyMarker != "empty"]
    
    dupMarker <- unlist(lapply(X = seq_len(ncol(dup)), FUN = function(i){
        dupI <- matrix(rep(dup[, i], times = ncol(x)),
                       nrow = nrow(x), ncol = ncol(x))
        colnames(x)[colSums(dupI == x) == nrow(x)]
    }))
    uniqueMarker <- colnames(x)[!colnames(x) %in% dupMarker]
    accMarker <- c(copyMarker, uniqueMarker)
    return(accMarker)
}
