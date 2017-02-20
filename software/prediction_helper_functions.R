# In the case of hybrid data, we still need to split all predictor matrices
# into Flint and Dent components first.
# In the case of pedigree data, we need to split the data by genotypes in the
# x- as well as the y-dimension because there are no features, which is
# different for other predictor matrices.
split_into_hetgroups <- function(x, y, pedigree = FALSE) {
  if (!pedigree) {
    y %>%
      map(., ~x[rownames(x) %in% ., ])

  } else if (isTRUE(pedigree)) {

    y %>%
      map(., ~x[rownames(x) %in% ., colnames(x) %in% .])
  }
}

# This function will add the names of the parental hybrids to the objects.
# This vector is ordered according to the agronomic data.
# The procedure is necessary for matching predictor data with agronomic data
# throughout all predictions.
add_parental_hybrid_names <- function(x) {
  x %>%
    map_if(.p = names(.) == "Dent",
           .f = function(mat) {
             append(mat, list(geno = hybrid_parents$Dent))
           }) %>%
    map_if(.p = names(.) == "Flint",
           .f = function(mat) {
             append(mat, list(geno = hybrid_parents$Flint))
           })
}


# Create list that names itself based on input object names.
# http://stackoverflow.com/a/16951524/2323832
create_named_list <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
    setNames(L, nm)
}
