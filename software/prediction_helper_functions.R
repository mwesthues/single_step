# Determine how much time (hh:mm:ss format) has elapsed since script initiation.
get_elapsed_time <- function(start_time, tz = "CEST") {
  start_time <- as.POSIXct(start_time)
  sec <- Sys.time() %>%
    difftime(time1 = ., time2 = start_time, units = "secs") %>%
    lubridate::seconds_to_period()
  paste0(sprintf("%02d", c(day(sec), hour(sec), minute(sec), second(sec))),
         collapse = ":")
}



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
    modify_if(.p = names(.) == "Dent",
           .f = function(mat) {
             append(mat, list(geno = hybrid_parents$Dent))
           }) %>%
    modify_if(.p = names(.) == "Flint",
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



# Function for the computation of the coefficient of variation for each scenario.
coefficient_of_variation <- function(x, y) {
  sqrt(mean(x)) / mean(y)
}




compute_totals_from_fractions <- function(nm, frac) {
  # ---
  # convert a fraction 'frac' * 100% of genotypes, that shall be included in 
  # the training set, into an integer specifying the number of genotypes to be
  # included in the training set.

  # nm: vector of names
  # frac: fractions of genotypes with both, mRNA and SNP information
  #---
  nm %>%
    length() %>%
    map(function(i) i * frac) %>%
    flatten_dbl() %>%
    ceiling()
}





sample_fraction <- function(nm, frac) {
  #---
  # generate a resampled training set provided a vector of genotype names and
  # the fraction of genotypes that shall be included in this training set.

  # nm: names of genotypes
  # frac: vector of genotypes fractions
  #---
  totals <- compute_totals_from_fractions(
    nm = nm,
    frac = frac
  )
  totals %>%
    map(.f = ~sample(
      x = nm,
      size = .,
      replace = FALSE
    )) %>%
  set_names(., nm = as.character(frac)) %>%
  stack()
}





create_t0_trn <- function(hybrids, split_char = "_") {
  #---
  # take a vector of hybrid names ('hybrids'), extract their parental components
  # and, separately for each hybrid, return the corresponding set of hybrids
  # that do not share a parent with the test set hybrid under evaluation

  # hybrids: names of hybrids
  # split_char: separator character for parental components
  #---

  # separate paternal and maternal genotypes and pair them with their offspring
  parent_df <- hybrids %>%
    as_data_frame() %>%
    rename(Hybrid = "value") %>%
    separate(
      Hybrid,
      into = c("Parent_A", "Parent_B"),
      sep = split_char,
      remove = FALSE
    )

  # returns the corresponding t0 hybrids for the test set hybrid under
  # consideration
  select_t0_trn_hybrids <- function(i) {
    parent_vec <- parent_df %>%
      slice(i) %>%
      select(-Hybrid) %>%
      flatten_chr()

    trn_hybrids <- parent_df %>%
      filter(Parent_A != parent_vec[1], Parent_B != parent_vec[2]) %>%
      pull(Hybrid)
    trn_hybrids
  }

  # get the set of t0 hybrid for each test set hybrid
  parent_df %>%
    nrow() %>%
    seq_len() %>%
    map(select_t0_trn_hybrids)
}




sample_loo_sets <- function(geno,
                            frac,
                            iter,
                            material = "Inbred",
                            split_char = "_") {
  #---
  # given a vector of genotype names, that shall be evaluated, generate 'iter'
  # randomly sampled training sets for each genotype, where 'frac' * 100% of the
  # genotypes in 'geno' shall be included in each training set

  # fraction: fraction of genotypes used for the training set
  # geno: names of genotypes
  # iter: number of replicates per genotype
  # material: 'Inbred' or 'Hybrid'
  # split_char: separator for hybrid parents
  #---

  # build the training set by excluding every genotype once from all others
  if (isTRUE(material == "Inbred")) {
    train_geno <- map(geno, .f = ~ setdiff(geno, .))
  } else if (isTRUE(material == "Hybrid")) {
    train_geno <- create_t0_trn(geno, split_char = split_char)
  }
  names(train_geno) <- geno

  # from each training set, sample 'frac' * 100% of genotypes at random and
  # declare them as the resampled training subset. repeat this 'iter' times.
  train_lst <- rerun(.n = iter, {
    sub_train_lst <- map(train_geno, .f = ~sample_fraction(., frac = frac))
    sub_train_df <- sub_train_lst %>% 
      bind_rows(.id = "TST_Geno") %>%
      as_data_frame()
    sub_train <- dplyr::rename(sub_train_df, TRN_Geno = values)
    sub_train
  })

  # concatenate the ten different training sets per test set genotype in one
  # data frame.
  iter_seq <- seq_len(iter)
  names(train_lst) <-  as.character(iter_seq)
  train_df <- bind_rows(train_lst, .id = "Iter")
  
  ### --- Tests ---
  ## test whether the test set genotype is never included in any of the 
  ## training sets
  if (isTRUE(material == "Inbred")) {
    isect <- map2(train_geno, geno, .f = intersect) %>%
      map_int(length) %>%
      unique()
  } else if (isTRUE(material == "Hybrid")) {
    isect <- train_geno %>%
      map(strsplit, split = split_char) %>%
      map(flatten_chr) %>%
      map2(.y = geno, .f = intersect) %>%
      map_int(length) %>%
      unique()
  }
  stopifnot(isect == 0)

  ## test whether the training sets differ in their size, which must be true for
  ## hybrids
  if (isTRUE(material == "Hybrid")) {
    t0_trn_set_size <- train_df %>%
      select(-id) %>%
      group_by(Iter, TST_Geno) %>%
      count() %>%
      pull(n) %>%
      unique() %>%
      length()
      if (isTRUE(t0_trn_set_size == 0)) {
        stop("Hybrid training sets should differ in their sizes")
      }
  }
  ## --- End of tests ---

  # return the data frame with the resampling scheme
  train_df
}
