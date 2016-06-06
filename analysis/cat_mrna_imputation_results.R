pacman::p_load("data.table", "reshape2")
fl_nms <- list.files("./data/derived/", pattern = "endo_imputation")

grp_lst <- lapply(seq_along(fl_nms), FUN = function(i) {
  fl_nm <- fl_nms[i]
  # Flint or Dent pool
  grp_pos <- regexpr("(?<=Pool=)[[:word:]]+(?=_Model)", text = fl_nm, 
                     perl = TRUE)
  grp_nm <- substring(fl_nm,
                      first = grp_pos,
                      last = grp_pos + attr(grp_pos, "match.length") - 1) 
  # Type of relationship matrix (SNP or Pedigree)
  rel_pos <- regexpr("(?<=RelSource=)[[:word:]]+(?=.RDS)", text = fl_nm,
                     perl = TRUE)
  rel_nm <- substring(fl_nm,
                      first = rel_pos,
                      last = rel_pos + attr(rel_pos, "match.length") - 1)
  dat <- readRDS(paste0("./data/derived/", fl_nm))
  mrna_nms <- sapply(dat, "[[", "mRNA")
  contents <- c("y", "yhat")
  y_lst <- lapply(seq_along(contents), FUN = function(j) {
    content <- contents[j]
    y_df <- as.data.frame(sapply(dat, FUN = "[[", content, 
                                 simplify = TRUE),
                                 stringsAsFactors = FALSE)
    colnames(y_df) <- mrna_nms
    y_df$G <- rownames(y_df)
    y_df$RelSource <- rel_nm
    molten_y <- melt(y_df,
                     id.vars = c("G", "RelSource"),
                     variable.name = "mRNA",    
                     value.name = contents[j],          
                     factorsAsStrings = FALSE)  
    molten_y
  })
  # Merge 'y' and 'yhat'.
  cmb_y_df <- Reduce(function(...) merge(..., all = TRUE), y_lst)
  cmb_y_df$Group <- grp_nm
  cmb_y_df
})
DT <- rbindlist(grp_lst)
setnames(DT, old = "mRNA", new = "Feature")
DT[, EndoClass := ifelse(grepl("MZ00", x = Feature), 
                         yes = "mRNA", no = "Roots"), 
   by = .(Feature)]
# Training set predictive ability
DT[!is.na(y), cor(y, yhat), by = .(Group, EndoClass, RelSource)]
# GTP-IDs of genotypes without mRNAs.
miss_mrna_geno <- DT[is.na(y), unique(G), ]
# GTP-IDs of genotypes with mRNAs.
with_mrna_geno <- DT[!is.na(y), unique(G), ]


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This data.table object still needs to be split according to 'RelSource' and
# 'EndoClass' before it can be transformed into a matrix!
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  
# Get the imputed mRNAs for genotypes without mRNA information only.
yhat_DT <- DT[is.na(y), .(G, mRNA, yhat), ]
yhat_DT <- dcast.data.table(yhat_DT,
                            formula = G ~ mRNA,
                            value.var = "yhat")
yhat_mat <- as.matrix(yhat_DT[, !"G", with = FALSE])  
rownames(yhat_mat) <- yhat_DT[, G, ]
saveRDS(yhat_mat, compress = FALSE,
        file = "./data/processed/imputed_mrna.RDS")
