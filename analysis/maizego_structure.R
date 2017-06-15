if (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "readr")
# Load function to save session info on the software used for the analysis in
# this script.
source("./software/session_info.R")
#


if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("K" = "2")
}

if (!interactive()) {
  job_id <- Sys.getenv("MOAB_JOBID")
} else job_id <- "interactive_00"

# Save session info
write_session_info(directory = "./data/derived/session_info/",
                   job_id = job_id)



# STRUCTURE ---------------------------------------------------------------
# Prepare SNP data
snp_txt_loc <- "./data/processed/maizego/imputed_snp_mat.txt"
snp <- "./data/processed/maizego/imputed_snp_mat.RDS" %>% 
  readRDS() 
genos <- rownames(snp)
snp %>% 
  as_data_frame() %>% 
  write_tsv(., path = snp_txt_loc, col_names = TRUE)

# number of loci
L <- ncol(snp)
# number of individuals
N <- nrow(snp)
# K: maximum number of populations
K <- as.integer(Sys.getenv("K"))


run_structure <- function(K) {
  
  # STRUCTURE parameters
  struc_lib <- paste0(getwd(), "/software/console")
  if (interactive()) {
    soft_loc <- paste0(struc_lib, "/structure")
  } else {
    soft_loc <- paste0(
      "/opt/bwhpc/common/bio/structure/2.3.4/console/structure"
    )
  }
  main_loc <- paste0(struc_lib, "/mainparams")
  extra_loc <- paste0(struc_lib, "/extraparams")
  data_loc <- paste0(getwd(), "/data/processed/maizego/imputed_snp_mat.txt")
  out_loc <- paste0(getwd(), "/data/processed/maizego/popstruc", K)

  struc_call <- paste(
    soft_loc,
    "-m", main_loc,
    "-e", extra_loc,
    "-i", data_loc,
    "-o", out_loc,
    "-K", K,
    "-L", L,
    "-N", N
  )
  # Run the STRUCTURE program.
  system(struc_call)
}

run_structure(K = K)
