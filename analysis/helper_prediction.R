# 300 seconds, one core, 10 rounds, 50,000 iterations
pacman::p_load("lubridate", "stringr", "dplyr", "tidyr")

# Determine the number of LOOCV rounds.
genos <- readRDS("./data/processed/common_genotypes.RDS")
if (length(setdiff(genos$Dent$mrna, genos$Dent$snp)) == 0) {
  dent <- genos$Dent$mrna
} else stop("Not all genotypes have SNP data")
if (length(setdiff(genos$Flint$mrna, genos$Flint$snp)) == 0) {
  flint <- genos$Flint$mrna
} else stop("Not all genotypes have SNP data")
n_rounds <- genos$Hybrid %>%
  as_data_frame() %>%
  separate(value, into = c("DF", "Dent", "Flint")) %>%
  mutate(Avail_Dent = ifelse(Dent %in% dent, yes = 1, no = 0),
         Avail_Flint = ifelse(Flint %in% flint, yes = 1, no = 0),
         Tested_Parents = Avail_Dent + Avail_Flint) %>%
  filter(Tested_Parents == 2) %>%
  nrow()

# The time for one prediction run
elapsed_time <- 300 / 10
n_cores <- 16L
n_traits <- 7L
tolerance <- 1.8
time_dt <- tibble(Elapsed_Time = elapsed_time,
                  Iter = 50000,
                  Trait = "",
                  Model = "BRR",
                  VCOV = "RadenII",
                  Pi = 0.5,
                  PriorPiCount = 10)
time_dt <- time_dt %>%
  mutate(Elapsed_Time = Elapsed_Time / n_cores * n_traits * n_rounds * 
         tolerance)
# Generate various combinations of the fraction of missing values in the Dent
# and Flint mRNA values.
na_cmbs <- expand.grid(dent_na_frac = seq(from = 0.1, to = 0.6, by = 0.1),
                       flint_na_frac = seq(from = 0.1, to = 0.6, by = 0.1))
na_cmbs <- na_cmbs %>%
  as_data_frame() %>%
  mutate(Iter = 50000) %>%
  left_join(time_dt, by = "Iter") %>%
  mutate(Predictor = "mrna")
no_na_cmbs <- na_cmbs %>%
  slice(seq_len(2)) %>%
  mutate(dent_na_frac = 0,
         flint_na_frac = 0,
         Predictor = c("mrna", "snp"))
all_cmbs <- rbind(na_cmbs, no_na_cmbs)

# Convert the seconds to 'day:hour:minute:second' format for the MOAB job
# scheduler.
convert_to_moab_time <- function(x) {
  # Purpose: Convert seconds into 'dd:hh:mm:ss' format.
  # sec: Time in seconds (numeric vector element)
  sec <- lubridate::seconds_to_period(x)
  paste0(sprintf("%02d", c(day(sec), hour(sec), minute(sec), second(sec))),
         collapse = ":")
}
all_cmbs <- all_cmbs %>%
  rowwise() %>%
  mutate(Moab_Time = convert_to_moab_time(Elapsed_Time))
all_cmbs$dent_na_frac <- paste0("DENT_NA_FRACTION=", all_cmbs$dent_na_frac)
all_cmbs$flint_na_frac <- paste0("FLINT_NA_FRACTION=", all_cmbs$flint_na_frac)

# Collapse multiple variables in order to generate MOAB commands.
cmb_mat <- all_cmbs %>%
  select(Trait, Iter, Model, VCOV, Pi, PriorPiCount, Predictor, 
         dent_na_frac, flint_na_frac, Moab_Time) %>%
  rowwise() %>%
  mutate(Trait = paste0("TRAIT=", Trait),
         Iter = paste0("ITER=", Iter),
         Model = paste0("MODEL=", Model),
         VCOV = paste0("VCOV=", VCOV),
         Pi = paste0("PI=", Pi),
         PriorPiCount = paste0("PRIOR_PI_COUNT=", PriorPiCount),
         Predictor = paste0("PREDICTOR=", Predictor)) %>%
  unite(Collapsed, Trait, Iter, Model , VCOV, Pi, PriorPiCount, Predictor,
        dent_na_frac, flint_na_frac, sep = ",") %>%
  select(Collapsed, Moab_Time) %>%
  as.matrix()

# Create the file with commands.
com_mat <- cbind(paste("msub -v", cmb_mat[, "Collapsed"], 
                       paste0("-l walltime=", cmb_mat[, "Moab_Time"],
                              ",pmem=3000mb"),
                       paste0("gamazon/moab_scripts/prediction.sh")))

com_args <- apply(X = com_mat, MARGIN = 1, paste, collapse = "")
write.table(x = com_args,
            file = "./moab_scripts/prediction_est_comp_time.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
