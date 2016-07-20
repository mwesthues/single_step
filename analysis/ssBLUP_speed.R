# Goal: Write MOAB-commands for hybrid prediction scripts that the timing from
# test runs into account to generate decent estimates of run time for any
# script.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "lubridate", "stringr")
(systime <- fread("./data/derived/speed_tests.txt", sep = "\t", header = TRUE))

# Unique combinations of data.frames.
expand.grid.df <- function(...) {
  Reduce(function(...) merge(..., by = NULL), list(...))
}
cv_nms <- c("cv800_Flint=80_Dent=100_TS=500.RDS",
            "cv800_Flint=90_Dent=100_TS=500.RDS")
scen_len <- length(cv_nms)
traits <- c("GTM", "GTS", "ADF", "FETT", "RPR", "STA", "XZ")
cv_df <- data.frame(CV_Scheme = cv_nms)
trait_df <- data.frame(Trait = traits)
elapsed_df <- data.frame(Elapsed = systime[, Elapsed, ])
mrg_dt <- data.table(expand.grid.df(cv_df, trait_df, elapsed_df))
systime[, Trait := NULL, ]
setkey(systime, Elapsed)
setkey(mrg_dt, Elapsed)
time_dt <- merge(systime, mrg_dt)
time_df <- data.frame(time_dt)
detach("package:data.table") # otherwise conflicts with 'lubridate' package


# Multiply the computation time of each predictor set by the number of runs 
# 100 (since 10 runs were simulated and 1,000 are intended), the number of
# scenarios (7), a tolerance of 200% (2) and divide it by the number of cores
# on which the computations will be runs (16).
time_df$Elapsed <- time_df$Elapsed * 100 * scen_len * 2 / 16
names(time_df)[names(time_df) == "Elapsed"] <- "Seconds"
runs <- 800

# Convert the seconds to 'day:hour:minute:second' format for the MOAB job
# scheduler.
moab.time.fun <- function(in_sec, tolerance) {
  # Purpose: Convert seconds into 'dd:hh:mm:ss' format.
  # sec: Time in seconds (numeric vector element)
  # tolerance: How much tolerance should be permitted? Default is '1'. 
  # 1.1 would mean that 10% time is added.
  if (isTRUE(length(in_sec) != 1)) {
    stop("Only a single vector element allowed for arg 'sec'!")
  }
  tol_time <- in_sec * tolerance
  sec <- lubridate::seconds_to_period(tol_time)
  paste0(sprintf("%02d", c(day(sec), hour(sec), minute(sec), second(sec))),
         collapse = ":")
}


# Define the maximum duration for a single job. Any job that requires more time
# on the MOAB server will be split into separate sub jobs based on CV-Runs.
time_lim <- 3600 * 72
# If the computation time is larger than one day, split the corresponding jobs
# into separate smaller jobs with less than one day computation time, each.
high_times <- time_df[time_df$Seconds > (time_lim), ]
time_fctr <- ceiling(high_times$Seconds / (time_lim))
high_time_lst <- lapply(seq_len(nrow(high_times)), FUN = function(i) {
  # Copy columns as many times as the 'time_fctr' dictates. This is necessary
  # to define the different sub jobs.
  dat <- high_times[i, ]
  dat <- dat[rep(1, times = time_fctr[i]), ]
  # Adjust the time requirements for each sub job by dividing the total time
  # requirements of the main job by the number of sub jobs.
  dat$Seconds <- dat$Seconds / time_fctr[i]
  # Add definitions for the CV-runs based on the number of sub jobs.
  subjob_cv_runs <- runs / time_fctr[i]
  dat$CV_Runs <- paste(seq(from = 1, to = runs, by = subjob_cv_runs), 
                       seq(from = subjob_cv_runs, to = runs, 
                           by = subjob_cv_runs), sep = "-")
  dat
})
high_time_df <- do.call("rbind", high_time_lst)
rownames(high_time_df) <- NULL

# Determine all jobs, which were split due to computational constraints in the
# original data frame and remove them since their information is now contained 
# in a different data frame with which the original will be merged afterwards.
low_time_df <- time_df[!time_df$Seconds %in%
                         intersect(time_df$Seconds, high_time_df$Seconds), ]
low_time_df$CV_Runs <- runs
full_time_df <- rbind(low_time_df, high_time_df)

# Convert the time from seconds to the 'dd:hh:mm:ss' format used by the MOAB
# scheduler.
full_time_df$MOAB <- unlist(lapply(seq_len(nrow(full_time_df)),
                                   FUN = function(i) {
  dat <- full_time_df[i, "Seconds"]
  moab.time.fun(in_sec = dat,
                tolerance = 1)
}))

# Based on experience, the jobs have different memory requirements. Hence, jobs
# which run for a long period of time will get more memory than other jobs.
full_time_df[full_time_df$CV_Runs != "800", "Proc.Memory"] <- "3500mb"
full_time_df[full_time_df$CV_Runs == "800", "Proc.Memory"] <- "1800mb"

pacman::p_load("data.table")
full_time_DT <- data.table(full_time_df)
setkey(full_time_DT, CV_Scheme, Trait, Imputation, Only_Profiled)
final_time <- data.frame(full_time_DT)


# Define the command line parameters.
moab_msub <- rep("msub -v", times = nrow(final_time))
moab_trait <- paste0("TRAIT=", final_time$Trait)
moab_iter <- paste0("ITER=", final_time$Iter)
moab_model <- paste0("MODEL=", final_time$Model)
moab_vcov <- paste0("VCOV=", final_time$VCOV)
moab_Pi <- paste0("PI=", final_time$Pi)
moab_PriorPiCount <- paste0("PRIOR_PI_COUNT=", final_time$PriorPiCount)
moab_imp <- paste0("IMPUTATION=", final_time$Imputation)
moab_method <- paste0("CV_METHOD=", final_time$CV_Method)
moab_cv_scheme <- paste0("CV_SCHEME=", final_time$CV_Scheme)
moab_speed <- paste0("SPEED_TEST=", rep("FALSE", times = nrow(final_time)))
moab_profile <- paste0("ONLY_PROFILED=", final_time$Only_Profiled)

fac_mat <- do.call("cbind", list(moab_trait, moab_iter, moab_model, moab_vcov,
                                 moab_Pi, moab_PriorPiCount, moab_imp, 
                                 moab_method, moab_cv_scheme, moab_speed,
                                 moab_profile))
fac_coll <- str_trim(apply(X = fac_mat, MARGIN = 1, paste, collapse = ","))

# Create the file with commands.
com_mat <- cbind(paste("msub -v", fac_coll, 
                       paste0("-l walltime=", final_time$MOAB, ",pmem=",
                              final_time$Proc.Memory),
                       paste0("gamazon/moab_scripts/fernando_ssBLUP.sh")))
com_args <- apply(X = com_mat, MARGIN = 1, paste, collapse = "")
write.table(x = com_args,
            file = "./data/derived/cv/moab_prediction_commands.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
