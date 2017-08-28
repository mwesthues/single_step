#############################################################################
# Estimate the time required for computations using Bayesian models
#############################################################################
pacman::p_load("tidyverse", "dtplyr", "data.table")

# Times (in fractions of 24 hours) for one cluster (`Interval`) to finish.
# Estimates are based on the variable `template_node_times` in the script
# `./analysis/prediction_template.R`.
cluster_time <- tibble::data_frame(
  Combi = c("CIA", "CIB", "FIA", "FIB", "CHN", "FHN"),
  Day_Fraction = c(0.5, rep(1, times = 4), 2.5)
)

# Load the prediction template to get all clusters.
pred_tmpl <- "./data/derived/prediction_runs/prediction_template.RDS" %>%
  readRDS() %>%
  tibble::as_data_frame() %>%
  dplyr::distinct(Interval) %>%
  tidyr::separate(
    col = Interval,
    into = c("Combi", "Number"),
    remove = FALSE
  ) %>%
  dplyr::inner_join(
    y = cluster_time,
    by = "Combi"
  )


# Convert 'x' days to seconds.
day_to_second <- function(x) {

  x * 24 * 60 * 60
}

# Convert the seconds to 'day:hour:minute:second' format for the MOAB job
# scheduler.
seconds_to_moab_time <- function(in_sec, tolerance) {
  # Purpose: Convert seconds into 'dd:hh:mm:ss' format.
  # sec: Time in seconds (numeric vector element)
  # tolerance: How much tolerance should be permitted? Default is '1'.
  # 1.1 would mean that 10% time is added.
  if (isTRUE(length(in_sec) != 1)) {
    stop("Only a single vector element allowed for arg 'sec'!")
  }
  tol_time <- in_sec * tolerance
  sec <- lubridate::seconds_to_period(tol_time)
  paste0(
    sprintf(
      "%02d",
      c(
        lubridate::day(sec),
        lubridate::hour(sec),
        lubridate::minute(sec),
        lubridate::second(sec)
      )
    ),
    collapse = ":"
  )
}

final_time <- pred_tmpl %>%
  dplyr::mutate(Seconds = day_to_second(Day_Fraction)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Moab_Time = seconds_to_moab_time(Seconds, tolerance = 1.1)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Iter = 30000)


# The computational load for scenario C should be almost identical to that of
# scenario B because they use the same set of genotypes.
# Therefore, we simply copy the parameters from scenario B, recode them and
# declare them as 'Scenario C'.
scenario_c <- final_time %>%
  dplyr::filter(Combi %in% c("CIB", "FIB")) %>%
  dplyr::mutate(Combi = dplyr::case_when(
    Combi == "CIB" ~ "CIC",
    Combi == "FIB" ~ "FIC"
  )) %>%
  dplyr::mutate(Interval = dplyr::case_when(
    grepl("CIB", x = Interval) ~ gsub("CIB", x = Interval, replacement = "CIC"),
    grepl("FIB", x = Interval) ~ gsub("FIB", x = Interval, replacement = "FIC")
  ))

final_time <- dplyr::bind_rows(final_time, scenario_c)

# Define the command line parameters.
moab_msub <- rep("msub -v", times = nrow(final_time))
moab_iter <- paste0("ITER=", final_time$Iter)
moab_interval <- paste0("INTERVAL=", final_time$Interval)
fac_mat <- do.call("cbind", list(moab_iter, moab_interval))
fac_coll <- stringr::str_trim(
  apply(X = fac_mat, MARGIN = 1, paste, collapse = ",")
)

# Create the file with commands.
com_mat <- cbind(
  paste(
    "msub -v",
    fac_coll,
    paste0(
      "-l walltime=",
      final_time$Moab_Time
    ),
    "moab_scripts/prediction.sh"
  )
)

com_args <- apply(X = com_mat, MARGIN = 1, paste, collapse = "")
write.table(
  x = com_args,
  file = "./analysis/moab_commands.txt",
  sep = "\n",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

