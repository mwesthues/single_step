if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "tidyverse", "dtplyr")
log_file <- fread("./data/derived/pred_log.txt")
 
convert_runs <- function(x) {
  x %>%
  str_split(., pattern = "-") %>%
  as_vector() %>%
  as.integer() %>%
  (function(y) seq(from = y[1], to = y[2])) %>%
  length()
}           

compute_runtime <- function(cu, ca, ru, ra, hu, tol = 1) {
  ### Input:
  # cu: number of cores used for development run
  # ca: number of cores used for production run
  # ru: number of development runs
  # ra: number of production runs
  # hu: number of hours needed in development
  # tol: multiplier for time tolerance (1 = no tolerance,
  #                                     1.2 i.e. 20% tolerance)
  ### Output:
  # number of hours needed in production in DD:HH:MM:SS format
  ha <- tol * (cu * ra * hu) / (ru * ca)
  td <- seconds_to_period(ceiling(ha))
  sprintf('%02d:%02d:%02d:%02d', day(td), td@hour, minute(td), second(td))
}

run_df <- log_file %>%
  filter(grepl("2016-11-09|2016-11-08", Date),
         Job_ID != "interactive_00",
         Cores < 16,
         Iter == 30000,
         Job_ID != 10076266,
         Job_ID != 10075993) %>%
  mutate(Elapsed = seconds(hms(Elapsed_Time)),
         Elapsed = str_replace(Elapsed, pattern = "S", replacement = ""),
         Elapsed = as.numeric(Elapsed)) %>%
  rowwise() %>%
  mutate(Runs = convert_runs(Runs),
         Needed = compute_runtime(cu = Cores,
                                  ca = 16,
                                  ru = Runs,
                                  ra = 1521,
                                  hu = Elapsed,
                                  tol = 2))

setup_moab_cmd <- function(pred1, pred2, pred3, trait, iter, model, varcov,
                           pi_param, prior_pi_cnt, walltime, pmem, nodes, 
                           cores) {
  base_cmd <- "msub -v"
  pred_lst <- list(pred1, pred2, pred3)
  pred_lst <- Filter(f = function(x) x != "none", pred_lst)
  pred_lst <- lapply(seq_along(pred_lst), FUN = function(i) {
    paste0("PRED", i, "=", pred_lst[[i]])
  })
  pred_vec <- paste(unlist(pred_lst), collapse = ",")
  trait <- paste0("TRAIT=", trait)
  iter <- paste0("ITER=", iter)
  model <- paste0("MODEL=", model)
  varcov <- paste0("VCOV=", varcov)
  pi_param <- paste0("PI=" , pi_param)
  prior_pi_cnt <- paste0("PRIOR_PI_COUNT=", prior_pi_cnt)
  params <- paste(pred_vec, trait, iter, model, varcov, pi_param, prior_pi_cnt,
                  sep = ",")

  l <- paste0("-l walltime=", walltime, ",", "pmem=", pmem, ",", 
              "nodes=", nodes, ":ppn=", cores)
  script <- "-q singlenode gamazon/moab_scripts/prediction.sh"

  paste(base_cmd, params, l, script)
}

moab_lst <- lapply(seq_len(nrow(run_df)), FUN = function(i) {
  dat <- run_df[i, ]
  setup_moab_cmd(
    pred1 = dat$Pred1, 
    pred2 = dat$Pred2,
    pred3 = dat$Pred3,
    trait = dat$Trait,
    iter = dat$Iter,
    model = dat$Model,
    varcov = "RadenII",
    pi_param = 0.5,
    prior_pi_cnt = 10,
    walltime = dat$Needed,
    pmem = "3900mb",
    nodes = 1,
    cores = 16
  )
})

moab_vec <- apply(X = cbind(moab_lst), MARGIN = 1, paste, collapse = "")
write.table(x = moab_vec,
            file = "./moab_scripts/predictions_comp_time.txt", 
            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)

