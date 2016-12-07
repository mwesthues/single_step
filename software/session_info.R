write_session_info <- function(directory, job_id) {
  # Goal: Store session info for the software environment in which this
  # function is being called.
  require("devtools")
  # Ensure that the path to the storage directory ends with a forward slash, 
  # otherwise add one to the path.
  path_length <- nchar(directory)
  if (!isTRUE(substr(directory, start = path_length,
                     stop = path_length) == "/")) {
    directory <- c(directory, "/")
  }
  writeLines(c(paste("Session info", paste0(rep("-", 68), collapse = "")),
               capture.output(session_info()[[1]]),
               "\n",
               paste("Packages", paste0(rep("-", 68), collapse = "")),
               capture.output(session_info()[[2]])),
             con = paste0(directory, job_id, ".txt"))
}


