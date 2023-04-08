#' @export
# remove infections that did not contribute to the focal sample of infections
prune.infection.record <- function(inf_record, focal_infections) {
  `%ni%` <- Negate(`%in%`)
  # Determine seed infections
  seed_infs <- inf_record[grep("seed", inf_record$infector), "inf_id"]
  keep_infs <- c()
  # Loop over each focal infection and trace back until a seed infection is reached
  for (i in focal_infections) {
    j <- i
    j_inf_history <- c()
    while (j %ni% c(seed_infs, keep_infs)) {
      # update j
      j <- inf_record[inf_record$inf_id==j, "origin_inf"]
      # keep this infection in pruned record
      j_inf_history[length(j_inf_history)+1] <- j
    }
    keep_infs[(length(keep_infs)+1):(length(keep_infs)+length(j_inf_history))] <- j_inf_history
  }
  # Keep all focal transmissions
  keep_infs[(length(keep_infs)+1):(length(keep_infs)+length(focal_infections))] <- focal_infections
  # Sort and remove duplicates
  keep_infs <- sort(unique(as.numeric(keep_infs)))
  pruned_inf_record <- inf_record[keep_infs,]
  return(pruned_inf_record)
}

# convert pruned output to format needed by SLiM-based pipeline
generate.SLiM.input <- function(pruned_inf_record, sim_params) {
  output <- as.data.frame(t(sapply(X=unique(pruned_inf_record$origin_inf),
                                   FUN=populate.SLiM.table,
                                   pruned_inf_record=pruned_inf_record)))
  colnames(output) <- c("infection_idx", "infection_source", "infector_type", "infector_id", "infected_ids", "RM_time_start", "inf_time_start")
  rownames(output) <- 1:nrow(output)
  # fill in infection start times
  for (infector_type in c("h", "v")) {
    seed_infs <- intersect(grep("NA", output$inf_time_start), which(output$infector_type==toupper(infector_type)))
    infected_type <- c("v", "h")[-match(infector_type, c("v", "h"))]
    trans_prob_distr <- (1-sim_params[paste0(infector_type, "_rec_rate"), 1])^(1:(sim_params[paste0(infector_type, "_max_duration"), 1])) * sim_params["bite_rate", 1] * sim_params[paste0("raw_", infector_type, infected_type, "_trans_rate"), 1]
    trans_prob_distr[1:(sim_params[paste0(infector_type, "_lag"), 1]-1)] <- 0
    trans_prob_distr <- as.data.frame(cbind(1:sim_params[paste0(infector_type, "_max_duration"), 1], trans_prob_distr));colnames(trans_prob_distr) <- c("x", "prob")
    trans_prob_distr <- pdqr::new_r(trans_prob_distr, type="discrete")
    for (i in seed_infs) {
      start_times <- unlist(strsplit(output[i, "inf_time_start"], split=";"))
      output[i, "inf_time_start"] <- paste(trans_prob_distr(length(start_times)), collapse=";")
    }
  }
  return(output)
}

# 'rewinds' a Ross-Macdonald simulation by removing any events beyond the specified runtime
### TBD - need to modify input/RM parameter table accordingly
rewind.RM <- function(dat, runtime) {
  dat$indiv_status <- dat$indiv_status[,1:runtime]
  dat$infection_record <- dat$infection_record[dat$infection_record$start_t <= runtime,]
  dat$infection_record[dat$infection_record$end_t > runtime, "end_t"] <- NA
  return(dat)
}

# sub-function of generate.SLiM.input()
populate.SLiM.table <- function(pruned_inf_record, inf_id) {
  # Filter to infections that resulted from inf_id
  dat <- pruned_inf_record[which(pruned_inf_record$origin_inf==inf_id),]
  # Determine infector ID
  infector <- unique(dat$infector)
  # Determine who was infected over the course of this infection
  infected <- dat$infected
  # Determine the infector type
  infector_type <- substring(infector,1,1)
  # Determine when these infections began in Ross-Macdonald sim time
  RM_start_times <- as.integer(dat$start_t)
  # If not a seed infection...
  if (length(grep("seed", inf_id))==0) {
    # Determine origin of inf_id
    origin <- pruned_inf_record[pruned_inf_record$inf_id==inf_id,"origin_inf"]
    # Determine infector id
    infector_id <- substring(infector, 2)
    # Determine when the infector was infected in RM time
    infector_start_time <- as.integer(pruned_inf_record[pruned_inf_record$inf_id==inf_id, "start_t"])
    # If a seed infection...
  } else {
    origin <- NA
    infector_id <- "seed"
    infector_start_time <- NA
  }
  # Determine when the infections began wrt to infector's status
  start_times <- suppressWarnings(RM_start_times - infector_start_time)
  # Return output
  output <- c(inf_id,
              origin,
              infector_type,
              infector_id,
              paste(infected, collapse=";"),
              paste(RM_start_times, collapse=";"),
              paste(start_times, collapse=";"))
  return(output)
}
