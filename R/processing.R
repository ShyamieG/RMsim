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

#' @export
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

#' @export
# adds sampling events to a Ross-Macdonald simulation
sample.RM <- function(RM_out,
                      time_step,
                      population=c("H", "V", "both"),
                      proportion=NULL,
                      number=NULL,
                      sample_post_lag=TRUE,
                      resample_possible=TRUE) {
  # Check values
  population <- match.arg(population)
  if (is.null(proportion) & is.null(number)) {
    stop("Either the 'proportion' of the infected population or the 'number' of samples must be specified.")
  }  else if (!is.null(proportion) & !is.null(number)) {
    stop("'proportion' and 'number' cannot both be specified.")
  }
  phase <- which(RM_out$input_parameters["t_start", ] <= time_step & RM_out$input_parameters["t_end", ] >= time_step)
  inf_record <- RM_out$infection_record
  if (sample_post_lag) {
    if (population == "both") {
      for (p in c("H", "V")) {
        h_lag <- RM_out$input_parameters["h_lag", phase]
        v_lag <- RM_out$input_parameters["v_lag", phase]
        dat <- rbind(inf_record[intersect(grep("H", inf_record$infected), which(inf_record$start_t + h_lag <= time_step  & inf_record$end_t >= time_step)),], inf_record[intersect(grep("V", inf_record$infected), which(inf_record$start_t + v_lag <= time_step  & inf_record$end_t >= time_step)),])
      }
    } else {
      lag <- RM_out$input_parameters[paste0(tolower(population), "_lag"), phase]
      dat <- inf_record[which(inf_record$start_t + lag <= time_step  & inf_record$end_t >= time_step),]
    }
  } else {
    dat <- inf_record[which(inf_record$start_t <= time_step & inf_record$end_t >= time_step),]
  }
  if (population != "both") {
    dat <- dat[grep(population, dat$infected),]
  }
  if (resample_possible == FALSE) {
    already.sampled <- inf_record[inf_record$infected=="sample", "origin_inf"]
    if (length(already.sampled) > 0) {
      dat <- dat[-which(dat$inf_id %in% already.sampled),]
    }
  }
  if (!is.null(proportion)) {
    sample.size <- as.integer(nrow(dat)*as.numeric(proportion))
  }
  if (!is.null(number)) {
    if (nrow(dat) < number & nrow(dat) > 0 & resample_possible == FALSE) {
      warning(paste0("Not enough infected individuals present to sample ", number, ". Sampling all ", nrow(dat), " infections instead."))
      sample.size <- nrow(dat)
    } else if (nrow(dat) == 0) {
      sample.size <- nrow(dat)
    } else {
      sample.size <- number
    }
  }
  if (sample.size < 1) {
    warning("Less than 1 infection chosen for sampling. Returning original RM_out object unchanged.")
    return(RM_out)
  } else {
    if (resample_possible == TRUE) {
      infs_to_sample <- dat[sample(dat$inf_id, size=sample.size, replace=T),]
    } else if (resample_possible == FALSE) {
      infs_to_sample <- dat[sample(dat$inf_id, size=sample.size, replace=F),]
    }
    sample_inf_record <- as.data.frame(cbind((nrow(inf_record)+1):(nrow(inf_record)+nrow(infs_to_sample)), infs_to_sample$inf_id, infs_to_sample$infected, rep("sample", nrow(infs_to_sample)), rep(time_step, nrow(infs_to_sample)), rep(time_step, nrow(infs_to_sample))))
    inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(infs_to_sample)),] <- sample_inf_record
    inf_record$start_t <- as.numeric(inf_record$start_t);inf_record$end_t <- as.numeric(inf_record$end_t)
    RM_out$infection_record <- inf_record
    return(RM_out)
  }
}

# 'rewinds' a Ross-Macdonald simulation by removing any events beyond the specified runtime
### TBD - need to modify input/RM parameter table accordingly
rewind.RM <- function(dat, runtime) {
  dat$indiv_status <- dat$indiv_status[,1:runtime]
  dat$infection_record <- dat$infection_record[dat$infection_record$start_t <= runtime,]
  dat$infection_record[dat$infection_record$end_t > runtime, "end_t"] <- NA
  return(dat)
}
