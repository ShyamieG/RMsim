#' @title Sample infected individuals from a Ross-Macdonald simulation
#' @description Adds events to the infection record representing samples taken from infected individuals, according to user-specified parameters.
#' @usage sample.RM(RM_out, time_step, population, proportion = NULL, number = NULL, sample_post_lag = TRUE, resample_possible = TRUE, sort_events = TRUE)
#' @details
#' The sampling of infected individuals is represented by a unique type of event that is included in the infection record of a Ross-Macdonald simulation. This function allows users to define the sampling scheme by specifying the time step when infected individuals should be sampled, which population (host, vector, or both) to draw samples from, whether a certain proportion of the infected population or a fixed number of individuals should be sampled, and if samples should only be drawn from individuals that are beyond the 'lag' period of their infection. Multiple samples can be drawn from the same infectious individual if `resample_possible` is set to TRUE. If `sort_events` is set to TRUE (which it is by default) the updated infection record will be sorted by `start_t` (infection start time or sample time).
#' @param RM_out object containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param time_step the time step at which infected samples should be drawn
#' @param population either 'H' (host) or 'V' (vector); the population to draw samples from
#' @param proportion if specified, the proportion of the total infected population to sample
#' @param number if specified, the number of infected individuals to sample
#' @param sample_post_lag if TRUE, only individuals that are beyond the infection lag phase can be sampled. Default is TRUE.
#' @param resample_possible if TRUE, infected individuals can be sampled multiple times at the same or different time points. Default is FALSE.
#' @param sort_events if TRUE, the updated infection record with new sampling events will be sorted by `start_t`. Default is TRUE.
#' @examples
#' ## add samples to the basic Ross-Macdonald simulation example
#'  set.seed(1234)
#'  RM_sampled <- sample.RM(RM_out = RMsim::RM_out_basic_sim, time_step = 1000, population = "H", number = 5, sort_events = FALSE)
#'  tail(RM_sampled$infection_record)
#' @export
sample.RM <- function(RM_out,
                      time_step,
                      population,
                      proportion = NULL,
                      number = NULL,
                      sample_post_lag = TRUE,
                      resample_possible = FALSE,
                      sort_events = TRUE) {
  # Check values
  `%ni%` <- Negate(`%in%`)
  if (is.null(population)) {
    stop("\'population\' must be either \'H\' (host), \'V\' (vector), or \'both\'")
  } else {
    population <- unique(tolower(population))
    if (any(population != "both" & population %ni% c("host", "vector") & population %ni% c("h", "v"))) {
      stop("\'population\' must be either \'H\' (host), \'V\' (vector), or \'both\'")
    }
  }

  if (is.null(proportion) & is.null(number)) {
    stop("Either the 'proportion' of the infected population or the 'number' of samples must be specified.")
  }  else if (!is.null(proportion) & !is.null(number)) {
    stop("'proportion' and 'number' cannot both be specified.")
  }

  input_parameters <- RM_out$input_parameters
  infection_record <- RM_out$infection_record

  if (input_parameters["t_end", ncol(input_parameters)] < time_step) {
    stop("Sample time ('time_step') must fall within the time course of the simulation.")
  }

  phase <- which(input_parameters["t_start", ] <= time_step & input_parameters["t_end", ] >= time_step)

  # Re-name population elements if necessary
  if (any(population == "both")){population <- c("h", "v")}
  if (any(population %in% c("host", "vector"))) {population[population=="host"] <- "h";population[population=="vector"] <- "v"}

  if (sample_post_lag) {
    dat <- NULL
    for (p in population) {
      lag <- input_parameters[paste0(p, "_lag"), phase]
      ids <- paste0(toupper(p), 1:input_parameters[paste0("N_", p), phase])
      dat <- rbind(dat, infection_record[infected %in% ids & start_t + lag < time_step & (end_t >= time_step | is.na(end_t))])
    }
  } else {
    dat <- infection_record[start_t < time_step & (end_t >= time_step | is.na(end_t))]
    if (length(population) == 1) {
      dat <- dat[grep(toupper(population), infected),]
    }
  }

  if (resample_possible == FALSE) {
    already.sampled <- infection_record[infected=="sample", "origin_inf"]
    if (length(already.sampled) > 0) {
      dat <- dat[!(dat$inf_id %in% already.sampled),]
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
    if (resample_possible) {
      infs_to_sample <- dat[inf_id %in% sample(dat$inf_id, size=sample.size, replace=T)]
    } else {
      infs_to_sample <- dat[inf_id %in% sample(dat$inf_id, size=sample.size, replace=F)]
    }
    inf_ids <- (nrow(infection_record)+1):(nrow(infection_record)+nrow(infs_to_sample))
    sample_infection_record <- cbind.data.frame(inf_ids, infs_to_sample$inf_id, infs_to_sample$infected, "sample", time_step, time_step);colnames(sample_infection_record) <- colnames(infection_record)
    infection_record <- data.table::rbindlist(list(infection_record, sample_infection_record))
    if (sort_events) {
      infection_record <- infection_record[order(start_t)]
    }
    RM_out$infection_record <- infection_record
    return(RM_out)
  }
}

#' @export
# remove any infections that did not contribute to the sample
prune.infection.record <- function(infection_record) {
  `%ni%` <- Negate(`%in%`)
  # Determine sample infections
  sample_infections <- infection_record[infected=="sample"]$inf_id
  if (length(sample_infections) == 0) {
    stop("There are no samples in this infection record. Did you run sample.RM() first?")
  }
  # Determine seed infections
  seed_infs <- infection_record[grep("seed", infector)]$inf_id
  keep_infs <- c()
  # Loop over each sample infection and trace back until a seed infection is reached
  for (i in sample_infections) {
    j_inf_history <- c()
    j <- i
    while (j %ni% c(seed_infs, keep_infs)) {
      # update j
      j <- infection_record[inf_id==j, "origin_inf"][[1]]
      # keep this infection in pruned record
      j_inf_history[length(j_inf_history)+1] <- j
    }
    keep_infs[(length(keep_infs)+1):(length(keep_infs)+length(j_inf_history))] <- j_inf_history
  }
  # Keep all focal transmissions
  keep_infs[(length(keep_infs)+1):(length(keep_infs)+length(sample_infections))] <- sample_infections
  # Remove duplicates
  pruned_infection_record <- infection_record[inf_id %in% unique(keep_infs)]
  return(pruned_infection_record)
}

#' @export
# convert pruned output to format needed by SLiM-based pipeline
generate.SLiM.input <- function(pruned_infection_record, sim_params) {
  output <- data.table::data.table(t(sapply(X=unique(pruned_infection_record$origin_inf),
                                            FUN=populate.SLiM.table,
                                            pruned_infection_record=pruned_infection_record)))
  colnames(output) <- c("infection_idx", "infection_source", "infector_type", "infector_id", "infected_ids", "RM_time_start", "inf_time_start")
  # fill in infection start times
  for (infector_type in c("h", "v")) {
    seed_infs <- intersect(grep("NA", output$inf_time_start), which(output$infector_type==toupper(infector_type)))
    if (length(seed_infs) > 0) {
      infected_type <- c("v", "h")[-match(infector_type, c("v", "h"))]
      trans_prob_distr <- (1-sim_params[paste0(infector_type, "_rec_rate"), 1])^(1:(sim_params[paste0(infector_type, "_max_duration"), 1])) * sim_params["bite_rate", 1] * sim_params[paste0("raw_", infector_type, infected_type, "_trans_rate"), 1]
      trans_prob_distr[1:(sim_params[paste0(infector_type, "_lag"), 1]-1)] <- 0
      trans_prob_distr <- as.data.frame(cbind(1:sim_params[paste0(infector_type, "_max_duration"), 1], trans_prob_distr));colnames(trans_prob_distr) <- c("x", "prob")
      trans_prob_distr <- pdqr::new_r(trans_prob_distr, type="discrete")
      for (i in seed_infs) {
        start_times <- unlist(strsplit(unlist(output[i, "inf_time_start"])[[1]], split=";"))
        output[i, inf_time_start := paste(trans_prob_distr(length(start_times)), collapse=";")]
        output[i, infection_idx := paste0(infector_type, '-seed')]
      }
    }
  }
  return(output)
}

# sub-function of generate.SLiM.input()
populate.SLiM.table <- function(pruned_infection_record, inf_id) {
  inf_id0 <- inf_id # to distinguish from column name
  # Filter to infections that resulted from inf_id
  dat <- pruned_infection_record[origin_inf==inf_id0]
  # Determine infector ID
  infector <- unique(dat$infector)
  # Determine who was infected over the course of this infection
  infected <- dat$infected
  # Determine the infector type
  infector_type <- substring(infector,1,1)
  # Determine when these infections began in Ross-Macdonald sim time
  RM_start_times <- dat$start_t
  # If not a seed infection...
  if (inf_id0 != 0) {
    # Determine origin of inf_id
    origin <- pruned_infection_record[inf_id==inf_id0,"origin_inf"]
    # Determine infector id
    infector_id <- substring(infector, 2)
    # Determine when the infector was infected in RM time
    infector_start_time <- as.integer(pruned_infection_record[inf_id==inf_id0, "start_t"])
    # If a seed infection...
  } else {
    origin <- NA
    infector_id <- "seed"
    infector_start_time <- NA
  }
  # Determine when the infections began wrt to infector's status
  start_times <- RM_start_times - infector_start_time
  # Return output
  output <- cbind.data.frame(inf_id0,
                             origin,
                             infector_type,
                             infector_id,
                             paste(infected, collapse=";"),
                             paste(RM_start_times, collapse=";"),
                             paste(start_times, collapse=";"))
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
