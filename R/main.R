#' @export
run.RM <- function(N_h,
                   N_h_t0=NULL,
                   N_v,
                   N_v_t0=NULL,
                   runtime,
                   bite_rate,
                   raw_hv_trans_rate=NULL,
                   raw_vh_trans_rate=NULL,
                   eff_hv_trans_rate=NULL,
                   eff_vh_trans_rate=NULL,
                   h_rec_rate,
                   v_rec_rate,
                   h_lag,
                   v_lag,
                   h_max_duration,
                   v_max_duration,
                   mean_hyp=0,
                   hyp_act_rate=NULL,
                   hyp_death_rate=NULL,
                   prev_sim_output=NULL,
                   verbose=TRUE) {
  `%nin%` <- Negate(`%in%`)

  # Input checks ----
  # Convert any NAs to NULLs
  for (arg in c("N_h_t0", "N_v_t0", "raw_hv_trans_rate", "raw_vh_trans_rate", "eff_hv_trans_rate", "eff_vh_trans_rate", "hyp_act_rate", "hyp_death_rate", "prev_sim_output")) {
    if (!is.null(get(arg))){if(length(get(arg))==1){if(is.na(get(arg))){assign(arg, NULL)}}}
  };rm(arg)

  # Check number of seed individuals
  if (is.null(prev_sim_output)) {
    if (is.null(N_h_t0) || is.null(N_v_t0)) {
      stop("Starting number of seed individuals must be specified.")
    }
    if (N_h < N_h_t0 || N_v < N_v_t0) {
      stop("Number of seed individuals exceeds total number of individuals.")
    }
  }

  # Check that all rates are between 0 and 1
  if (h_rec_rate < 0 | h_rec_rate > 1 | v_rec_rate < 0 | v_rec_rate > 1 ) {
    stop("Recovery rates must be between 0 and 1.")
  }

  # Calculate raw and effective transmission probabilities and check that only one of each is provided
  if (is.null(eff_hv_trans_rate) & !is.null(raw_hv_trans_rate)) {
    eff_hv_trans_rate <- calc.trans.rates(lag=h_lag,
                                          max_duration=h_max_duration,
                                          recovery_rate=h_rec_rate,
                                          raw_trans_rate=raw_hv_trans_rate)
  } else if (!is.null(eff_hv_trans_rate) & is.null(raw_hv_trans_rate)) {
    raw_hv_trans_rate <- calc.trans.rates(lag=h_lag,
                                          max_duration=h_max_duration,
                                          recovery_rate=h_rec_rate,
                                          eff_trans_rate=eff_hv_trans_rate)
  } else if (!is.null(eff_hv_trans_rate) & !is.null(raw_hv_trans_rate)) {
    stop("'raw_hv_trans_rate' and 'eff_hv_trans_rate' cannot both be specified.")
  }  else {
    stop("Either 'raw_hv_trans_rate' or 'eff_hv_trans_rate' must be specified.")
  }
  if (is.null(eff_vh_trans_rate) & !is.null(raw_vh_trans_rate)) {
    eff_vh_trans_rate <- calc.trans.rates(lag=v_lag,
                                          max_duration=v_max_duration,
                                          recovery_rate=v_rec_rate,
                                          raw_trans_rate=raw_vh_trans_rate)
  } else if (!is.null(eff_vh_trans_rate) & is.null(raw_vh_trans_rate)) {
    raw_vh_trans_rate <- calc.trans.rates(lag=v_lag,
                                          max_duration=v_max_duration,
                                          recovery_rate=v_rec_rate,
                                          eff_trans_rate=eff_vh_trans_rate)
  } else if (!is.null(eff_vh_trans_rate) & !is.null(raw_vh_trans_rate)) {
    stop("'raw_vh_trans_rate' and 'eff_vh_trans_rate' cannot both be specified.")
  }  else {
    stop("Either 'raw_vh_trans_rate' or 'eff_vh_trans_rate' must be specified.")
  }

  # Check that maximum duration is at least 1 day
  if (h_max_duration < 1 | v_max_duration < 1) {
    stop("h_max_duration and v_max_duration must be at least 1 day.")
  }

  # Check hypnozoite parameters
  if (mean_hyp < 0) {
    stop("Mean number of hypnozoites cannot be less than 0.")
  } else if (mean_hyp > 0) {
    if (is.null(hyp_act_rate) | is.null(hyp_death_rate)) {
      stop("If implementing hypnozoites, 'hyp_act_rate' and 'hyp_death_rate' must be defined.")
    }
    if (hyp_act_rate < 0 | hyp_act_rate > 1 | hyp_act_rate < 0 | hyp_act_rate > 1) {
      stop("Hypnozoite activation and death rates must be between 0 and 1.")
    }
  }
  if (mean_hyp == 0) {
    if (!is.null(hyp_act_rate)) {
      if (hyp_act_rate != -9) {
        hyp_act_rate = NULL
        message("Setting 'hyp_act_rate' to NULL because 'mean_hyp' is 0.")
      }
    }
    if (!is.null(hyp_death_rate)) {
      if (hyp_death_rate != -9) {
        hyp_death_rate = NULL
        message("Setting 'hyp_death_rate' to NULL because 'mean_hyp' is 0.")
      }
    }
  }

  # Check expected booleans
  if (!is.logical(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }

  # If previous input was provided, do some format checks
  # TO DO

  # Store parameters and initialize objects ----
  # Calculate equation parameters
  RM_parameters_temp <- calc.RM.params(N_h=N_h,
                                   N_v=N_v,
                                   bite_rate=bite_rate,
                                   eff_hv_trans_rate=eff_hv_trans_rate,
                                   eff_vh_trans_rate=eff_vh_trans_rate,
                                   h_rec_rate=h_rec_rate,
                                   v_rec_rate=v_rec_rate)
  for (param in names(RM_parameters_temp)) {
    assign(param, RM_parameters_temp[[param]])
  };rm(param, RM_parameters_temp)

  # Create vectors of host and vector names
  hIDs <- paste0(rep("H", N_h), 1:N_h)
  vIDs <- paste0(rep("V", N_v), 1:N_v)

  # Determine how many rows to add to infection record at a time
  n_infs_per_time_step <- ceiling(bite_rate * V_eq * N_v * (1-H_eq) * eff_vh_trans_rate + # infected vectors transmit to susceptible hosts
                                  bite_rate * (1 - V_eq) * N_v * H_eq * eff_hv_trans_rate) # infected hosts transmit to susceptible vectors
  if (mean_hyp > 0) {
    n_infs_per_time_step <- ceiling(n_infs_per_time_step + (1-H_eq) * H_eq * hyp_act_rate)
  }
  n_add_rows <- max(c(1e4, max(rpois(1e6, n_infs_per_time_step))));rm(n_infs_per_time_step)

  # Create data frame to store infection record
  if (is.null(prev_sim_output)) {
    infection_record <- data.table::as.data.table(matrix(ncol=6, nrow=n_add_rows))
    colnames(infection_record) <- c("inf_id","origin_inf","infector", "infected","start_t","end_t")
    col_types <- c("integer", "integer", "character", "character", "integer", "integer");names(col_types) <- colnames(infection_record)
    for (col in names(col_types)) {
      data.table::set(infection_record, j=col, value=as(infection_record[[col]], col_types[col]))
    };rm(col, col_types)
    data.table::setkey(infection_record, inf_id)
    current_row <- 1
  } else {
    infection_record <- prev_sim_output$infection_record
    current_row <- nrow(infection_record) + 1
    empty_rows <- data.table::data.table(matrix(NA, nrow = n_add_rows, ncol = ncol(infection_record)));colnames(empty_rows) <- colnames(infection_record)
    infection_record <- data.table::rbindlist(list(infection_record, empty_rows), fill = TRUE)
  }
  total_rows <- nrow(infection_record)

  # Create a vector to store individual infection age
  if (is.null(prev_sim_output)) {
    inf_age <- rep(NA, N_h + N_v)
    names(inf_age) <- c(hIDs, vIDs)
  } else {
    old_hIDs <- names(prev_sim_output$inf_age)[grep("H", names(prev_sim_output$inf_age))]
    old_vIDs <- names(prev_sim_output$inf_age)[grep("V", names(prev_sim_output$inf_age))]
    if (prev_sim_output$input_parameters["N_h", ncol(prev_sim_output$input_parameters)] == N_h & prev_sim_output$input_parameters["N_v", ncol(prev_sim_output$input_parameters)] == N_v) {
      inf_age <- prev_sim_output$inf_age
    } else {
      max_N_h = max(c(N_h, as.double(prev_sim_output$input_parameters["N_h",])))
      max_N_v = max(c(N_v, as.double(prev_sim_output$input_parameters["N_v",])))
      inf_age <- rep(NA, max_N_h+max_N_v)
      names(inf_age) <- c(paste0(rep("H", max_N_h), 1:max_N_h), paste0(rep("V", max_N_v), 1:max_N_v))
      inf_age[c(old_hIDs, old_vIDs)] <- prev_sim_output$inf_age
    }
  }

  # If specified, create objects to store information about hypnozoite reservoirs
  if (mean_hyp > 0) {
    if (is.null(prev_sim_output)) {
      hyp_reservoir <- list();length(hyp_reservoir) <- N_h;names(hyp_reservoir) <- hIDs
      n_hypno <- rep(0, N_h)
      names(n_hypno) <- hIDs
    } else if (prev_sim_output$input_parameters["mean_hyp", ncol(prev_sim_output$input_parameters)] == 0) {
      hyp_reservoir <- list();length(hyp_reservoir) <- N_h;names(hyp_reservoir) <- hIDs
      n_hypno <- rep(0, N_h)
      names(n_hypno) <- hIDs
    } else {
      if (prev_sim_output$input_parameters["N_h", ncol(prev_sim_output$input_parameters)] == N_h) {
        hyp_reservoir <- prev_sim_output$hyp_reservoir
        n_hypno <- prev_sim_output$n_hypno
      } else {
        old_hIDs <- names(prev_sim_output$inf_age)[grep("H", names(prev_sim_output$inf_age))]
        max_N_h = max(c(N_h, as.double(prev_sim_output$input_parameters["N_h",])))
        hyp_reservoir <- list();length(hyp_reservoir) <- N_h;names(hyp_reservoir) <- hIDs
        hyp_reservoir[[old_hIDs]] <- prev_sim_output$hyp_reservoir
        n_hypno <- rep(0, max_N_h)
        names(n_hypno) <- paste0(rep("H", max_N_h), 1:max_N_h)
        n_hypno[old_hIDs, 1:ncol(n_hypno)] <- prev_sim_output$n_hypno
      }
    }
  }

  # Store parameters
  if (is.null(prev_sim_output)) {
    # create data frame to store input parameters
    input_parameters <- as.data.frame(matrix(ncol=1, nrow=21))
    current_phase <- 1
    t_start <- 1
    t_end <- runtime
    rownames(input_parameters) <- c("t_start", "t_end", "runtime", "N_h", "N_h_t0", "N_v", "N_v_t0", "bite_rate", "raw_hv_trans_rate", "raw_vh_trans_rate", "eff_hv_trans_rate", "eff_vh_trans_rate", "h_rec_rate", "v_rec_rate", "h_max_duration", "v_max_duration", "h_lag", "v_lag", "mean_hyp", "hyp_death_rate", "hyp_act_rate")
    for (i in rownames(input_parameters)) {
      val <- get(i)
      if (!is.null(val)) {
        input_parameters[i, current_phase] <- val
      } else {
        input_parameters[i, current_phase] <- -9 # represents NULL
      }
    }
    # create data frame to store Ross-Macdonald parameters
    RM_parameters <- as.data.frame(matrix(ncol=1, nrow=4))
    rownames(RM_parameters) <- c("h_inf_rate", "v_inf_rate", "H_eq", "V_eq")
    for (i in rownames(RM_parameters)) {
      val <- get(i)
      if (!is.null(val)) {
        RM_parameters[i, current_phase] <- val
      }
    }
  } else {
    input_parameters <- prev_sim_output$input_parameters
    previous_phase <- ncol(input_parameters)
    current_phase <- ncol(input_parameters) + 1
    t_start <- input_parameters["t_end",  previous_phase] + 1
    t_end <- input_parameters["t_end",  previous_phase] + runtime
    # store input parameters
    for (i in rownames(input_parameters)) {
      val <- get(i)
      if (!is.null(val)) {
        input_parameters[i, current_phase] <- val
      } else {
        input_parameters[i, current_phase] <- -9 # represents NULL
      }
    }
    # if two subsequent phases have identical parameters, merge them
    param_identity_check <- 0
    for (i in 4:21) {
      if (input_parameters[i, previous_phase] == input_parameters[i, current_phase]) {
        param_identity_check <- param_identity_check + 1
      }
    }
    if (param_identity_check == 18) {
      t_end <- input_parameters["t_end", current_phase]
      new_runtime <- input_parameters["runtime", previous_phase] + input_parameters["runtime", current_phase]
      new_input_parameters <- data.frame(input_parameters[,-current_phase]);colnames(new_input_parameters) <- colnames(input_parameters)[-current_phase];rownames(new_input_parameters) <- row.names(input_parameters);input_parameters <- new_input_parameters
      current_phase <- previous_phase
      input_parameters[c("t_end", "runtime"), current_phase] <- c(t_end, new_runtime)
    }
    # store Ross-Macdonald parameters for new phase
    RM_parameters <- prev_sim_output$RM_parameters
    if (current_phase > previous_phase) {
      for (i in rownames(RM_parameters)) {
        val <- get(i)
        if (!is.null(val)) {
          RM_parameters[i, current_phase] <- val
        }
      }
    };rm(i)
  }

  # Create a data frame to store proportion of population infected over time
  time.steps <- (input_parameters["t_end", current_phase]-runtime+1):input_parameters["t_end", current_phase]
  if (is.null(prev_sim_output)) {
    proportion_infected <- as.data.frame(matrix(nrow=2, ncol=runtime))
    rownames(proportion_infected) <- c("H", "V")
  } else {
    proportion_infected <- prev_sim_output$proportion_infected
  }
  proportion_infected[,time.steps] <- NA

  ## Run simulation sequentially over each time step
  ## -----------------------------------------------
  time.steps <- (input_parameters["t_end", current_phase]-runtime+1):input_parameters["t_end", current_phase]
  for (x in time.steps) {
    if (x == 1) {
      # -- Step 0 - Initialize w/ infected hosts/vectors ----
      # choose infected individuals
      new_inf_indivs <- c(sample(hIDs, size=N_h_t0), sample(vIDs, size=N_v_t0))
      inf_age[new_inf_indivs] <- 0 # initiate infection
      proportion_infected[,x] <- 0
      # add extra rows to infection record if running low
      while ((total_rows - current_row) < length(new_inf_indivs)) {
        empty_rows <- data.table::data.table(matrix(NA, nrow = n_add_rows, ncol = ncol(infection_record)));colnames(empty_rows) <- colnames(infection_record)
        infection_record <- data.table::rbindlist(list(infection_record, empty_rows), fill = TRUE)
        total_rows <- nrow(infection_record)
      }
      # add these events to the infection record
      for (type in c("H", "V")) {
        if (length(grep(type, new_inf_indivs)) > 0) {
          inf_ids <- current_row:(current_row+length(grep(type, new_inf_indivs))-1)
          current_row <- max(inf_ids) + 1
          new_infections <- cbind.data.frame(inf_ids, 0, paste0(c("H","V")[-match(type, c("H","V"))], "-seed"), new_inf_indivs[grep(type, new_inf_indivs)], x, NA)
          infection_record[inf_ids, names(infection_record) := new_infections]
        }
      };rm(type, inf_ids)
      # create hypnozoite reservoir for newly infected hosts
      if (mean_hyp > 0) {
        new_inf_hosts <- new_inf_indivs[new_inf_indivs %in% hIDs]
        if (length(new_inf_hosts) > 0) {
          n_hyps <- rgeom(length(new_inf_hosts), prob = 1/(mean_hyp+1))
          # remove hosts that did not generate any hypnozoites
          new_inf_hosts <- new_inf_hosts[n_hyps != 0];n_hyps <- n_hyps[n_hyps != 0]
          origins_infs <- infection_record[!is.na(inf_id) & is.na(end_t)][match(infected, new_inf_hosts), "origin_inf"]
          # populate reservoir
          hyp_reservoir[new_inf_hosts] <- mapply(rep, origins_infs, n_hyps)
          n_hypno[new_inf_hosts] <- n_hypno[new_inf_hosts] - n_hyps
          rm(origins_infs, n_hyps)
        };rm(new_inf_hosts)
      };rm(new_inf_indivs)
    } else {
      # -- Step 1 - All vectors bite hosts, transmitting and/or becoming infected ----
      # increment infection age
      inf_age[!is.na(inf_age)] <- inf_age[!is.na(inf_age)] + 1
      # update proportion infected
      for (type in c("h", "v")) {
        n_infected <- sum(!is.na(inf_age[get(paste0(type, "IDs"))]))
        proportion_infected[toupper(type), x] <- n_infected/get(paste0("N_", type))
      };rm(type, n_infected)
      # create new_infections placeholder
      new_infections <- NULL
      # simulate vector biting
      n_bites <- rpois(n=N_v, lambda=bite_rate)
      biting_events <- as.data.frame(matrix(nrow=sum(n_bites), ncol=8))
      colnames(biting_events) <- c("V", "H", "V_infected", "V_infectious", "H_infected", "H_infectious", "vh_trans_event", "hv_trans_event")
      # store all bites by vectors
      biting_events$V <- unlist(mapply(rep, paste0("V", 1:N_v), n_bites))
      # are biting vectors infected or infectious?
      biting_events$V_infected <- !is.na(inf_age[biting_events$V])
      biting_events$V_infectious <- !is.na(inf_age[biting_events$V]) & inf_age[biting_events$V] > v_lag # when lag is 0, infection can be transmitted once it is 1 day old
      # sample hosts to bite
      biting_events$H <- sample(hIDs, nrow(biting_events), replace=TRUE)
      # are bitten hosts infected or infectious?
      biting_events$H_infected <- !is.na(inf_age[biting_events$H])
      biting_events$H_infectious <- !is.na(inf_age[biting_events$H]) & inf_age[biting_events$H] > h_lag
      # filter down to biting events between an infectious host and susceptible vector OR infectious vector and susceptible host
      biting_events <- biting_events[(biting_events$V_infectious & !biting_events$H_infected) | (biting_events$H_infectious & !biting_events$V_infected),]
      # if there are consequential biting events...
      if (nrow(biting_events) > 0) {
        # determine if transmission occurs
        biting_events[biting_events$V_infectious, "vh_trans_event"] <- rbinom(n=sum(biting_events$V_infectious), size=1, prob=raw_vh_trans_rate);biting_events$vh_trans_event[is.na(biting_events$vh_trans_event)] <- 0
        biting_events[biting_events$H_infectious, "hv_trans_event"] <- rbinom(n=sum(biting_events$H_infectious), size=1, prob=raw_hv_trans_rate);biting_events$hv_trans_event[is.na(biting_events$hv_trans_event)] <- 0
        transmission_events <- biting_events[na.omit(biting_events$vh_trans_event==1 | biting_events$hv_trans_event==1),]
        # if there are successful transmission events...
        if (nrow(transmission_events) > 0) {
          infectors <- c(transmission_events[transmission_events$vh_trans_event==1,"V"], transmission_events[transmission_events$hv_trans_event==1,"H"])
          active_infs <- infection_record[is.na(end_t) & !is.na(inf_id)]
          origin_infs <- active_infs$inf_id[match(infectors, active_infs$infected)]
          infecteds <- c(transmission_events[transmission_events$vh_trans_event==1,"H"], transmission_events[transmission_events$hv_trans_event==1,"V"])
          # simplify multiple infection events to single host/vector
          n_infs <- length(infecteds)
          infectors <- infectors[match(unique(infecteds), infecteds)]
          origin_infs <- origin_infs[match(unique(infecteds), infecteds)]
          infecteds <- infecteds[match(unique(infecteds), infecteds)]
          n_infs_simplified <- length(infecteds)
          n_infs_removed <- n_infs - n_infs_simplified
          if (n_infs_removed > 0) {
            if (verbose == T) {
              if (n_infs_removed>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
              message(paste0(n_infs_removed, " superinfection", msg_plural, " in time step ", x, " ", msg_were, " suppressed"))
            }
          };rm(n_infs, n_infs_simplified, n_infs_removed)
          # add extra rows to infection record if running low
          while ((total_rows - current_row) < length(infecteds)) {
            empty_rows <- data.table::data.table(matrix(NA, nrow = n_add_rows, ncol = ncol(infection_record)));colnames(empty_rows) <- colnames(infection_record)
            infection_record <- data.table::rbindlist(list(infection_record, empty_rows), fill = TRUE)
            total_rows <- nrow(infection_record)
          }
          # add these events to the infection record
          inf_ids <- current_row:(current_row+length(infecteds)-1)
          new_infections <- cbind.data.frame(inf_ids, origin_infs, infectors, infecteds, x, NA)
          infection_record[inf_ids, names(infection_record) := new_infections]
          current_row <- max(inf_ids) + 1
          # update individual statuses
          inf_age[infecteds] <- 0
          rm(inf_ids, active_infs, origin_infs, infectors, infecteds)
        };rm(transmission_events)
      };rm(n_bites, biting_events)
      # -- Step 2 - Dormant hypnozoites die/activate ----
      if (mean_hyp > 0) {
        if (sum(n_hypno) > 0) {
          # determine actions fate of each hypnozoite
          hyp_events <- as.data.frame(matrix(nrow=sum(n_hypno), ncol=5))
          colnames(hyp_events) <- c("H", "H_infected", "origin_inf", "dies", "activates")
          # which hosts have hypnozoites? what is their infection status?
          dormant_hosts <- names(n_hypno[n_hypno > 0])
          hyp_events$H <- unlist(mapply(rep, dormant_hosts, n_hypno[n_hypno > 0]))
          hyp_events$H_infected <- unlist(mapply(rep, !is.na(inf_age[dormant_hosts]), n_hypno[n_hypno > 0]))
          hyp_events$origin_inf <- as.integer(unlist(hyp_reservoir[dormant_hosts])) #####
          # simulate hypnozoite death
          hyp_events$dies <- rbinom(n=nrow(hyp_events), size=1, prob=hyp_death_rate)
          # simulate hypnozoite activation in hosts that do not have an active infection
          hyp_events$activates[!hyp_events$H_infected] <- rbinom(n=sum(!hyp_events$H_infected), size=1, prob=hyp_act_rate)
          hyp_events[is.na(hyp_events$activates), "activates"] <- 0
          hyp_events <- hyp_events[sample(1:nrow(hyp_events), size=nrow(hyp_events)),];rownames(hyp_events) <- 1:nrow(hyp_events)
          # prune events where no death or activation occurs
          hyp_events <- hyp_events[hyp_events$dies == 1 | hyp_events$activates == 1,]
          # simplify events where both death and activation occur
          hyp_events_to_simplify <- which(hyp_events$dies == 1 & hyp_events$activates == 1)
          if (length(hyp_events_to_simplify) > 0) {
            if (verbose == T) {
              if (length(hyp_events_to_simplify) > 1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
              message(paste0(length(hyp_events_to_simplify), " hypnozoite event", msg_plural, " in time step ", x, " ", msg_were, " simplified"))
            }
            hyp_events[hyp_events_to_simplify, "activates"] <- sample(c(0,1), size=length(hyp_events_to_simplify), replace=T)
            hyp_events[intersect(hyp_events_to_simplify, which(hyp_events$activates == 1)), "dies"] <- 0
          };rm(hyp_events_to_simplify)
          # simplify events where multiple hypnozoites activate
          multi_activation <- names(which(table(hyp_events[hyp_events$activates==1, "H"])>1))
          if (length(multi_activation) > 0) {
            if (verbose == T) {
              if (length(multi_activation)>1) {msg_plural="s";msg_host="each of these hosts"} else {msg_plural="";msg_host="this host"}
              message(paste0(length(multi_activation)," host", msg_plural," in time step ", x, " activated multiple hypnozoites - for ", msg_host, ", only 1 hypnozoite was actually activated"))
            }
            keep_idx <- rownames(hyp_events[hyp_events$H %in% multi_activation & hyp_events$activates == 1,])[match(multi_activation, hyp_events[hyp_events$H %in% multi_activation & hyp_events$activates == 1, "H"])]
            keep_idx <- c(keep_idx, rownames(hyp_events)[rownames(hyp_events) %nin% rownames(hyp_events[hyp_events$H %in% multi_activation & hyp_events$activates == 1,])])
            hyp_events <- hyp_events[keep_idx,]
            rm(keep_idx)
          };rm(multi_activation)
          # activate hypnozoites
          activation_events <-  hyp_events[hyp_events$activates==1, c("H", "origin_inf")]
          if (nrow(activation_events) > 0) {
            inf_ids <- current_row:(current_row+nrow(activation_events)-1)
            new_infections <- cbind.data.frame(inf_ids, activation_events$origin_inf, activation_events$H, activation_events$H, x, NA)
            # add extra rows to infection record if running low
            while ((total_rows - current_row) < length(inf_ids)) {
              empty_rows <- data.table::data.table(matrix(NA, nrow = n_add_rows, ncol = ncol(infection_record)));colnames(empty_rows) <- colnames(infection_record)
              infection_record <- data.table::rbindlist(list(infection_record, empty_rows), fill = TRUE)
              total_rows <- nrow(infection_record)
            }
            # add these events to the infection record
            infection_record[inf_ids, names(infection_record) := new_infections]
            current_row <- max(inf_ids) + 1
            hyp_reservoir[activation_events$H] <- mapply(function(x, y) x[-match(y, x)], hyp_reservoir[activation_events$H], y=activation_events$origin_inf, SIMPLIFY = F) ###
            n_hypno[activation_events$H] <- n_hypno[activation_events$H]-1
            inf_age[activation_events$H] <- 0
            rm(inf_ids)
          };rm(activation_events)
          # kill hypnozoites
          hyp_death_events <- hyp_events[hyp_events$dies==1, c("H", "origin_inf")]
          if (nrow(hyp_death_events) > 0) {
            hyp_reservoir[unique(hyp_death_events$H)] <- lapply(unique(hyp_death_events$H), function(x) {hyps <- hyp_reservoir[[x]];to_remove <- hyp_death_events[hyp_death_events$H==x, "origin_inf"];while (length(to_remove) > 0) {hyps <- hyps[-match(unique(to_remove), hyps)];to_remove <- to_remove[-match(unique(to_remove), to_remove)]};return(hyps)})
            n_hypno[names(table(hyp_death_events$H))] <- n_hypno[names(table(hyp_death_events$H))] - as.integer(table(hyp_death_events$H))
          };rm(hyp_death_events, hyp_events)
          # generate hypnozoite reservoir for hosts newly infected by a vector in Step 1 (not re-activations)
          new_inf_hosts <- new_infections[grep("H", new_infections$infecteds), "infecteds"]
          if (length(new_inf_hosts) > 0) {
            n_hyps <- rgeom(length(new_inf_hosts), prob = 1/(mean_hyp+1))
            # remove hosts that did not generate any hypnozoites
            new_inf_hosts <- new_inf_hosts[n_hyps != 0];n_hyps <- n_hyps[n_hyps != 0]
            origins_infs <- infection_record[!is.na(inf_id) & is.na(end_t)][match(infected, new_inf_hosts), "origin_inf"]
            # populate reservoir
            hyp_reservoir[new_inf_hosts] <- mapply(rep, origins_infs, n_hyps)
            n_hypno[new_inf_hosts] <- n_hypno[new_inf_hosts] - n_hyps
            rm(origins_infs, n_hyps)
          };rm(new_inf_hosts)
        }
      }
      # set all empty or NA entries in the hypnozoite reservoir to NULL
      no_hyp <- unique(c(names(which(lapply(hyp_reservoir, length) == 0)), names(which(is.na(unlist(hyp_reservoir))))))
      hyp_reservoir[no_hyp] <- NULL;names(hyp_reservoir[no_hyp]) <- no_hyp;rm(no_hyp)
      # -- Step 3 - Infected individuals clear infections ----
      no_longer_infected <- c()
      for (type in c("h", "v")) {
      	# Determine who is infected
        ids <- get(paste0(type, "IDs"))
        currently_infected <- names(na.omit(inf_age[inf_age>0 & names(inf_age) %in% ids])) # infections can only be cleared once they are at least 1 day old
        if (length(currently_infected) > 0) {
          # Individuals recover with a 'recovery rate' probability
          no_longer_infected <- c(no_longer_infected, currently_infected[which(rbinom(n=length(currently_infected), size=1, prob=get(paste0(type,"_rec_rate")))==1)])
          # Clear infections from individuals that have exceeded their maximum duration
          too_long <- names(which(inf_age[currently_infected] > get(paste0(type,"_max_duration"))))
          if (length(too_long) > 0) {
            type.full <- c("host", "vector")[match(type, c("h", "v"))]
            if (verbose == T){
              if (length(too_long)>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
              message(paste0(length(too_long), " ", type.full, " infection", msg_plural, " ", msg_were, " forcibly terminated in time step ", x))
            };rm(type.full)
            no_longer_infected <- c(no_longer_infected, too_long)
          }
          # Identify individuals that have recently emigrated
          if (!is.null(prev_sim_output) & x == time.steps[1]) {
            old_ids <- get(paste0("old_", type, "IDs"))
            if (length(old_ids) > length(ids)) {
              emigrants <- old_ids[old_ids %nin% ids]
              no_longer_infected <- c(no_longer_infected, emigrants)
              if (type == "h") {
                hyp_reservoir[[emigrants]] <- NULL
                n_hypno[emigrants] <- 0
              };rm(emigrants)
            };rm(old_ids)
          }
        }
      };rm(type, ids, currently_infected)
      if (length(no_longer_infected) > 0) {
        inf_age[no_longer_infected] <- NA
        infection_record[infected %in% no_longer_infected & is.na(end_t), end_t := x]
      };rm(no_longer_infected)

      # If no active infections remain, end simulation
      if (sum(is.na(inf_age)) == N_h + N_v) {
        if (mean_hyp == 0) {
          message(paste("Simulation ended on day", x, "- no active infections remaining."))
          break
        } else if (mean_hyp > 0 & sum(n_hypno) == 0) {
          message(paste("Simulation ended on day", x, "- no active or dormant infections remaining."))
          break
        }
      }
    }

    # Fix, store, and return output at the end of the simulation
    if (x == max(time.steps)) {
      # process infection record
      infection_record <- infection_record[!is.na(inf_id)]
      output <- list(input_parameters=input_parameters, RM_parameters=RM_parameters, infection_record=infection_record, inf_age=inf_age, proportion_infected=proportion_infected)
      if (mean_hyp > 0) {
        output <- append(output, list(hyp_reservoir=hyp_reservoir, n_hypno=n_hypno))
      }
      return(output)
    }
  }
}
