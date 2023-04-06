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
                   sampling_df=NULL,
                   prev_sim_output=NULL) {
  `%nin%` <- Negate(`%in%`)

  # Check number of seed individuals
  if (is.null(prev_sim_output)) {
    if (is.null(N_h_t0) || is.null(N_v_t0)) {
      stop("Starting number of seed individuals must be specified.")
    }
    if (N_h < N_h_t0 || N_v < N_v_t0) {
      stop("Number of seed individuals exceeds total number of individuals.")
    }
  }

  # Check hypnozoite parameters
  if (mean_hyp < 0) {
    stop("Mean number of hypnozoites cannot be less than 0")
  } else if (mean_hyp > 0) {
    if (is.null(hyp_act_rate) | is.null(hyp_death_rate)) {
      stop("If implementing hypnozoites, 'hyp_act_rate' and 'hyp_death_rate' must be defined.")
    }
    if (hyp_act_rate < 0 | hyp_act_rate > 1 | hyp_act_rate < 0 | hyp_act_rate > 1) {
      stop("Hypnozoite activation and death rates must be between 0 and 1.")
    }
  }

  # Check that all rates are between 0 and 1
  if (h_rec_rate < 0 | h_rec_rate > 1 | v_rec_rate < 0 | v_rec_rate > 1 ) {
    stop("Recovery rates must be between 0 and 1.")
  }

  # Calculate raw and effective transmission probabilities
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

  # Calculate equation parameters
  RM_params_temp <- calc.RM.params(N_h=N_h,
                                   N_v=N_v,
                                   bite_rate=bite_rate,
                                   eff_hv_trans_rate=eff_hv_trans_rate,
                                   eff_vh_trans_rate=eff_vh_trans_rate,
                                   h_rec_rate=h_rec_rate,
                                   v_rec_rate=v_rec_rate)
  for (param in names(RM_params_temp)) {
    assign(param, RM_params_temp[[param]])
  };rm(RM_params_temp)

  # Create vectors of host and vector names
  hIDs <- paste0(rep("H", N_h), 1:N_h)
  vIDs <- paste0(rep("V", N_v), 1:N_v)

  # Create dataframe to store infection record
  if (is.null(prev_sim_output)) {
    inf_record <- as.data.frame(matrix(ncol=6, nrow=0))
    colnames(inf_record) <- c("inf_id","origin_inf","infector", "infected","start_t","end_t")
  } else {
    inf_record <- prev_sim_output$infection_record
  }

  # Create matrix to store individual infection states over time
  if (is.null(prev_sim_output)) {
    indiv_status <- as.data.frame(matrix(0, nrow=N_h+N_v, ncol=0))
    rownames(indiv_status) <- c(hIDs, vIDs)
  } else {
    if (prev_sim_output$input_parameters["N_h", ncol(prev_sim_output$input_parameters)] == N_h & prev_sim_output$input_parameters["N_v", ncol(prev_sim_output$input_parameters)] == N_v) {
      indiv_status <- prev_sim_output$indiv_status
    } else {
      old_hIDs <- rownames(prev_sim_output$indiv_status)[grep("H", rownames(prev_sim_output$indiv_status))]
      old_vIDs <- rownames(prev_sim_output$indiv_status)[grep("V", rownames(prev_sim_output$indiv_status))]
      max_N_h = max(c(N_h, as.double(prev_sim_output$input_parameters["N_h",])))
      max_N_v = max(c(N_v, as.double(prev_sim_output$input_parameters["N_v",])))
      indiv_status <- as.data.frame(matrix(0, nrow=max_N_h+max_N_v, ncol=ncol(prev_sim_output$indiv_status)))
      rownames(indiv_status) <- c(paste0(rep("H", max_N_h), 1:max_N_h), paste0(rep("V", max_N_v), 1:max_N_v))
      colnames(indiv_status) <- paste0("T", 1:ncol(indiv_status))
      indiv_status[c(old_hIDs, old_vIDs),] <- prev_sim_output$indiv_status
    }
  }

  # If specified, create objects to store information about hypnozoite reservoirs
  if (mean_hyp > 0) {
    if (is.null(prev_sim_output) || prev_sim_output$input_parameters["mean_hyp", ncol(prev_sim_output$input_parameters)] == 0) {
      hyp_reservoir <- list();length(hyp_reservoir) <- N_h;names(hyp_reservoir) <- hIDs
      n_hypno <- as.data.frame(matrix(0, nrow=N_h, ncol=runtime));rownames(n_hypno) <- hIDs;colnames(n_hypno) <- paste0("T",1:runtime)
    } else {
      hyp_reservoir <- prev_sim_output$hyp_reservoir
      n_hypno <- prev_sim_output$n_hypno
    }
  }

  # Store parameters
  if (is.null(prev_sim_output)) {
    # store input parameters
    input_params <- as.data.frame(matrix(ncol=1, nrow=19))
    rownames(input_params) <- c("runtime", "N_h", "N_h_t0", "N_v", "N_v_t0", "bite_rate", "raw_hv_trans_rate", "raw_vh_trans_rate", "eff_hv_trans_rate", "eff_vh_trans_rate", "h_rec_rate", "v_rec_rate", "h_max_duration", "v_max_duration", "h_lag", "v_lag", "mean_hyp", "hyp_death_rate", "hyp_act_rate")
    for (i in rownames(input_params)) {
      if (!is.null(get(i))) {
        input_params[i, ncol(input_params)] <- get(i)
      }
    }
    # store Ross-Macdonald parameters
    RM_params <- as.data.frame(matrix(ncol=1, nrow=6))
    rownames(RM_params) <- c("t_start", "t_end", "h_inf_rate", "v_inf_rate", "H_eq", "V_eq")
    RM_params[c("t_start", "t_end"), ncol(RM_params)] <- c(1, runtime)
    for (i in rownames(RM_params)[-(1:2)]) {
      if (!is.null(get(i))) {
        RM_params[i, ncol(RM_params)] <- get(i)
      }
    }
  } else {
    input_params <- prev_sim_output$input_parameters
    # store input parameters
    new_col = ncol(input_params)+1
    for (i in rownames(input_params)) {
      if (!is.null(get(i))) {
        input_params[i, new_col] <- get(i)
      }
    }
    # store Ross-Macdonald parameters
    RM_params <- prev_sim_output$RM_parameters
    RM_params[c("t_start", "t_end"), new_col] <- c(ncol(indiv_status)+1, ncol(indiv_status)+runtime)
    for (i in rownames(RM_params)[-(1:2)]) {
      if (!is.null(get(i))) {
        RM_params[i, ncol(RM_params)] <- get(i)
      }
    }
  }

  ## Run simulation sequentially over each time step
  ## -----------------------------------------------
  if (is.null(prev_sim_output)) {
    time.steps <- 1:runtime
  } else {
    time.steps <- (ncol(indiv_status)+1):(ncol(indiv_status)+runtime)
  }

  for (t in time.steps) {
    if (t == 1) {
      # -- Step 0 - Initialize w/ infected hosts/vectors
      indiv_status[,t] <- 0 # new time step (adds column to indiv_status)
      colnames(indiv_status)[t] <- paste0("T",t)
      # Choose infected individuals
      new_inf_indivs <- c(sample(hIDs, size=N_h_t0), sample(vIDs, size=N_v_t0))
      indiv_status[new_inf_indivs,t] <- 1 # infect these
      for (type in c("H", "V")) {
        if (length(grep(type, new_inf_indivs)) > 0) {
          entry_rows <- (nrow(inf_record)+1):(nrow(inf_record)+length(grep(type, new_inf_indivs)))
          inf_record[entry_rows,] <- cbind(entry_rows, paste0(c("H","V")[-match(type, c("H","V"))], "-seed"), paste0(c("H","V")[-match(type, c("H","V"))], "-seed"), new_inf_indivs[grep(type, new_inf_indivs)], t, NA)
        }
      }
      # Create hypnozoite reservoir for newly infected hosts
      if (mean_hyp > 0) {
        new_inf_hosts <- new_inf_indivs[new_inf_indivs %in% hIDs]
        if (length(new_inf_hosts) > 0) {
          hyp_reservoir[new_inf_hosts] <- lapply(X=new_inf_hosts, populate.hypno.reservoir, mean_hyp=mean_hyp, inf_record=inf_record)
        }
        # Record total number of hypnozoites per host
        n_hypno[, "T1"] <- as.numeric(unlist(lapply(hyp_reservoir, length)))
      }
    } else {
      # -- Step 1 - All vectors bite hosts, transmitting and/or becoming infected
      indiv_status[,t] <- indiv_status[,t-1] # new time step (adds column to indiv_status)
      colnames(indiv_status)[t] <- paste0("T",t)
      # Simulate vector biting
      dat <- lapply(X=vIDs, FUN=sim.vector.biting,
                    hIDs=hIDs, t=t,
                    bite_rate=bite_rate,
                    vh_trans_rate=raw_vh_trans_rate,
                    hv_trans_rate=raw_hv_trans_rate,
                    h_lag=h_lag,
                    v_lag=v_lag,
                    indiv_status=indiv_status[,paste0("T",(t-1):t)],
                    inf_record=inf_record)
      names(dat) <- vIDs
      # Update individual statuses
      indiv_status[,t] <- summarize.status(dat, t)
      # Record infections
      new_infections <- as.data.frame(do.call(rbind, lapply(dat, function(l) l[[2]])))
      if (nrow(new_infections) > 0) {
        # If an individual is infected more than once in this time step, only keep one
        co_infected <- which(table(new_infections$infected) > 1)
        if (length(co_infected) > 0) {
          co_infected <- names(table(new_infections$infected))[co_infected]
          if (length(co_infected)>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
          message(paste0("WARNING: ", length(co_infected), " co-infection", msg_plural, " in time step ", t, " ", msg_were, " simplified"))
          new_infections <- rbind(new_infections[-which(new_infections$infected %in% co_infected),], new_infections[match(co_infected, new_infections$infected),])
        }
        rownames(new_infections) <- (nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections))
        inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)),] <- cbind((nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)), new_infections)
      }
      # Create hypnozoite reservoir for hosts newly infected by a vector (not re-activations)
      if (mean_hyp > 0) {
        new_inf_hosts <- hIDs[indiv_status[hIDs,t] == 1 & indiv_status[hIDs,t-1] == 0]
        new_inf_hosts <- inf_record[intersect(which(inf_record$infected %in% new_inf_hosts & is.na(inf_record$end_t)), grep("V", inf_record$infector)), "infected"]
        if (length(new_inf_hosts) > 0) {
          new_hyp_reservoir <- lapply(X=new_inf_hosts, populate.hypno.reservoir, mean_hyp=mean_hyp, inf_record=inf_record)
          names(new_hyp_reservoir) <- new_inf_hosts
        }
      }
      # -- Step 2 - Dormant hypnozoites die/activate
      if (mean_hyp > 0) {
        # Simulate hypnozoite death
        hyp_reservoir[hIDs] <- lapply(X=hIDs, sim.hypno.death, hyp_reservoir=hyp_reservoir, hyp_death_rate=hyp_death_rate)
        # Identify hosts with hypnozoites but no active infection
        dormant_hosts <- rownames(indiv_status[which(indiv_status[hIDs,t]==0 & lapply(hyp_reservoir, FUN=length)>0),])
        if (length(dormant_hosts) > 0) {
          # Simulate hypnozoite reactivation
          dat <- lapply(X=dormant_hosts, t=t,
                        FUN=sim.hypno.activation,
                        hyp_reservoir=hyp_reservoir,
                        hyp_act_rate=hyp_act_rate,
                        indiv_status=indiv_status[,paste0("T",(t-1):t)])
          names(dat) <- dormant_hosts
          # Update individual statuses
          indiv_status[,t] <- summarize.status(dat, t)
          # Record infections
          new_infections <- as.data.frame(do.call(rbind, lapply(dat, function(l) l[[2]])))
          if (nrow(new_infections) > 0) {
            rownames(new_infections) <- (nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections))
            inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)),] <- cbind((nrow(inf_record)+1):(nrow(inf_record)+nrow(new_infections)), new_infections)
          }
          # Update hypnozoite reservoir
          hyp_reservoir[dormant_hosts] <- lapply(dat, function(l) l[[3]])
        }
        # Combine new_hyp_reservoir with hyp_reservoir
        if (length(new_inf_hosts) > 0) {
          hyp_reservoir <- setNames(mapply(c, hyp_reservoir[hIDs], new_hyp_reservoir[hIDs]), hIDs)
        }
        # Record total number of hypnozoites per host
        n_hypno[, paste0("T",t)] <- as.numeric(unlist(lapply(hyp_reservoir, length)))
      }

      # -- Step 3 - Infected individuals clear infections
      for (type in c("h", "v")) {
        # Determine who is infected
        ids <- get(paste0(type, "IDs"))
        inf_indivs <- ids[which(indiv_status[ids,t-1]==1)]
        if (length(inf_indivs) > 0) {
          # Clear infections with 'recovery rate' probability
          indiv_status[inf_indivs,t] <- indiv_status[inf_indivs,t-1] - rbinom(n=length(inf_indivs), size=1, prob=get(paste0(type,"_rec_rate")))
          # Clear infections from individuals that have exceeded their maximum duration
          active_infs <- inf_record[inf_record$infected %in% inf_indivs & is.na(inf_record$end_t),]
          too_long <- active_infs[t - as.numeric(active_infs$start_t) == get(paste0(type,"_max_duration"))-1, "infected"]
          if (length(too_long) > 0) {
            type.full <- c("host", "vector")[match(type, c("h", "v"))]
            if (length(too_long)>1) {msg_plural="s";msg_were="were"} else {msg_plural="";msg_were="was"}
            message(paste0("WARNING: ", length(too_long), " ", type.full, " infection", msg_plural, " ", msg_were, " forcibly terminated in time step ", t))
            indiv_status[too_long,t] <- 0
          }
          # Clear infections from individuals that have recently emigrated
          if (!is.null(prev_sim_output) & t == time.steps[1]) {
            if (input_params["N_h",ncol(input_params)-1] > input_params["N_h",ncol(input_params)]) {
              emigrants <- old_hIDs[old_hIDs %nin% hIDs]
              indiv_status[emigrants, t] <- 0
              inf_record[inf_record$infected %in% emigrants & is.na(inf_record$end_t), "end_t"] <- t
            }
            if (input_params["N_v",ncol(input_params)-1] > input_params["N_v",ncol(input_params)]) {
              emigrants <- old_vIDs[old_vIDs %nin% vIDs]
              indiv_status[emigrants, t] <- 0
              inf_record[inf_record$infected %in% emigrants & is.na(inf_record$end_t), "end_t"] <- t
            }
          }
          # Update infection record
          cleared <- inf_indivs[indiv_status[inf_indivs,t]==0]
          inf_record[inf_record$infected %in% cleared & is.na(inf_record$end_t),"end_t"] <- t
        }
      }
      # If no active infections remain, end simulation
      if (sum(indiv_status[,t])==0) {
        if (mean_hyp==0) {
          message(paste("Simulation ended on day", t, "- no active infections remaining."))
          break
        } else if (mean_hyp>0 & sum(n_hypno[,t])==0) {
          message(paste("Simulation ended on day", t, "- no active or dormant infections remaining."))
        }
      }

      # Sample infected individuals, if specified
      if (!is.null(sampling_df)) {
        if (t %in% sampling_df$t) {
          sampling_t <- sampling_df[sampling_df$t==t,]
          # Extract individuals that are currently infected
          dat <- inf_record[is.na(inf_record$end_t), ]
          for (i in 1:nrow(sampling_t)) {
            matching_rows <- grep(sampling_t[i, "pop"], dat$infected)
            samples <- dat[sample(matching_rows, size=length(matching_rows)*as.numeric(sampling_t$proportion[i])),]
            if (nrow(samples) > 0) {
              dat.new <- as.data.frame(cbind((nrow(inf_record)+1):(nrow(inf_record)+nrow(samples)), samples$inf_id, samples$infected, rep("sample", nrow(samples)), rep(t, nrow(samples)), rep(t, nrow(samples))))
              inf_record[(nrow(inf_record)+1):(nrow(inf_record)+nrow(samples)),] <- dat.new
            }
          }
        }
      }
    }

    # Fix, store, and return output at the end of the simulation
    if (t == max(time.steps)) {
      inf_record$start_t <- as.numeric(inf_record$start_t)
      inf_record$end_t <- as.numeric(inf_record$end_t)
      output <- list(input_parameters=input_params, RM_parameters=RM_params, indiv_status=indiv_status, infection_record=inf_record)
      if (mean_hyp > 0) {
        output <- append(output, list(hyp_reservoir=hyp_reservoir, n_hypno=n_hypno))
      }
      if (!is.null(sampling_df)) {
        output <- append(output, list(sampling=sampling_df))
      }
      return(output)
    }
  }
}
