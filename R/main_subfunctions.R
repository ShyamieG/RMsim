# simulates vectors biting hosts
sim.vector.biting <- function(v, hIDs, t,
                              bite_rate,
                              v_lag,
                              h_lag,
                              hv_trans_rate,
                              vh_trans_rate,
                              indiv_status,
                              inf_record) {
  # Determine vector's infection status
  v_inf_status <- indiv_status[v, paste0("T",t-1)]
  # Determine number of bites
  n_bites <- rpois(n=1, lambda=bite_rate)
  if (n_bites > 0) {
    # Choose hosts to bite
    hosts_bitten <- sample(x=hIDs, size=n_bites, replace=T)
    # Bite hosts and transmit or become infected
    dat_v <- lapply(X=hosts_bitten, FUN=sim.transmission,
                    v=v, t=t,
                    v_inf_status=v_inf_status,
                    v_lag=v_lag,
                    h_lag=h_lag,
                    hv_trans_rate=hv_trans_rate,
                    vh_trans_rate=vh_trans_rate,
                    indiv_status=indiv_status,
                    inf_record=inf_record)
    names(dat_v) <- hosts_bitten
    # Update individual status
    indiv_status[,paste0("T",t)] <- summarize.status(dat_v, t)
    # Record new infections
    new_infections <- as.data.frame(do.call(rbind,lapply(dat_v, function(l) l[[2]])))
    if (ncol(new_infections)==0) {
      new_infections <- as.data.frame(matrix(ncol=5, nrow=0))
    }
  } else {
    new_infections <- as.data.frame(matrix(ncol=5, nrow=0))
  }
  colnames(new_infections) <- c("origin_inf","infector", "infected","start_t","end_t")
  return(list(indiv_status=indiv_status, new_infections=new_infections))
}

# simulates transmission from host to vector and vice versa
sim.transmission <- function(h, v, t,
                             v_inf_status,
                             v_lag,
                             h_lag,
                             hv_trans_rate,
                             vh_trans_rate,
                             indiv_status,
                             inf_record) {
  infections <- NULL
  infection_occurs <- 0
  # Check host infection status
  h_inf_status <- indiv_status[h, paste0("T",t-1)]
  # If vector is infected and host is susceptible...
  if (v_inf_status == 1 & h_inf_status == 0) {
    # Does the host become infected?
    inf_age <- t - as.numeric(inf_record[inf_record$infected==v & is.na(inf_record$end_t), "start_t"])
    if (inf_age > v_lag) {
      infection_occurs <- rbinom(n=1, size=1, prob=vh_trans_rate)
    }
    if (infection_occurs == 1) {
      # Update individual infection status
      indiv_status[h, paste0("T",t)] <- 1
      # Determine index/source of the vector's current infection
      origin_inf <- rownames(inf_record[inf_record$infected==v & is.na(inf_record$end_t),])
      # Record transmission event
      infections <- c(origin_inf, v, h, t, NA)
    }
  }
  # If host is infected and vector is susceptible...
  if (v_inf_status == 0 & h_inf_status == 1) {
    # Does the vector become infected?
    inf_age <- t - as.numeric(inf_record[inf_record$infected==h & is.na(inf_record$end_t), "start_t"])+1
    if (inf_age > h_lag) {
      infection_occurs <- rbinom(n=1, size=1, prob=hv_trans_rate)
    }
    if (infection_occurs == 1) {
      # Update individual status
      indiv_status[v, paste0("T",t)] <- 1
      # Determine index of the host's current infection
      origin_inf <- rownames(inf_record[inf_record$infected==h & is.na(inf_record$end_t),])
      # Record transmission event
      infections <- c(origin_inf, h, v, t, NA)
    }
  }
  return(list(indiv_status=indiv_status, new.infections=infections))
}

# summarize results of vector biting and transmission events on individual status
summarize.status <- function(dat, t) {
  dat0 <- lapply(dat, function(l) l[[1]][,paste0("T",t)])
  new_status <- apply(X=do.call(cbind,dat0), FUN=sum, MARGIN=1)
  if (any(new_status>1)) {
    new_status[new_status>1] <- 1
  }
  return(new_status)
}

# generates a hypnozoite reservoir upon host infection
populate.hypno.reservoir <- function(hID,
                                     mean_hyp,
                                     hyp_distr,
                                     inf_record) {
  # Identify active infection
  origin_inf <- which(inf_record$infected==hID & is.na(inf_record$end_t))
  # Choose number of hypnozoites
  n_hypnozoites <- rpois(1, lambda=mean_hyp)
  # Generate hypnozoites
  if (length(n_hypnozoites) > 0) {
    output <- rep(origin_inf, n_hypnozoites)
  } else {
    output <- NULL
  }
  return(output)
}

# simulates hypnozoites activating
sim.hypno.activation <- function(h, t,
                                 hyp_reservoir,
                                 hyp_act_rate,
                                 indiv_status) {
  infections <- NULL
  # Decide which hypnozoites activate
  hypnozoites <- hyp_reservoir[[h]]
  activated_hypnozoite <- hypnozoites[rbinom(length(hypnozoites), 1, hyp_act_rate) == 1]
  # If more than one activates, choose just one
  if (length(activated_hypnozoite) > 1) {
    message(paste0("WARNING: ", length(activated_hypnozoite)," hypnozoites were activated in time step ", t, " - only 1 was kept"))
    activated_hypnozoite <- sample(activated_hypnozoite, size=1)
  }
  # Activate hypnozoites
  if (length(activated_hypnozoite) == 1) {
    # Remove hypnozoite from reservoir
    hypnozoites <- hypnozoites[-match(activated_hypnozoite, hypnozoites)]
    # Update host infection status
    indiv_status[h, paste0("T", t)] <- 1
    # Record transmission event
    infections <- c(activated_hypnozoite, h, h, t, NA)
  }
  return(list(indiv_status=indiv_status, new.infections=infections, hypnozoites=hypnozoites))
}

# simulates hypnozoites dying
sim.hypno.death <- function(h,
                            hyp_reservoir,
                            hyp_death_rate) {
  # Determine while hypnozoites die
  hypnozoites <- hyp_reservoir[[h]]
  if (length(hypnozoites)>0) {
    hyp_deaths <- rbinom(n=length(hypnozoites), size=1, prob=hyp_death_rate)
    # Remove these hypnozoites from the reservoir
    if (sum(hyp_deaths) > 0) {
      hypnozoites <- hypnozoites[-which(hyp_deaths==1)]
    }
  }
  return(hypnozoites)
}
