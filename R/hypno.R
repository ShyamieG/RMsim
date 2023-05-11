# generates a hypnozoite reservoir upon host infection
populate.hypno.reservoir <- function(hID,
                                     hyp_reservoir,
                                     mean_hyp,
                                     infection_record) {
  hypnozoites <- hyp_reservoir[[hID]]
  # Identify active infection
  origin_inf <- which(infection_record$infected==hID & is.na(infection_record$end_t))
  # Choose number of hypnozoites
  n_hypnozoites <- rgeom(1, prob=1/(mean_hyp+1))
  # Generate hypnozoites
  if (length(n_hypnozoites) > 0) {
    return(c(hypnozoites, rep(origin_inf, n_hypnozoites)))
  } else {
    return(hypnozoites)
  }
}

# simulates hypnozoites activating
sim.hypno.activation <- function(hID, hyp_reservoir) {
  hypnozoites <- hyp_reservoir[[hID]]
  activated_hypnozoite <- sample(1:length(hypnozoites), 1)
  # Record transmission event
  infections <- c(hypnozoites[activated_hypnozoite], hID, hID, x, NA)
  # Remove hypnozoite from reservoir
  hypnozoites <- hypnozoites[-activated_hypnozoite]
  return(list(new.infections=infections, hypnozoites=hypnozoites))
}

# simulates hypnozoites dying
sim.hypno.death <- function(hID, n_to_kill, hyp_reservoir) {
  hypnozoites <- hyp_reservoir[[hID]]
  if (length(hypnozoites)>0 & n_to_kill > 0) {
    # Remove hypnozoites from the reservoir
    hypnozoites <- hypnozoites[-sample(1:length(hypnozoites), size=n_to_kill)]
  }
  return(hypnozoites)
}
