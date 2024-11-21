#' @title Calculate raw and effective transmission rates
#' @description Converts 'raw' transmission rates to 'effective' transmission rates and vice versa.
#' @usage calc.trans.rates(lag, max_duration, recovery_rate, raw_trans_rate = NULL, eff_trans_rate = NULL)
#' @details
#' Because RMsim allows for a lag between infection and infectiousness, as well as an upper limit on the duration of infection, we distinguish between 'raw' transmission rates and 'effective' transmission rates. This function converts between these two types of transmission rate.\cr\cr
#' The 'raw' transmission rate is the probability that an infectious bite leads to a transmission event. The 'raw' transmission rate is used by the main Ross-Macdonald simulation function `run.RM()`.\cr\cr
#' The 'effective' transmission rate represents the overall probability that an infection will be transmitted given that transmission cannot occur during the `lag` period or after the `max_duration` period. The 'effective' transmission rate is used to calculate the expected Ross-Macdonald parameters, including the expected proportion of infected hosts and vectors at equilibrium.\cr\cr
#' Either 'raw_trans_rate' or 'eff_trans_rate' must be specified, but not both.
#' @param lag number of days it takes for an individual to become contagious post infection
#' @param max_duration maximum number of days that an infection can last
#' @param recovery_rate the per day probability that an infected host will recover (or die and be replaced: i.e. transition to the uninfected and susceptible state)
#' @param raw_trans_rate probability that an infectious bite leads to a transmission event (used in the main simulation function `run.RM()`)
#' @param eff_trans_rate the effective transmisson rate, accounting for the initial lag period when the individual is not infectious (`lag`) and the maximum duration of any infection (`max_duration`) (used for Ross-Macdonald parameter calculations)
#' @returns The corresponding 'raw' or 'effective' transmission rate (float).
#' @examples
#' ## calculate the 'effective' transmission rate, given a 'raw' rate of 0.2, infection lag time of 7 days, maximum infection duration of 60 days, and a per-day recovery rate of 0.05
#'  calc.trans.rates(lag = 7, max_duration = 60, recovery_rate = 0.05, raw_trans_rate = 0.2)
#'  # returns a value of 0.1302141
#'
#' ## calculate the 'raw' transmission rate, given an 'effective' rate of 0.1, infection lag time of 14 days, maximum infection duration of 180 days, and a per-day recovery rate of 0.1
#'  calc.trans.rates(lag = 14, max_duration = 180, recovery_rate = 0.1, eff_trans_rate = 0.1)
#'  # returns a value of 0.4371222
#' @export
calc.trans.rates <- function(lag,
                             max_duration,
                             recovery_rate,
                             raw_trans_rate=NULL,
                             eff_trans_rate=NULL) {
  if (is.null(raw_trans_rate) & is.null(eff_trans_rate)) {
    stop("Either 'raw_trans_rate' or 'eff_trans_rate' must be specified.")
  }
  if (!is.null(raw_trans_rate) & !is.null(eff_trans_rate)) {
    stop("'raw_trans_rate' and 'eff_trans_rate' cannot both be specified.")
  }

  trans_prob_vector <- c(rep(0, lag), rep(1, max_duration-lag))
  eff_trans_distr <- trans_prob_vector * dgeom(0:(max_duration-1), recovery_rate) # probability that infection resolves when it is idx days old
  if (!is.null(raw_trans_rate)) {
    eff_trans_rate <- raw_trans_rate * MESS::auc(1:max_duration,
                                                 eff_trans_distr,
                                                 type="spline")
    if (eff_trans_rate > 1) {
      stop(paste("Calculated effective transmission probability is", signif(eff_trans_rate, digits=3),"- cannot be greater than 1."))
    } else {
      return(eff_trans_rate)
    }
  }
  if (!is.null(eff_trans_rate)) {
    raw_trans_rate <- eff_trans_rate/MESS::auc(1:max_duration,
                                               eff_trans_distr,
                                               type="spline")
    if (raw_trans_rate > 1) {
      stop(paste("Calculated  transmission probability is", signif(raw_trans_rate, digits=3),"- cannot be greater than 1."))
    } else {
      return(raw_trans_rate)
    }
  }
}

#' @title Calculate Ross-Macdonald parameters
#' @description Calculates the Ross-Macdonald host and vector infection rates, as well as the host and vector infected proportions, for a given set of epidemiological parameters.
#' @usage calc.RM.params(N_h, N_v, bite_rate, eff_hv_trans_rate, eff_vh_trans_rate, h_rec_rate, v_rec_rate, h_relapse_rate)
#' @details
#' This function calculates the proportion of the host and vector populations that are infected at equilibrium from Ross-Macdonald equations. It is used by the main `run.RM()` function as well as plotting functions `plot.inf.over.time()` and `plot.inf.histogram()`. It is also useful for quickly checking if a given set of epidemiological parameters is plausible/possible.
#' @param N_h host population size (number of individuals)
#' @param N_v vector population size (number of individuals)
#' @param bite_rate mean number of hosts bitten per vector per day
#' @param eff_hv_trans_rate 'effective' host-to-vector transmission rate (see [calc.trans.rates()])
#' @param eff_vh_trans_rate 'effective' vector-to-host transmission rate (see [calc.trans.rates()])
#' @param h_rec_rate the per day probability that an infected host will recover (or die and be replaced)
#' @param v_rec_rate the per day probability that an infected vector will recover (or die and be replaced)
#' @param h_relapse_rate the per day probability of a relapse infection
#' @returns A named list where `h_inf_rate` is the host infection rate, `v_inf_rate` is the vector infection rate, `H_eq` is the expected proportion of the host population that is infected at equilibrium, and `V_eq` is the expected proportion of the vector population that is infected at equilibrium.
#' @examples
#' ## calculate Ross-Macdonald parameters corresponding to a host to vector ratio of 1:5, a per day vector bite rate of 0.3, an 'effective' host to vector transmission rate of 0.1, an 'effective' vector to host transmission rate of 0.2, a per day host recovery rate of 0.01, and a per day vector death rate of 0.05
#'  calc.RM.params(N_h = 1000, N_v = 5000, bite_rate = 0.2, eff_hv_trans_rate = 0.1, eff_vh_trans_rate = 0.2, h_rec_rate = 0.01, v_rec_rate = 0.05, h_relapse_rate=0.0001)
#'
#'  # returns the following values:
#'
#'  $h_inf_rate
#'  [1] 0.2001
#'
#'  $v_inf_rate
#'  [1] 0.02
#'
#'  $H_eq
#'  [1] 0.8334127
#'
#'  $V_eq
#'  [1] 0.2500178
#' @export
calc.RM.params <- function(N_h,
                           N_v,
                           bite_rate,
                           eff_hv_trans_rate,
                           eff_vh_trans_rate,
                           h_rec_rate,
                           v_rec_rate,
                           h_relapse_rate=NULL) {
  if (is.null(h_relapse_rate)) {
    h_relapse_rate = 0
  }
  h_inf_rate = bite_rate * (N_v/N_h) * eff_vh_trans_rate
  v_inf_rate = bite_rate * eff_hv_trans_rate
  H_eq = (h_inf_rate*v_inf_rate - h_rec_rate*v_rec_rate)/(h_inf_rate*v_inf_rate + v_inf_rate*h_rec_rate)
  V_eq = (h_inf_rate*v_inf_rate - h_rec_rate*v_rec_rate)/(h_inf_rate*v_inf_rate + h_inf_rate*v_rec_rate)
  if (H_eq < 0) {
    stop("Calculated host equilibrium is < 0. Please choose different starting parameters.")
  } else if (V_eq < 0) {
    stop("Calculated vector equilibrium is < 0. Please choose different starting parameters.")
  } else {
    return(list(h_inf_rate=h_inf_rate,
                v_inf_rate=v_inf_rate,
                H_eq=H_eq,
                V_eq=V_eq))
  }
}

#' @title Calculate the basic reproductive number (R0)
#' @description Calculates the basic reproductive number (R0) for a given Ross-Macdonald simulation for either hosts or vectors.
#' @usage calc.R0(RM_out, population, epoch)
#' @details
#' This function returns the average number of hosts that an infected host goes on to infect OR the average number of vectors that an infected vector goes on to infect. If specified, the function will only calculate this value for infections that started within a particular time range (epoch).
#' @param RM_out named list containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param population either "H" or "V" for the host or vector population, respectively
#' @param epoch a vector of length 2 specifying the time range that should be considered - defaults to the entire simulated time period
#' @returns The R0 value
#' @examples
#' ## calculate R0 value for hosts in the example dataset
#'  calc.R0(RM_out = RMsim::sim3, population="H")
#'
#'  # returns:
#'
#'  [1] 1.000079

#'  calc.R0(RM_out = RMsim::sim3, population="H", epoch=c(1, 300))
#'
#'  # returns:
#'
#'  [1] 2.201474
#' @export
calc.R0 <- function(RM_out,
                    population = NULL,
                    epoch = NULL) {
  `%ni%` <- Negate(`%in%`)
  if (is.null(population)) {
    stop("\'population\' must be either \'H\' (host) or \'V\' (vector)")
  } else {
    if (population %ni% c("H", "V")) {
      stop("\'population\' must be either \'H\' (host) or \'V\' (vector)")
    }
  }
  if (is.null(epoch)) {
    epoch <- range(RM_out$input_parameters["t_start",], RM_out$input_parameters["t_end",])
  } else {
    if (length(epoch) != 2 | !is.numeric(epoch)) {
      stop("\'epoch\' must be a numeric vector of length 2, specifying the start and end time step.")
    }
  }
  infection_record <- RM_out$infection_record
  infected <- infection_record[intersect(which(start_t > min(epoch) & start_t < max(epoch)), grep(population, infected)), "inf_id"]
  R0 <- mean(sapply(infected$inf_id, n.infs, infection_record=infection_record))
  return(R0)
}

# calculates how many 2nd degree infections result from a particular infection
n.infs <- function(infection_record, inf) {
  first.degree <- infection_record[origin_inf==as.integer(inf), "inf_id"]
  if (length(first.degree$inf_id) > 0) {
    n <-  nrow(infection_record[origin_inf %in% first.degree$inf_id, "inf_id"])
  } else {
    n <- 0
  }
  return(n)
}
