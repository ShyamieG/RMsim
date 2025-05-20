#' @title Plot number or proportion of infected individuals over time
#' @description Generates a plot of the actual and expected numbers or proportions of hosts that are infected over the simulated time period.
#' @usage plot_inf_over_time(RM_out, type = "proportion", scale_y_axis = TRUE, plot_hosts = TRUE, plot_vectors = TRUE, h_col = "#9A5EA1", v_col = "#98823C")
#' @details
#' Parses the result of `run.RM()` to plot the proportion of host and vector population that are infected at each time step of the simulation. The expected values of these parameters, calculated from the underlying Ross-Macdonald equations, are plotted as dashed lines.
#' @param RM_out object containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param type what kind of y-values to plot; either `number` or `proportion`
#' @param scale_y_axis if TRUE, the y-axis is scaled according to the values in RM_out. If FALSE, the y-axis range is set to c(0,1). Default is TRUE.
#' @param plot_hosts if TRUE, the host values are plotted. Default is TRUE.
#' @param plot_vectors if TRUE, the vector values are plotted. Default is TRUE.
#' @param h_col desired colour for plotting the host data
#' @param v_col desired colour for plotting the vector data
#' @examples
#' ## plot the results of the basic Ross-Macdonald simulation example
#'  plot_inf_over_time(RMsim::sim3, type="number")
#'  plot_inf_over_time(RMsim::sim3, type="proportion")
#'  plot_inf_over_time(RMsim::sim3, type="proportion", scale_y_axis = FALSE)
#'  plot_inf_over_time(RMsim::sim3, type="proportion", plot_hosts = FALSE)
#' @export
plot_inf_over_time <- function(RM_out,
                               type=NULL,
                               scale_y_axis=T,
                               plot_hosts=T,
                               plot_vectors=T,
                               h_col="#9A5EA1",
                               v_col="#98823C") {
  # Check that something will be plotted
  if (plot_hosts==F & plot_vectors==F) {
    stop("At least one of 'plot_hosts' and 'plot_vectors' must be TRUE.")
  }
  if (is.null(type)) {
    stop("'type' of plot to generate must be specified, either 'number' or 'proportion'.")
  }

  input_parameters <- RM_out$input_parameters
  RM_params <- RM_out$RM_parameters
  dat <- RM_out$proportion_infected
  if (type == "number") {
    for (i in 1:ncol(input_parameters)) {
      start <- input_parameters["t_start", i]
      end <- input_parameters["t_end", i]
      N_h <- input_parameters["N_h", i]
      N_v <-input_parameters["N_v", i]
      dat["H", start:end] <- dat["H", start:end] * N_h
      dat["V", start:end] <- dat["V", start:end] * N_v
    }
  }
  time.steps <- input_parameters["t_start", 1]:input_parameters["t_end", ncol(input_parameters)]

  # Calculate y limit
  if (scale_y_axis==T) {
    pops <- NULL
    if (plot_hosts==T) {pops <- c(pops, "H")}
    if (plot_vectors==T) {pops <- c(pops, "V")}
    Ylim <- c(min(dat[pops,]), max(dat[pops,]))
    Ylim <- c(floor(Ylim[1] * 20)/20, ceiling(Ylim[2] * 20)/20)
  } else {
    Ylim <- c(0,1)
  }

  # Plot observed values over entire simulation
  par(mfrow=c(1,1))
  par(mgp=c(2.5,0.75,0))
  plot(c(1, max(time.steps)), c(0,1), type="n", main="", ylim=c(min(Ylim), max(Ylim)), xlab="Time step", ylab="Proportion of individuals infected", las=1)
  if (plot_hosts==T) {
    points(time.steps, dat["H",], type="l", col=h_col, lwd=1.2)
  }
  if (plot_vectors==T) {
    points(time.steps, dat["V",], type="l", col=v_col, lwd=1.2)
  }

  # Plot expected equilibrium values for each phase of the simulation
  for (i in 1:ncol(input_parameters)) {
    # Load equilibrium values for ith phase
    for (j in 1:nrow(RM_params)) {
      assign(rownames(RM_params)[j], RM_params[j,i])
    }
    for (j in 1:nrow(input_parameters)) {
      assign(rownames(input_parameters)[j], input_parameters[j,i])
    }
    if (plot_hosts==T) {
      if (type=="number") {H_eq <- H_eq * N_h}
      segments(x0=t_start, x1=t_end, y0=H_eq, y1=H_eq, lty=2, col=h_col, lwd=1.5)
      if (i==1) {
        text(x=0.001*length(time.steps), y=H_eq + 0.04*(Ylim[2]-Ylim[1]), "hosts", col=h_col, pos=4)
      }
    }
    if (plot_vectors==T) {
      if (type=="number") {V_eq <- V_eq * N_v}
      segments(x0=t_start, x1=t_end, y0=V_eq, y1=V_eq, lty=2, col=v_col, lwd=1.5)
      if (i==1) {
        text(x=0.001*length(time.steps), y=V_eq + 0.04*(Ylim[2]-Ylim[1]), "vectors", col=v_col, pos=4)
      }
    }
  }
}

#' @title Plot infected proportion at equilibrium
#' @description Generates a histogram of the actual proportions of the host and vector populations that are infected following a specified burn-in period.
#' @usage plot_inf_histogram(RM_out, days_of_burn_in, scale_x_axis = TRUE, plot_hosts = TRUE, plot_vectors = TRUE, phases, observed_col = "red")
#' @details
#' Parses the result of `run.RM()` to plot a histogram of the proportion of host and vector population that are infected at each time step of the simulation. Note that these data will be substantially auto-correlated. Nevertheless, this function allows users to visualize how well the results of the stochastic simulation match with the theoretical expectation.\cr\cr
#' The `days_of_burn_in` and `phases` arguments allow users to specify the time period to be plotted. `days_of_burn_in` should be chosen so as to exclude the initial phase of a simulation when equilibrium has not yet been reached. Viewing the result of [plot_inf_over_time()] can help users choose this value.\cr\cr
#' Values related to the empirical distribution are plotted in the colour specified by `observed_col`; the dashed vertical line represents the mean, and the solid lines represent 1 standard deviation above and below the mean. The expected infected proportion of the population, calculated from the underlying Ross-Macdonald equations, is plotted in black.
#' @param RM_out object containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param days_of_burn_in the number of time steps to discard from the start of the simulation phase in calculating the empirical mean and standard deviation
#' @param scale_x_axis if TRUE, the x-axis is scaled according to the values in RM_out. If FALSE, the x-axis range is set to c(0,1). Default is TRUE.
#' @param plot_hosts if TRUE, values for the host population will be plotted
#' @param plot_vectors if TRUE, values for the vector population will be plotted
#' @param phases which phases of the Ross-Macdonald simulation, as defined by a unique set of epidemiological parameters, should be plotted. Each phase is plotted separately.
#' @param phases_per_plot the number of phases to plot on the same plotting space, allowing multiple phases to be visualized side by side. Default is 1.
#' @param observed_col the colour to plot values related to the empirical distribution
#' @examples
#' ## plot the results of the basic Ross-Macdonald simulation example
#'  plot_inf_histogram(RMsim::sim3, days_of_burn_in = 500, phases_per_plot=3)
#'  plot_inf_histogram(RMsim::sim3, days_of_burn_in = 500, phases_per_plot=3, scale_x_axis = FALSE)
#'  plot_inf_histogram(RMsim::sim3, days_of_burn_in = 500, phases_per_plot=3, plot_hosts = FALSE)
#' @export
plot_inf_histogram <- function(RM_out,
                               days_of_burn_in=NULL,
                               scale_x_axis=T,
                               plot_hosts=T,
                               plot_vectors=T,
                               phases=NULL,
                               phases_per_plot=1,
                               observed_col="red") {

  input_parameters <- RM_out$input_parameters
  RM_params <- RM_out$RM_parameters
  proportion_infected <- RM_out$proportion_infected
  time.steps <- input_parameters["t_start", 1]:input_parameters["t_end", ncol(input_parameters)]

  # Determine which populations will be plotted
  pops <- c()
  if (plot_hosts==T) {pops <- c(pops, "H")}
  if (plot_vectors==T) {pops <- c(pops, "V")}
  if (length(pops) == 0) {
    stop("At least one of 'plot_hosts' and 'plot_vectors' must be TRUE.")
  }

  # Check that a burn-in period has been defined and completed
  if (is.null(days_of_burn_in)) {
    stop("User must specify 'days_of_burn_in'.")
  } else if (max(time.steps) < days_of_burn_in) {
    stop("Burn-in period was not completed for this simulation.")
  }

  # Check which phases to plot
  if (is.null(phases)) {
    phases <- 1:ncol(input_parameters)
  } else {
    if (!is.numeric(phases)) {
      stop("\'phases\' must be a numeric vector.")
    }
  }

  # Plot each phase of the simulation
  par(mfrow=c(length(pops), phases_per_plot))
  for (i in phases) {
    # Load parameters for ith phase
    for (j in 1:nrow(RM_params)) {
      assign(rownames(RM_params)[j], RM_params[j,i])
    }
    for (j in 1:nrow(input_parameters)) {
      assign(rownames(input_parameters)[j], input_parameters[j,i])
    }

    # Make sure burn-in is not too long
    if (t_start+days_of_burn_in > t_end) {
      stop("Burn-in period cannot be longer than the length of the phase")
    }

    # Specify data to plot for this phase
    dat <- proportion_infected[pops, (t_start+days_of_burn_in):t_end]

    # Scale x axis, if specified
    if (scale_x_axis==T) {
      Xlim <- c(min(dat), max(dat))
      Xlim <- c(floor(Xlim[1] * 20)/20, ceiling(Xlim[2] * 20)/20)
    } else {
      Xlim <- c(0,1)
    }
    for (pop in pops) {
      type <- c("host", "vector")[match(pop, c("H","V"))]
      X <- as.numeric(dat[pop,])
      hist(X, breaks=seq(min(Xlim), max(Xlim), 0.01), main=paste("Phase", i), xlab=paste0("Proportion of ", type, "s infected"), las=1)
      abline(v=mean(X), lty=2, col=observed_col, lwd=1.5)
      abline(v=c(mean(X)-sd(X), mean(X)+sd(X)), col=observed_col, lwd=1.5)
      abline(v=get(paste0(pop, "_eq")), lty=2, lwd=1.5)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(X, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.95, "empirical", col=observed_col)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(X, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.85, "expected")
    }
  }
}
