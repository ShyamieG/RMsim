#' @title Plot proportion infected over time
#' @description Generates a plot of the actual and expected proportions of the host and vector populations that are infected over the duration of the simulation given the output of `run.RM()`.
#' @usage plot.inf.over.time(RM_out, scale_y_axis = T, plot_hosts = T, plot_vectors = T, h_col = "#9A5EA1", v_col = "#98823C")
#' @details
#' asdf.\cr\cr
#' asdf.
#' @param RM_out object containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param scale_y_axis if TRUE, the y-axis is scaled according to the values in RM_out. If FALSE, the y-axis range is set to c(0,1). Default is TRUE.
#' @param plot_hosts if TRUE, the host values are plotted. Default is TRUE.
#' @param plot_vectors if TRUE, the vector values are plotted. Default is TRUE.
#' @param h_col desired colour for plotting the host data
#' @param v_col desired colour for plotting the vector data
#' @examples
#' ## plot the results of the basic Ross-Macdonald simulation example
#'  plot.inf.over.time(RM_out_basic_sim)
#'  plot.inf.over.time(RM_out_basic_sim, scale_y_axis = F)
#'  plot.inf.over.time(RM_out_basic_sim, plot_hosts = F)
#' @export
plot.inf.over.time <- function(RM_out, scale_y_axis=T, plot_hosts=T, plot_vectors=T, h_col="#9A5EA1", v_col="#98823C") {
  indiv_status <- RM_out$indiv_status
  input_params <- RM_out$input_parameters
  RM_params <- RM_out$RM_parameters

  # Check that something will be plotted
  if (plot_hosts==F & plot_vectors==F) {
    stop("At least one of 'plot_hosts' and 'plot_vectors' must be TRUE.")
  }

  # Calculate y limit
  if (scale_y_axis==T) {
    Ylim <- c(NA, NA)
    pops <- c()
    if (plot_hosts==T) {pops <- c(pops, "H")}
    if (plot_vectors==T) {pops <- c(pops, "V")}
    for (i in 1:ncol(RM_params)) {
      for (j in pops) {
        Ylim <- c(min(c(Ylim[1], apply(indiv_status[grep(j, rownames(indiv_status)), RM_params["t_start", i]:RM_params["t_end", i]], MARGIN=2, FUN=sum)/input_params[paste0("N_", tolower(j)), i]), na.rm=T), max(c(Ylim[2], apply(indiv_status[grep(j, rownames(indiv_status)), RM_params["t_start", i]:RM_params["t_end", i]], MARGIN=2, FUN=sum)/input_params[paste0("N_", tolower(j)), i]), na.rm=T))
      }
    }
    Ylim <- c(floor(Ylim[1] * 5)/5, ceiling(Ylim[2] * 5)/5)
  } else {
    Ylim <- c(0,1)
  }

  # Plot each phase of the simulation
  par(mfrow=c(1,1))
  par(mgp=c(2.5,0.75,0))
  plot(c(1, ncol(indiv_status)), c(0,1), type="n", main="", ylim=c(min(Ylim), max(Ylim)), xlab="Time step", ylab="Proportion of individuals infected", las=1)
  for (i in 1:ncol(input_params)) {
    # Load parameters for ith phase
    for (j in 1:nrow(RM_params)) {
      assign(rownames(RM_params)[j], RM_params[j,i])
    }
    for (j in 1:nrow(input_params)) {
      assign(rownames(input_params)[j], input_params[j,i])
    }
    if (plot_hosts==T) {
      RM_out.h <- apply(indiv_status[grep("H", rownames(indiv_status)), t_start:t_end], MARGIN=2, FUN=sum)/N_h
      segments(x0=t_start, x1=t_end, y0=H_eq, y1=H_eq, lty=2, col=h_col)
      points(t_start:t_end, RM_out.h, type="l", col=h_col)
      if (i==1) {
        text(x=0.001*ncol(indiv_status), y=0.95*max(Ylim), "hosts", col=h_col, pos=4)
      }
    }
    if (plot_vectors==T) {
      RM_out.v <- apply(indiv_status[grep("V", rownames(indiv_status)), t_start:t_end], MARGIN=2, FUN=sum)/N_v
      segments(x0=t_start, x1=t_end, y0=V_eq, y1=V_eq, lty=2, col=v_col);par(new=T)
      points(t_start:t_end, RM_out.v, type="l", col=v_col)
      if (i==1) {
        text(x=0.001*ncol(indiv_status), y=0.95*max(Ylim) - 0.05*max(Ylim), "vectors", col=v_col, pos=4)
      }
    }
  }
}

#' @title Plot proportion infected at equilibrium
#' @description asdf
#' @usage asdf
#' @details asdf
#' asdf.\cr\cr
#' asdf.
#' @param RM_out object containing the results of a Ross-Macdonald simulation (output of `run.RM()`)
#' @param days_of_burn_in description
#' @param scale_x_axis description
#' @param plot_hosts description
#' @param plot_vectors description
#' @param phases description
#' @param observed_col description
#' @examples
#' ## plot the results of the basic Ross-Macdonald simulation example
#'  plot.inf.histogram(RM_out_basic_sim, days_of_burn_in = 400)
#'  plot.inf.histogram(RM_out_basic_sim, scale_y_axis = F)
#'  plot.inf.histogram(RM_out_basic_sim, plot_hosts = F)
#' @export
plot.inf.histogram <- function(RM_out,
                               days_of_burn_in=NULL,
                               scale_x_axis=T,
                               plot_hosts=T,
                               plot_vectors=T,
                               phases=NULL,
                               observed_col="red") {
  indiv_status <- RM_out$indiv_status
  input_params <- RM_out$input_parameters
  RM_params <- RM_out$RM_parameters

  # Check that something will be plotted
  if (plot_hosts==F & plot_vectors==F) {
    stop("At least one of 'plot_hosts' and 'plot_vectors' must be TRUE.")
  }

  # Check that a burn-in period has been defined
  if (is.null(days_of_burn_in)) {
    stop("User must specify 'days_of_burn_in'.")
  } else if (ncol(indiv_status) < days_of_burn_in) {
    stop("Burn-in period was not completed for this simulation.")
  }

  # Specify how to scale x-axis, if necessary
  if (scale_x_axis==T) {
    pops <- c()
    if (plot_hosts==T) {pops <- c(pops, "H")}
    if (plot_vectors==T) {pops <- c(pops, "V")}
  }

  # Check which phases to plot
  if (is.null(phases)) {
    phases <- 1:ncol(input_params)
  }

  # Plot each phase of the simulation
  for (i in phases) {
    # Load parameters for ith phase
    for (j in 1:nrow(RM_params)) {
      assign(rownames(RM_params)[j], RM_params[j,i])
    }
    for (j in 1:nrow(input_params)) {
      assign(rownames(input_params)[j], input_params[j,i])
    }

    # Make sure burn-in is not too long
    if (t_start+days_of_burn_in > t_end) {
      stop("Burn-in period cannot be longer than the length of the phase")
    }

    # Set up plot format
    if (plot_hosts==T & plot_vectors==T) {
      par(mfrow=c(2,1))
    } else {
      par(mfrow=c(1,1))
    }

    # Calculate x limit
    if (plot_hosts==T) {RM_out.h <- apply(indiv_status[grep("H", rownames(indiv_status)), (t_start+days_of_burn_in):t_end], sum, MARGIN=2)/N_h}
    if (plot_vectors==T) {RM_out.v <- apply(indiv_status[grep("V", rownames(indiv_status)), (t_start+days_of_burn_in):t_end], sum, MARGIN=2)/N_v}
    if (scale_x_axis==T) {
      Xlim <- c(NA, NA)
      pops <- c()
      if (plot_hosts==T) {pops <- c(pops, "H")}
      if (plot_vectors==T) {pops <- c(pops, "V")}
      for (j in pops) {
        Xlim <- c(min(c(Xlim[1], get(paste0("RM_out.", tolower(j)))), na.rm=T), max(c(Xlim[2], get(paste0("RM_out.", tolower(j)))), na.rm=T))
      }
      Xlim <- c(floor(Xlim[1] * 5)/5, ceiling(Xlim[2] * 5)/5)
    } else {
      Xlim <- c(0,1)
    }
    if (plot_hosts==T) {
      hist(RM_out.h, breaks=seq(min(Xlim), max(Xlim), 0.01), main=paste("Phase", i), xlab="Proportion of hosts infected", las=1)
      abline(v=mean(RM_out.h), lty=2, col=observed_col)
      abline(v=c(mean(RM_out.h)-sd(RM_out.h), mean(RM_out.h)+sd(RM_out.h)), col=observed_col)
      abline(v=H_eq, lty=2)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(RM_out.h, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.95, "empirical", col=observed_col)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(RM_out.h, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.85, "expected")
    }
    if (plot_vectors==T) {
      hist(RM_out.v, seq(min(Xlim), max(Xlim), 0.01), main=paste("Phase", i), xlab="Proportion of vectors infected")
      abline(v=mean(RM_out.v), lty=2, col=observed_col)
      abline(v=c(mean(RM_out.v)-sd(RM_out.v), mean(RM_out.v)+sd(RM_out.v)), col=observed_col)
      abline(v=V_eq, lty=2)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(RM_out.v, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.95, "empirical", col=observed_col)
      text(x=Xlim[1]+0.1*(Xlim[2]-Xlim[1]), y=max(hist(RM_out.v, breaks=seq(min(Xlim), max(Xlim), 0.01), plot=F)$counts)*0.85, "expected")
    }
  }
}
