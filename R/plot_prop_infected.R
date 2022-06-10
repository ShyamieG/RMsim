#' @export
plot.prop.infected <- function(dat) {
  # Load parameters
  for (i in 1:length(dat$model_parameters)) {
    assign(names(dat$model_parameters)[i], dat$model_parameters[[i]])
  }
  for (i in 1:length(dat$call)) {
    assign(names(dat$call)[i], dat$call[[i]])
  }
  indiv_status <- dat$indiv_status
  # Create vectors of host and vector names for convenience
  hIDs <- paste0(rep("H", N_h), 1:N_h)
  vIDs <- paste0(rep("V", N_v), 1:N_v)
  par(mgp=c(2.5,0.75,0))
  plot(c(1,runtime), c(0,1), type="n", main="", xlab="Time step", ylab="Proportion of individuals infected")
  abline(h=H_eq, lty=2, col="#9A5EA1")
  abline(h=V_eq, lty=2, col="#98823C");par(new=T)
  plot(1:runtime, apply(indiv_status[hIDs,], MARGIN=2, FUN=sum)/N_h, type="l", ylim=c(0,1), axes=F, ann=F, col="#9A5EA1");par(new=TRUE)
  plot(1:runtime, apply(indiv_status[vIDs,], MARGIN=2, FUN=sum)/N_v, type="l", ylim=c(0,1), axes=F, ann=F, col="#98823C")
  text(x=rep(-3,2), y=c(0.95, 0.90), c("hosts", "vectors"), col=c("#9A5EA1","#98823C"), pos=4)
  text(x=runtime*0.25, y=c(H_eq, V_eq)+0.02, labels=c(paste("host equil =", signif(H_eq, 3)), paste("vector equil =", signif(V_eq, 3))), pos=2, cex=0.85)
}
