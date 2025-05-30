% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{sample.RM}
\alias{sample.RM}
\title{Sample infected individuals from a Ross-Macdonald simulation}
\usage{
sample.RM(RM_out, time_step, population, proportion = NULL, number = NULL, sample_post_lag = TRUE, resample_possible = TRUE, sort_events = TRUE)
}
\arguments{
\item{RM_out}{object containing the results of a Ross-Macdonald simulation (output of \code{run.RM()})}

\item{time_step}{the time step (i.e. simulation day) at which infected samples should be drawn}

\item{population}{either 'H' (host), 'V' (vector), or 'both'; the population to draw samples from}

\item{proportion}{if specified, the proportion of the total infected population to sample}

\item{number}{if specified, the number of infected individuals to sample}

\item{sample_post_lag}{if TRUE, only individuals that are beyond the infection lag phase can be sampled. Default is TRUE.}

\item{resample_possible}{if TRUE, the same individual can be sampled multiple times. Default is FALSE.}

\item{sort_events}{if TRUE, the updated infection record with new sampling events will be sorted by \code{start_t}. Default is TRUE.}
}
\description{
Adds events to the infection record representing samples taken from infected individuals.
}
\details{
The sampling of infected individuals is represented by a unique type of event that is included in the infection record. This function allows users to define the sampling scheme by specifying the time step when infected individuals should be sampled, which population (host, vector, or both) to draw samples from, whether a certain proportion of the infected population or a fixed number of individuals should be sampled, and if samples should only be drawn from individuals that are beyond the 'lag' period of their infection. Multiple samples can be drawn from the same individual if \code{resample_possible} is set to TRUE. If \code{sort_events} is set to TRUE (which it is by default) the updated infection record will be sorted by \code{start_t} (infection start time or sample time).
}
\examples{
## add samples to the basic Ross-Macdonald simulation example
 set.seed(123456)
 RM_sampled <- sample.RM(RM_out = RMsim::sim3, time_step = 1000, population = "H", number = 5, sort_events = FALSE) # sample 5 host infections on day 1000
 tail(RM_sampled$infection_record) # sample events are added to the end of the infection record because sort_events was set to FALSE
}
