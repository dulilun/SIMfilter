% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{SIM_threshold_for_p}
\alias{SIM_threshold_for_p}
\title{Function: threshold value for SIM p-valueS}
\usage{
SIM_threshold_for_p(
  alpha,
  p_theta,
  hat_pi_0_theta,
  SIM_F_0_method,
  mu_theta,
  sigma_theta,
  option
)
}
\arguments{
\item{alpha}{error rate}

\item{p_theta}{vector of SIM p-values}

\item{hat_pi_0_theta}{estimated \eqn{\pi_0}}

\item{SIM_F_0_method}{1: parametric method; 2: nonparametric method.}

\item{mu_theta}{\eqn{\mu_0(\theta)}}

\item{sigma_theta}{\eqn{\sigma_0(\theta)}}

\item{option}{options for choosing the grid}
}
\value{
t_alpha_SIM: threshold value for SIM p-values
}
\description{
Function: threshold value for SIM p-valueS
}
