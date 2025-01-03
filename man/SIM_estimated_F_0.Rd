% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{SIM_estimated_F_0}
\alias{SIM_estimated_F_0}
\title{Parametric method for estimating the null distribution}
\usage{
SIM_estimated_F_0(t, p_theta, SIM_F_0_method, mu_0_theta, sigma_0_theta)
}
\arguments{
\item{t}{vector of thresholds for SIM p-values p(theta)}

\item{p_theta}{vector of 1D SIM p-values; here is p(theta)}

\item{SIM_F_0_method}{1: parametric method; 2: nonparametric method}

\item{mu_0_theta}{the mean of the single-index p-value}

\item{sigma_0_theta}{the sd of the single-index p-valu}
}
\value{
the null distribution, same size as t
}
\description{
Parametric method for estimating the null distribution
}
