% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_diffnet_from_exp}
\alias{get_diffnet_from_exp}
\title{The differential matrix}
\usage{
get_diffnet_from_exp(X, Y, adj.beta = -1, trans = FALSE, location = NULL)
}
\arguments{
\item{X}{n1-by-p matrix for samples from the first population.
Rows are samples/observations, while columns are the features.}

\item{Y}{n2-by-p matrix for samples from the second population.
Rows are samples/observations, while columns are the features.}

\item{adj.beta}{Power to transform correlation matrices to weighted adjacency matrices
by \eqn{N_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson's correlation.
When \code{adj.beta=0}, the correlation marix is used.
When \code{adj.beta<0}, the covariance matrix is used.
The default value is \code{adj.beta=-1}.}

\item{trans}{logic variable, whether to only consider the trans-correlation (between genes from two different chromosomes or regions); see "details"}

\item{location}{vector, the (chromosome) locations for items}
}
\value{
The p-by-p differential matrix \eqn{D = N(Y) - N(X)}.
}
\description{
Given observations from two populations X and Y,
compute the differential matrix
\deqn{D = N(Y) - N(X)}
where N() is the covariance matrix, or the weighted adjacency matrices defined as
\deqn{N_{ij} = |corr(i, j)|^beta}
for some constant beta > 0, 1 <= i, j <= p.
Let N represent the regular correlation matrix when beta=0, and covariance matrix when beta<0.
}
