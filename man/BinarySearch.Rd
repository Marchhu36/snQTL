% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/symmPMD.R
\name{BinarySearch}
\alias{BinarySearch}
\title{Search soft threshold}
\usage{
BinarySearch(argv, sumabsv, maxiter = 150)
}
\arguments{
\item{argv}{the vector to be soft thresholded}

\item{sumabsv}{upperbound of the L_1 norm of sv}

\item{maxiter}{max iteration to perform binary search}
}
\value{
the proper threshold level \code{lamv}.
}
\description{
A binary search to find proper soft threshold \code{lamv} such that
\deqn{sv = soft(argv, lamv) / ||soft(argv, lamv)||_2, ||sv||_1 = sumabsv}
}
\seealso{
\code{symmPMD()}.
}
