# Bricks for snQTL

# normalization -------
#' Normalization
#'
#' @description
#' Normalize the vector to have L2 norm 1
#'
#' @param x vector, the numerical vector need to be normalized
#'
#' @return the normalized vector

normalize <- function(x) {
  return(x / sqrt(sum(x^2)))
}

# permute indices -----------

#' Permute indices
#'
#' @param n1 number of observations from first population
#' @param n2 number of observations from second population
#' @param n3 number of observations from third population
#'
#' @return A list with the following components:
#'  \item{i1}{a numeric vector with length \code{n1}, the indices for permuted first population}
#'  \item{i2}{a numeric vector with length \code{n2}, the indices for permuted second population}
#'  \item{i3}{a numeric vector with length \code{n3}, the indices for permuted second population}
#'

multi_permuteIndex <- function(n1, n2, n3){
  i.sample <- sample(1:(n1+n2+n3), replace=FALSE)
  return(list(i1=i.sample[1:n1], i2=i.sample[(n1+1):(n1+n2)], i3 = i.sample[(n1+n2+1):(n1+n2+n3)]))
}


# from sLED ------------

# get statistics --------------

#' Calculate of sLME for matrices
#'
#' @description
#' Calculate the sLME given a matrix \eqn{D}.
#' For any symmetric matrix \eqn{D}, sLME test statistic is defined as
#' \deqn{max{ sEig(D), sEig(-D) }}
#' where \code{sEig()} is the sparse leading eigenvalue, defined as
#' \deqn{$max_{v} v^T A v$}{max_{v} t(v)*A*v}
#' subject to
#' \eqn{$||v||_2 \leq 1, ||v||_1 \leq s$}{||v||_2 <= 1, ||v||_1 <= s}.
#'
#' @param Dmat p-by-p numeric matrix, the differential matrix
#' @param rho a large positive constant such that \eqn{D+diag(rep(rho, p))} and \eqn{-D+diag(rep(rho, p))}
#'        are positive definite.
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters, each between \eqn{1/sqrt(p)} and 1.
#'        Each sumabs*\eqn{$sqrt(p)$}{sqrt(p)} is the upperbound of the L_1 norm of leading sparse eigenvector \eqn{v}.
#' @param niter the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#'
#' @return A list containing the following components:
#'  \item{sumabs.seq}{the sequence of sparsity parameters}
#'  \item{rho}{a positive constant to augment the diagonal of the differential matrix
#'            such that \eqn{D + rho*I} becomes positive definite.}
#'  \item{stats}{a numeric vector of test statistics when using different sparsity parameters
#'            (corresponding to \code{sumabs.seq}).}
#'  \item{sign}{a vector of signs when using different sparsity parameters (corresponding to \code{sumabs.seq}).
#'          Sign is "pos" if the test statistic is given by sEig(D), and "neg" if is given by sEig(-D),
#'          where \code{sEig} denotes the sparse leading eigenvalue.}
#'  \item{v}{the sequence of sparse leading eigenvectors, each row corresponds to one sparsity
#'          parameter given by \code{sumabs.seq}.}
#'  \item{leverage}{the leverage score for genes (defined as \eqn{v^2} element-wise) using
#'          different sparsity parameters. Each row corresponds to one sparsity
#'          parameter given by \code{sumabs.seq}.}
#'
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Covariance Matrices,
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.

sLEDTestStat <- function(Dmat, rho=1000, sumabs.seq=0.2,
                         niter=20, trace=FALSE) {
  ndim <- 1 ## only consider the first sparse eigenvector
  p <- ncol(Dmat)
  ntest <- length(sumabs.seq)

  results <- list()
  results$sumabs.seq <- sumabs.seq
  results$rho <- rho

  results$stats <- rep(NA, ntest)
  results$sign <- rep(NA, ntest)
  results$v <- matrix(NA, nrow=ntest, ncol=p)
  results$leverage <- matrix(NA, nrow=ntest, ncol=p)

  ## for each sparsity parameter
  for (i in 1:ntest) {
    sumabs <- sumabs.seq[i]
    pos.out <- symmPMD(Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)
    neg.out <- symmPMD(- Dmat + rho * diag(p),
                       sumabs=sumabs, trace=trace, niter=niter)

    if (pos.out$d >= neg.out$d) {
      results$sign[i] <- "pos"
      results$stats[i] <- pos.out$d - rho
      results$v[i, ] <- pos.out$v
      results$leverage[i, ] <- (pos.out$v)^2
    } else {
      results$sign[i] <- "neg"
      results$stats[i] <- neg.out$d - rho
      results$v[i, ] <- neg.out$v
      results$leverage[i, ] <- (neg.out$v)^2
    }
  }

  return(results)
}

# get difference matrix ----------

#' The differential matrix
#'
#' @description Given observations from two populations X and Y,
#' compute the differential matrix
#' \deqn{D = N(Y) - N(X)}
#' where N() is the covariance matrix, or the weighted adjacency matrices defined as
#' \deqn{N_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p.
#' Let N represent the regular correlation matrix when b=0, and covariance matrix when b<0.
#'
#' @param X n1-by-p matrix for samples from the first population.
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population.
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta Power to transform correlation matrices to weighted adjacency matrices
#'        by \eqn{N_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson's correlation.
#'        When \code{adj.beta=0}, the correlation marix is used.
#'        When \code{adj.beta<0}, the covariance matrix is used.
#'        The default value is \code{adj.beta=-1}.
#'
#' @return The p-by-p differential matrix \eqn{D = N(Y) - N(X)}

getDiffMatrix <- function(X, Y, adj.beta=-1) {
  if (adj.beta < 0) {
    Dmat <- cov(Y) - cov(X)
  } else if (adj.beta == 0) {
    Dmat <- cor(Y) - cor(X)
  } else {
    Dmat <- abs(cor(Y))^adj.beta - abs(cor(X))^adj.beta
  }
  return(Dmat)
}

# get trans differential matrix

#' The differential matrix with trans-correlation
#'
#' @description Given observations from two populations X and Y,
#' compute the differential matrix
#' \deqn{D = N(Y) - N(X)}
#' where N() is the covariance matrix, or the weighted adjacency matrices defined as
#' \deqn{N_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p.
#' Let N represent the regular correlation matrix when b=0, and covariance matrix when b<0.
#'
#' Then, we let \eqn{N_{ij} = 0} if \eqn{i, j} are from the same region.
#' Under gene co-expression context, trans-correlation usually refer to the correlation between
#' two genes \eqn{i, j} from two chromosomes.
#'
#' @param X n1-by-p matrix for samples from the first population.
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population.
#'        Rows are samples/observations, while columns are the features.
#' @param location vector, the (chromosome) locations for items
#' @param adj.beta Power to transform correlation matrices to weighted adjacency matrices
#'        by \eqn{N_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson's correlation.
#'        When \code{adj.beta=0}, the correlation marix is used.
#'        When \code{adj.beta<0}, the covariance matrix is used.
#'        The default value is \code{adj.beta=-1}.
#'
#' @return The p-by-p differential matrix \eqn{D = N(Y) - N(X)} with only trans-correlation


get_trans_diffmatrix = function(X, Y, location, adj.beta = -1){

  # location: the locations of each gene, only positive numbers refer to a valid location

  if (adj.beta < 0) {
    Dmat <- cov(Y) - cov(X)
  } else if (adj.beta == 0) {
    Dmat <- cor(Y) - cor(X)
  } else {
    Dmat <- abs(cor(Y))^adj.beta - abs(cor(X))^adj.beta
  }

  # delete the cis-correlation
  locs = unique(location)[unique(location) >= 0]
  for (i in locs) {
    ind_vec = location == i
    Dmat[ind_vec, ind_vec] = 0
  }

  return(Dmat)
}
