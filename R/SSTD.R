# Sparse Symmetric Rank-1 Tensor Decomposition
# Feb 26 2023 Jiaxin
# Last edition: Jul 17 add option to hide message

# Modified algorithm of ``Provable Sparse Tensor Decomposition" Sun et.al (2019)

# # dependencies -------
# # library(rTensor)
# source("tensor_class.R")
# source("tensor_product.R")
# source("symmPMD.R")
# source("bricks.R")


# main function ---------

#' Sparse Symmetric Tensor Decomposition (SSTD)
#'
#' @description
#' SSTD solves the rank-1 approximation to the a p-by-p-by-q sparse symmetric tensor \eqn{\mathcal{D}}:
#' \deqn{ \min_{\Lambda, v, u} ||\mathcal{D} - \Lambda v \circ v \circ u||_F^2}
#' subject to
#' \deqn{ \Lambda > 0, v \in R^p, u \in R^q, ||v||_2 = ||u||_2 = 1, ||v||_0 <= R}
#'
#' The solution \eqn{\Lambda} is the sparse leading tensor eigenvalue (sLTE),
#' \eqn{v} is the sparse leading tensor eigenvector, and \eqn{u} is the loading vector.
#'
#' The Symmetric Penalized Matrix Decomposition \code{symmPMD()} is used in the iterative algorithm.
#'
#'
#' @param T_obs array, a p-by-p-by-q tensor; each p-by-p layer in \code{T_obs} should be symmetric
#' @param u_ini vector, with length q; the random initialization for loading vector
#' @param v_ini vector, with length p; the random initialization for tensor eigenvector
#' @param max_iter integer, the maximal iteration number
#' @param sumabs number, the number specify the sparsity level in the matrix/tensor eigenvector; \code{sumabs} takes value between \eqn{1/sqrt(p)} and 1, where \eqn{p} is the dimension; \code{sumabs}\eqn{*sqrt(p)} is the upperbound of the L1 norm of the leading matrix/tensor eigenvector (see \code{symmPMD()})
#' @param niter integer, the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param rho number, a large positive constant adding to the diagonal elements to ensure positive definiteness in symmetric matrix spectral decomposition
#' @param tol number, the tolerance threshold for SSTD convergence; if the error difference between two iterations is smaller than \code{tol}, then we stop the iteration and consider the algorithm converges
#' @param verbose logic variable, whether to print the progress during permutation tests
#'
#' @return a list containing the following:
#'
#' \item{u_hat}{vector, with length q; the estimated loading vector}
#'
#' \item{v_hat}{vector, with length p; the estimated tensor eigenvector}
#'
#' \item{gamma_hat}{number, the estimated sLTE \eqn{\Lambda}}
#'
#' @references Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2024).
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." bioRxiv, 2024-03.
#' @references Sun, W. W., Lu, J., Liu, H., & Cheng, G. (2017). "Provable sparse tensor decomposition."
#' Journal of the Royal Statistical Society Series B: Statistical Methodology, 79(3), 899-916.
#'
#' @seealso \code{symmPMD()}
#'
#' @export

SSTD_R1 <- function(T_obs, u_ini, v_ini, max_iter = 20, sumabs = 0.5, niter = 20, rho = 1000, tol = 10^(-3), verbose = FALSE) {
  # sparse symmetric rank-1 tensor decomposition
  # pre-specified initialization

  # T_obs: observed tensor
  # max_iter: maximal number of iteration
  # sumabs, niter, rho: sparse parameters, iteration number, and diagonal element in covariance for symmPMD()
  # tol: tolerance to end the iteration before max_iter loops

  t <- dim(T_obs)[3]
  p <- dim(T_obs)[1]

  u <- normalize(u_ini)
  v <- normalize(v_ini)

  err <- c()

  if(verbose){
    cat("iteration: ")
  }

  for (j in 1:max_iter) { # iteration
    if(verbose){
      cat(j, ", ")
    }
    # update u
    u_new <- normalize(as.vector(ttl(as.tensor(T_obs), list(t(as.matrix(v)), t(as.matrix(v))), ms = c(1, 2))@data))

    # update v
    symmT <- ttm(as.tensor(T_obs), t(as.matrix(u_new)), m = 3)@data[, , 1]
    pos.out <- symmPMD(symmT + rho * diag(p),
      sumabs = sumabs, trace = F, niter = niter
    )
    neg.out <- symmPMD(-symmT + rho * diag(p),
      sumabs = sumabs, trace = F, niter = niter
    )

    if (pos.out$d >= neg.out$d) {
      v_new <- pos.out$v
    } else {
      v_new <- neg.out$v
    }

    v_new <- as.vector(normalize(v_new))

    err <- c(err, sum((v %o% v %o% u - v_new %o% v_new %o% u_new)^2))

    if (abs(err[j]) <= tol) {
      v <- v_new
      u <- u_new

      if(verbose){
        cat("converge, ")
      }
      break
    }

    v <- v_new
    u <- u_new
  } # end j

  # calculate tensor singular value
  gamma = as.numeric(ttl(as.tensor(T_obs), list(t(as.matrix(v)), t(as.matrix(v)), t(as.matrix(u))), ms = c(1,2,3))@data)

  if(verbose){
    cat("stop \n")
  }


  return(list(u_hat = u, v_hat = v, gamma_hat = gamma))
}

