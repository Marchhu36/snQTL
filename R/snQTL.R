# snQTL -- main functions

# # dependencies -------
# source("bricks.R")
# source("symmPMD.R")
# source("SSTD.R")
# source("tensor_class.R")
# source("tensor_product.R")

# Obtain test statistics with multiple networks -------------

#' Test statistics for snQTL
#'
#' @description
#' Generate snQTL test statistics from a given list of differential networks.
#' This function takes a list of differential networks, the choice of test statistics, and other computational tuning parameters as inputs.
#' Outputs include the calculated statistics, recall of the choice, and the decomposition components associated with the statistics.
#'
#' @param network_list list, a list of p-by-p differential networks
#' @param method character, the choice of test statistics; see "details"
#' @param rho number, a large positive constant adding to the diagonal elements to ensure positive definiteness in symmetric matrix spectral decomposition
#' @param sumabs number, the number specify the sparsity level in the matrix/tensor eigenvector; \code{sumabs} takes value between \eqn{1/sqrt(p)} and 1, where \eqn{p} is the dimension; \code{sumabs}\eqn{*sqrt(p)} is the upperbound of the L1 norm of the leading matrix/tensor eigenvector (see \code{symmPMD()})
#' @param niter integer, the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace logic variable, whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#' @param tensor_iter integer, the maximal number of iteration in SSTD algorithm (see \code{max_iter} in \code{SSTD()})
#' @param tensor_tol number, a small positive constant for error difference to indicate the SSTD convergence (see \code{tol} in \code{SSTD()})
#' @param tensor_seed number, the seed to generate random initialization for SSTD algorithm
#'
#' @return  a list containing the following:
#'
#' \item{method}{character, recall of the choice of test statistics}
#'
#' \item{stats}{number, the calculated test statistics with given network list and choices}
#'
#' \item{decomp_result}{list, if \code{method = c("sum", "sum_square", "max")}, the matrix decomposition components for all pairwise differential networks are recorded; if \code{method = "tensor"}, the tensor decomposition components for the differential tensor are recorded}
#'
#' @details
#'
#' The list \code{network_list} records the pairwise differential networks \eqn{D_{AB}, D_{AH}, D_{AB}}. This package provides four options for test statistics:
#'
#' \enumerate{
#' \item{sum}, the sum of sparse leading matrix eigenvalues (sLMEs) of all pairwise differential networks:
#'
#' \deqn{Stat_sum =  \lambda(D_{AB}) + \lambda(D_{AH}) + \lambda(D_{BH}),}
#'
#' where \eqn{\lambda} refers to the sLME operation with given sparsity level set up by \code{sumabs}.
#'
#' \item{sum_square}, the sum of squared sLMEs:
#'
#' \deqn{Stat_sumsquare =  \lambda^2(D_{AB}) + \lambda^2(D_{AH}) + \lambda^2(D_{BH}).}
#'
#' \item{max}, the maximal of sLMEs:
#'
#' \deqn{Stat_max = \max(\lambda(D_{AB}), \lambda(D_{AH}), \lambda(D_{BH})).}
#'
#' \item{tensor}, the sparse leading tensor eigenvalue (sLTE) of the differential tensor:
#'
#' \deqn{Stat_tensor = \Lambda(\mathcal{D}),}
#'
#' where \eqn{\Lambda} refers to the sLTE operation with given sparsity level set up by \code{sumabs},
#' and \eqn{\mathcal{D}} is the differential tensor composed by stacking three pairwise differential networks.
#' }
#'
#' The sparse symmetric matrix decomposition is implemented by \code{symmPMD()} with parameters \code{rho, sumabs, niter, trace}.
#' The sparse symmetric tensor decomposition is implemented by \code{SSTD()}.
#' Since \code{symmPMD()} is used in \code{SSTD()}, parameters for \code{symmPMD()} are used for \code{SSTD()}.
#' While parameters \code{tensor_iter, tensor_tol, tensor_seed} should be uniquely defined for \code{tensor} method.
#'
#' @references Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2024).
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." bioRxiv, 2024-03.
#'
#' @export



net_to_stats <- function(network_list, method = c("sum", "sum_square" , "max", "tensor"),
                         rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE,
                         tensor_iter = 20, tensor_tol = 10^(-3), tensor_seed = NULL) {

  # given multiple networks, obtain the test statistics by methods
  # this function can be applied to for general types of networks, not only for covariance networks

  # rho, trace, niter, trace for symmPMD(), all applied for tensor method
  # tensor_iter, tensor_tol for SSTD()

  # if method == "sum" or method == "max", decomp_result becomes the list of sLED results for every network
  # if method == "tensor", decomp_result becomes the SSTD outputs of the tensor

  n_num <- length(network_list)

  if (method == "sum" | method == "max" | method == "sum_square") {
    decomp_result <- lapply(network_list, sLEDTestStat,
      rho = rho, sumabs.seq = sumabs, niter = niter, trace = trace
    )

    # collect eigenvalue for each network
    eigen_vec <- c()
    for (i in 1:n_num) {
      eigen_vec <- c(eigen_vec, decomp_result[[i]]$stats)
    }

    # get stats
    if (method == "sum") {
      stats <- sum(abs(eigen_vec))
    } else if (method == "max") {
      stats <- max(eigen_vec)
    } else if (method == "sum_square"){
      stats <- sum((abs(eigen_vec))^2)
    }
  } # sum & max & sum_square

  if (method == "tensor") {
    p <- dim(network_list[[1]])[1]

    T_obs <- array(0, dim = c(p, p, n_num))
    for (i in 1:n_num) { # networks to tensor
      T_obs[, , i] <- network_list[[i]]
    }

    # random initialization
    u_ini <- runif(n_num, -1, 1)
    v_ini <- runif(p, -1, 1)

    decomp_result <- SSTD_R1(T_obs, u_ini, v_ini,
      max_iter = tensor_iter, tol = tensor_tol,
      rho = rho, sumabs = sumabs, niter = niter
    )
    stats <- decomp_result$gamma_hat
  }


  return(list(method = method, stats = stats, decomp_result = decomp_result))
}

# Obtain a single test statistics from expression data ----------------

#' Generate one single snQTL test statistics from expression data
#'
#' Generate one single snQTL test statistics from a given list of expression data.
#' This function takes a list of expression data, the choice of test statistics, the choice to permute or not,
#' the choice of considering trans-correlation or not, and other computational tuning parameters as inputs.
#' Outputs include the calculated statistics, recall of the choices, and the decomposition components associated with the statistics.
#'
#' @param seed number, the random seed to shuffle the expression data if \code{permute = TRUE} and for \code{SSTD()} initialization if \code{method = "tensor"}
#' @param permute logic variable, whether to shuffle the samples in expression data; see "details"
#' @param exp_list list, a list of expression data from samples with different genotypes; see "details"
#' @param method character, the choice of test statistics (see \code{net_to_stats()})
#' @param rho number, a large positive constant adding to the diagonal elements to ensure positive definiteness in symmetric matrix spectral decomposition
#' @param sumabs number, the number specify the sparsity level in the matrix/tensor eigenvector; \code{sumabs} takes value between \eqn{1/sqrt(p)} and 1, where \eqn{p} is the dimension; \code{sumabs}\eqn{*sqrt(p)} is the upperbound of the L1 norm of the leading matrix/tensor eigenvector (see \code{symmPMD()})
#' @param niter integer, the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace logic variable, whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#' @param adj.beta number, the power transformation to the correlation matrices (see \code{getDiffMatrix()}); particularly, when \code{adj.beta=0}, the correlation matrix is used, when \code{adj.beta<0}, the covariance matrix is used.
#' @param tensor_iter integer, the maximal number of iteration in SSTD algorithm (see \code{max_iter} in \code{SSTD()})
#' @param tensor_tol number, a small positive constant for error difference to indicate the SSTD convergence (see \code{tol} in \code{SSTD()})
#' @param trans logic variable, whether to only consider the trans-correlation (between genes from two different chromosomes or regions); see "details"
#' @param location vector, the (chromosome) locations for genes if \code{trans = TRUE}
#'
#' @return  a list containing the following:
#'
#' \item{method}{character, recall of the choice of test statistics}
#'
#' \item{permute}{logic variable, recall of the choice of permutation}
#'
#' \item{stats}{number, the calculated test statistics with given expression list and choices}
#'
#' \item{decomp_result}{list, if \code{method = c("sum", "sum_square", "max")}, the matrix decomposition components for all pairwise differential networks are recorded; if \code{method = "tensor"}, the tensor decomposition components for the differential tensor are recorded}
#'
#' @details
#'
#' In \code{exp_list}, the dimensions for data matrices are n1-by-p, n2-by-p, and n3-by-p, respectively.
#' The expression data is usually normalized. We use expression data to generate the Pearson's correlation co-expression networks.
#'
#'
#' If \code{permute = TRUE}, we shuffle the samples in three expression matrices while keeping the same dimensions.
#' The test statistics from randomly shuffled data are considered as the statistics from null distribution.
#'
#' If \code{trans = TRUE}, we only consider the trans-correlation between the genes from two different chromosomes or regions in co-expression networks.
#' The entries in correlation matrices \eqn{N_{ij} = 0} if gene i and gene j are from the same chromosome or region.
#'
#' @references Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2024).
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." bioRxiv, 2024-03.
#'
#' @export


single_QTL_test_stats <- function(seed = NULL, permute = FALSE,
                                  exp_list, method = c("sum", "sum_square","max", "tensor"),
                                  rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE, adj.beta = -1,
                                  tensor_iter = 20, tensor_tol = 10^(-3),
                                  trans = FALSE, location = NULL) {

  # expression data in exp_list follows order A, B, and H, sample x genes
  # adj.beta adjust the way to get differential networks
  # if permute = TRUE, then we permute the labels for expression samples,
  # and obtain single permuted test statistics
  # recommend run permutation test with this function,
  # submit parallel jobs in which each job calculates the stats for a single permutation

  # set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # obtain network_list from exp_list
  if (!permute) { # no permute
    network_list <- vector(mode = "list", length = 3)

    if(!trans){ # include cis-correlation
      network_list[[1]] <- getDiffMatrix(exp_list[[1]], exp_list[[2]], adj.beta) # A B
      network_list[[2]] <- getDiffMatrix(exp_list[[1]], exp_list[[3]], adj.beta) # A H
      network_list[[3]] <- getDiffMatrix(exp_list[[2]], exp_list[[3]], adj.beta) # B H
    }else if(trans){ # only trans-correlation
      network_list[[1]] <- get_trans_diffmatrix(exp_list[[1]], exp_list[[2]], location, adj.beta) # A B
      network_list[[2]] <- get_trans_diffmatrix(exp_list[[1]], exp_list[[3]], location, adj.beta) # A H
      network_list[[3]] <- get_trans_diffmatrix(exp_list[[2]], exp_list[[3]], location, adj.beta) # B H
    }


  } else if (permute) {

    # collect all expression data
    Z <- rbind(exp_list[[1]], exp_list[[2]], exp_list[[3]])
    n1 <- dim(exp_list[[1]])[1]
    n2 <- dim(exp_list[[2]])[1]
    n3 <- dim(exp_list[[3]])[1]

    # shuffle labels

    i.permute <- multi_permuteIndex(n1, n2, n3)

    network_list <- vector(mode = "list", length = 3)

    if(!trans){ # include cis-correlation
      network_list[[1]] <- getDiffMatrix(Z[i.permute$i1, ], Z[i.permute$i2, ], adj.beta) # A B
      network_list[[2]] <- getDiffMatrix(Z[i.permute$i1, ], Z[i.permute$i3, ], adj.beta) # A H
      network_list[[3]] <- getDiffMatrix(Z[i.permute$i2, ], Z[i.permute$i3, ], adj.beta) # B H
    }else if(trans){ # only trans-correlation
      network_list[[1]] <- get_trans_diffmatrix(Z[i.permute$i1, ], Z[i.permute$i2, ], location, adj.beta) # A B
      network_list[[2]] <- get_trans_diffmatrix(Z[i.permute$i1, ], Z[i.permute$i3, ], location, adj.beta) # A H
      network_list[[3]] <- get_trans_diffmatrix(Z[i.permute$i2, ], Z[i.permute$i3, ], location, adj.beta) # B H
    }

  }

  # obtain stats
  res <- net_to_stats(network_list,
    method = method,
    rho = rho, sumabs = sumabs, niter = niter, trace = trace,
    tensor_iter = tensor_iter, tensor_tol = tensor_tol
  )


  return(list(method = method, permute = permute, stats = res$stats, decomp_result = res$decomp_result, network_list = network_list))
}

# Obtain empirical p-value via permutation from expression data --------

#' Spectral network quantitative trait loci (snQTL) test
#'
#' @description
#' Spectral framework to detect network QTLs affecting the co-expression networks. This is the main function for snQTL test.
#'
#' Given a list of expression data matrices from samples with different gentoypes, we test whether there are significant difference among three co-expression networks.
#' Statistically, we consider the hypothesis testing task:
#'
#' \deqn{H_0: N_A = N_B = N_H,}
#'
#' where \eqn{A,B,H} refer to different genotypes, \eqn{N} refers to the adjacency matrices corresponding to the co-expression network.
#'
#' We provide four options for the test statistics, composed by sparse matrix/tensor eigenvalues.
#' We perform permutation test to obtain the empirical p-values for the hypothesis testing.
#'
#'
#' @param exp_list list, a list of expression data from samples with different genotypes; the dimensions for data matrices are n1-by-p, n2-by-p, and n3-by-p, respectively; see "details"
#' @param method character, the choice of test statistics; see "details"
#' @param npermute number, the number of permutations to obtain empirical p-values
#' @param seeds vector, the random seeds for permutation; length of the vector is equal to the \code{npermute}
#' @param stats_seed number, the random seed for test statistics calculation with non-permuted data
#' @param rho number, a large positive constant adding to the diagonal elements to ensure positive definiteness in symmetric matrix spectral decomposition
#' @param sumabs number, the number specify the sparsity level in the matrix/tensor eigenvector; \code{sumabs} takes value between \eqn{1/sqrt(p)} and 1, where \eqn{p} is the dimension; \code{sumabs}\eqn{*sqrt(p)} is the upperbound of the L1 norm of the leading matrix/tensor eigenvector (see \code{symmPMD()})
#' @param niter integer, the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace logic variable, whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#' @param adj.beta number, the power transformation to the correlation matrices (see \code{getDiffMatrix()}); particularly, when \code{adj.beta=0}, the correlation matrix is used, when \code{adj.beta<0}, the covariance matrix is used.
#' @param tensor_iter integer, the maximal number of iteration in SSTD algorithm (see \code{max_iter} in \code{SSTD()})
#' @param tensor_tol number, a small positive constant for error difference to indicate the SSTD convergence (see \code{tol} in \code{SSTD()})
#' @param trans logic variable, whether to only consider the trans-correlation (between genes from two different chromosomes or regions); see "details"
#' @param location vector, the (chromosome) locations for genes if \code{trans = TRUE}
#'
#' @return  a list containing the following:
#'
#' \item{method}{character, recall of the choice of test statistics}
#'
#' \item{res_original}{list, test result for non-permuted data, including the recall of method choices, test statistics, and decomposition components}
#'
#' \item{res_permute}{list, test results for each permuted data, including the recall of method choices, test statistics, and decomposition components}
#'
#' \item{emp_p_value}{number, the empirical p-value from permutation test}
#'
#' @details
#'
#' In \code{exp_list}, the data matrices are usually ordered with marker's genotypes AA, BB, and AB.
#' The expression data is usually normalized. We use expression data to generate the Pearson's correlation co-expression networks.
#'
#' Given the list of co-expression networks, we generate pairwise differential networks
#' \deqn{D_{AB} = N_A - N_B, D_{AH} = N_H - N_A, D_{BH} = N_H - N_B.}
#' We use pairwise differential networks to generate the snQTL test statistics.
#'
#' We provide four options of test statistics with different choices of \code{method}:
#' \enumerate{
#' \item{sum}, the sum of sparse leading matrix eigenvalues (sLMEs) of all pairwise differential networks:
#'
#' \deqn{Stat_sum =  \lambda(D_{AB}) + \lambda(D_{AH}) + \lambda(D_{BH}),}
#'
#' where \eqn{\lambda} refers to the sLME operation with given sparsity level set up by \code{sumabs}.
#'
#' \item{sum_square}, the sum of squared sLMEs:
#'
#' \deqn{Stat_sumsquare =  \lambda^2(D_{AB}) + \lambda^2(D_{AH}) + \lambda^2(D_{BH}).}
#'
#' \item{max}, the maximal of sLMEs:
#'
#' \deqn{Stat_max = \max(\lambda(D_{AB}), \lambda(D_{AH}), \lambda(D_{BH})).}
#'
#' \item{tensor}, the sparse leading tensor eigenvalue (sLTE) of the differential tensor:
#'
#' \deqn{Stat_tensor = \Lambda(\mathcal{D}),}
#'
#' where \eqn{\Lambda} refers to the sLTE operation with given sparsity level set up by \code{sumabs},
#' and \eqn{\mathcal{D}} is the differential tensor composed by stacking three pairwise differential networks.
#' }
#'
#' Additionally, if \code{trans = TRUE}, we only consider the trans-correlation between the genes from two different chromosomes or regions in co-expression networks.
#' The entries in correlation matrices \eqn{N_{ij} = 0} if gene i and gene j are from the same chromosome or region.
#' The gene location information is required if \code{trans = TRUE}.
#'
#'
#' @references Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2024).
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." bioRxiv, 2024-03.
#'
#' @export
#'
#' @examples
#' ### artificial example
#' n1 = 50
#' n2 = 60
#' n3 = 100
#'
#' p = 200
#'
#' location = c(rep(1,20), rep(2, 50), rep(3, 100), rep(4, 30))
#'
#' ## expression data from null
#' set.seed(0416) # random seeds for example data
#' exp1 = matrix(rnorm(n1*p, mean = 0, sd = 1), nrow = n1)
#' exp2 = matrix(rnorm(n2*p, mean = 0, sd = 1), nrow = n2)
#' exp3 = matrix(rnorm(n3*p, mean = 0, sd = 1), nrow = n3)
#'
#' exp_list = list(exp1, exp2, exp3)
#'
#' result = network_QTL_test(exp_list = exp_list, method = 'tensor',
#'                           npermute = 100, seeds = 1:100, stats_seed = 0416,
#'                           trans = TRUE, location = location)
#'
#' result$emp_p_value
#'
#' ## expression data from alternative
#' Sigma = diag(p)
#' # trans-correlation
#' Sigma[20:50, 20:50] = Sigma[20:50, 20:50] + 0.5
#'
#' set.seed(0416) # random seeds for example data
#' exp1 = matrix(rnorm(n1*p, mean = 0, sd = 1), nrow = n1)
#' exp2 = matrix(rnorm(n2*p, mean = 0, sd = 1), nrow = n2)
#' exp3 = MASS::mvrnorm(n3, mu = rep(0,p), Sigma = Sigma)
#'
#' exp_list = list(exp1, exp2, exp3)
#'
#' result = network_QTL_test(exp_list = exp_list, method = 'tensor',
#'                           npermute = 100, seeds = 1:100, stats_seed = 0416,
#'                           trans = TRUE, location = location)
#'
#' result$emp_p_value
#'
#' # leverage
#' leverage = (result$res_original$decomp_result$v_hat)^2
#'


network_QTL_test <- function(exp_list, method = c("sum", "sum_square", "max", "tensor"),
                             npermute = 100, seeds = 1:100, stats_seed = NULL,
                             rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE, adj.beta = -1,
                             tensor_iter = 20, tensor_tol = 10^(-3),
                             trans = FALSE, location = NULL) {

  # main network QTL test function
  # no parallel version
  # do not recommend run this function with large expression data, large number of permutations on local PC
  # expression data in exp_list follows order A, B, and H

  # warning
  if(npermute != length(seeds)){
    cat("number of seeds should be equal to the number of permutations")
    return()
  }

  # original stats
  res_original <- single_QTL_test_stats(
    seed = stats_seed, permute = FALSE,
    exp_list = exp_list, method = method,
    rho = rho, sumabs = sumabs, niter = niter, trace = trace, adj.beta = adj.beta,
    tensor_iter = tensor_iter, tensor_tol = tensor_tol,
    trans, location
  )
  # permutation results
  res_permute <- lapply(seeds, single_QTL_test_stats,
    permute = TRUE, exp_list = exp_list, method = method,
    rho = rho, sumabs = sumabs, niter = niter, trace = trace, adj.beta = adj.beta,
    tensor_iter = tensor_iter, tensor_tol = tensor_tol,
    trans, location
  )

  # calculate empirical p value
  stats_obs <- res_original$stats
  stats_permute <- rep(0, npermute)
  for (i in 1:npermute) {
    stats_permute[i] <- res_permute[[i]]$stats
  }
  emp_p_value <- sum(stats_permute > stats_obs) / npermute

  return(list(method = method, res_original = res_original, res_permute = res_permute, emp_p_value = emp_p_value))
}
