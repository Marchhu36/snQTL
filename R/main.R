##### snQTL (Generalized) -- main functions and bricks
##### Jiaxin Hu, Apr 5, 2025
##### Generalized snQTL for multiple network comparison

# Note:
# Current version only supports the empirical p-value calculation for correlation networks,
# since the permutation step shuffles the sample labels and need to RE-GENERATE networks
# to get a new statistic.

# People may extend this test once the network generation machanism is known.

# # Dependencies (Just for construction)
# source("symmPMD.R")
# source("SSTD.R")
# source("tensor_class.R")
# source("tensor_product.R")

############################# Main functions ##################################

# Key function 1: Differential networks to statistics -------------------------

#' Test statistics for snQTL
#'
#' @description
#' Generate snQTL test statistics from a given list of differential networks.
#' This function takes a list of differential networks, the choice of test statistics, and other computational tuning parameters as inputs.
#' Outputs include the calculated statistics, recall of the choice, and the decomposition components associated with the statistics.
#'
#' @param diffnet_list list, a list of p-by-p differential networks
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
#' The list \code{diffnet_list} records the pairwise differential networks \eqn{D_{AB}, D_{AH}, D_{AB}}. This package provides four options for test statistics:
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
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." PLOS Computational Biology.
#'
#' @export

diffnet_to_snQTL_stats = function(diffnet_list,
                                  method = c("sum", "sum_square" , "max", "tensor"),
                                  rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE,
                                  tensor_iter = 20, tensor_tol = 10^(-3), tensor_seed = NULL){

  # given multiple networks, obtain the test statistics by methods
  # this function can be applied to for general types of networks, not only for covariance networks

  # rho, trace, niter, trace for symmPMD(), all applied for tensor method
  # tensor_iter, tensor_tol for SSTD()

  # if method == "sum" or method == "max", decomp_result becomes the list of sLED results for every network
  # if method == "tensor", decomp_result becomes the SSTD outputs of the tensor


  n_num <- length(diffnet_list)

  if (method == "sum" | method == "max" | method == "sum_square") {
    decomp_result <- lapply(diffnet_list, sLME,
                            rho = rho, sumabs.seq = sumabs, niter = niter, trace = trace)

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
    p <- dim(diffnet_list[[1]])[1]

    T_obs <- array(0, dim = c(p, p, n_num))
    for (i in 1:n_num) { # networks to tensor
      T_obs[, , i] <- diffnet_list[[i]]
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






# Key function 2 (Correlation-based): Single test statistics calculation from expression list -------------------------

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
#' If \code{permute = TRUE}, we shuffle the samples in three expression matrices while keeping the same dimensions.
#' The test statistics from randomly shuffled data are considered as the statistics from null distribution.
#'
#' If \code{trans = TRUE}, we only consider the trans-correlation between the genes from two different chromosomes or regions in co-expression networks.
#' The entries in correlation matrices \eqn{N_{ij} = 0} if gene i and gene j are from the same chromosome or region.
#'
#' @references Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2024).
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." PLOS Computational Biology.
#'
#' @export


single_exp_to_snQTL_stats = function(seed = NULL, permute = FALSE,
                                     exp_list, method = c("sum", "sum_square","max", "tensor"),
                                     rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE, adj.beta = -1,
                                     tensor_iter = 20, tensor_tol = 10^(-3),
                                     trans = FALSE, location = NULL){

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

  if(!permute){ # no permute

    diffnet_list = get_diffnet_list_from_exp(exp_list, adj.beta, trans, location)

  }else if (permute){ # permute

    k = length(exp_list)
    Z = c()
    n_vec = rep(NA,k)
    for (i in 1:k) {
      Z = rbind(Z, exp_list[[i]])
      n_vec[i] = dim(exp_list[[i]])[1]
    }

    n_list = multi_permuteIndex(n_vec)

    new_exp_list = vector(mode = "list", length = k)
    for (i in 1:k) {
      new_exp_list[[i]] = Z[n_list[[i]], ]
    }

    diffnet_list = get_diffnet_list_from_exp(new_exp_list, adj.beta, trans, location)

  }# end diffnet_list

  # obtain stats
  res <- diffnet_to_snQTL_stats(diffnet_list,
                      method = method,
                      rho = rho, sumabs = sumabs, niter = niter, trace = trace,
                      tensor_iter = tensor_iter, tensor_tol = tensor_tol)

  return(list(method = method, permute = permute, stats = res$stats, decomp_result = res$decomp_result, diffnet_list = diffnet_list))
}

# Key function 3 (Correlation-based): snQTL test for correlation networks from expression list -------------------------

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
#' NOTE: This function is also applicable for generalized cases to compare multiple (K > 3) biological networks.
#' Instead of separating the samples by genotypes, people can separate the samples into K groups based on other interested metrics, e.g., locations, treatments.
#' The generalized hypothesis testing problem becomes
#' \deqn{H_0: N_1 = ... = N_K,}
#' where \eqn{N_k} refers to the correlation-based network corresponding to the group k.
#' For consistency, we stick with the original genotype-based setting in this help document.
#' See details and examples for the generalization on the Github manual https://github.com/Marchhu36/snQTL.
#'
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
#' "A spectral framework to map QTLs affecting joint differential networks of gene co-expression." PLOS Computational Biology.
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
#' # loading
#' loading = result$res_original$decomp_result$u_hat
#'# joint differential network
#'joint_diff_network = result$res_original$decomp_result$v_hat %*% t(result$res_original$decomp_result$v_hat)


snQTL_test_corrnet = function(exp_list, method = c("sum", "sum_square", "max", "tensor"),
                              npermute = 100, seeds = 1:100, stats_seed = NULL,
                              rho = 1000, sumabs = 0.2, niter = 20, trace = FALSE, adj.beta = -1,
                              tensor_iter = 20, tensor_tol = 10^(-3),
                              trans = FALSE, location = NULL){

  # main network QTL test function
  # no parallel version
  # do not recommend run this function with large expression data, large number of permutations on local PC

  # warning
  if(npermute != length(seeds)){
    cat("number of seeds should be equal to the number of permutations")
    return()
  }

  # original stats
  res_original <- single_exp_to_snQTL_stats(
    seed = stats_seed, permute = FALSE,
    exp_list = exp_list, method = method,
    rho = rho, sumabs = sumabs, niter = niter, trace = trace, adj.beta = adj.beta,
    tensor_iter = tensor_iter, tensor_tol = tensor_tol,
    trans, location
  )
  # permutation results
  res_permute <- lapply(seeds, single_exp_to_snQTL_stats,
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

############################# Bricks ##########################################

# Brick function 1: sLME for matrices (from sLED test) --------------------------
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
#'
#' @export

sLME = function(Dmat, rho=1000, sumabs.seq=0.2,
                niter=20, trace=FALSE){
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

# Brick function 2: get a list of differential matrices between multiple networks ------------------

# get difference matrix ----------

#' Get the list of differential matrix from a list of expression data
#'
#' @description Given a list of expression data, \eqn{X_1, ..., X_K}, compute the list of differential matrix
#' \deqn{D^{(k,l)} = N(X_l) - N(X_k), k < l, }
#' where N() is the covariance matrix, or the weighted adjacency matrices defined as
#' \deqn{N_{ij} = |corr(i, j)|^beta}
#' for some constant beta > 0, 1 <= i, j <= p.
#' Let N represent the regular correlation matrix when beta=0, and covariance matrix when beta<0.
#' In total, we will have K*(K-1)/2 pairwise differential networks in the list.
#'
#'
#' If \code{trans = TRUE}, we let \eqn{N_{ij} = 0} if \eqn{i, j} are from the same region based on \code{location}.
#' Under gene co-expression context, trans-correlation usually refer to the correlation between
#' two genes \eqn{i, j} from two chromosomes.
#'
#' @param exp_list a list of nk-by-p matrices from the K populations.
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta Power to transform correlation matrices to weighted adjacency matrices
#'        by \eqn{N_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson's correlation.
#'        When \code{adj.beta=0}, the correlation marix is used.
#'        When \code{adj.beta<0}, the covariance matrix is used.
#'        The default value is \code{adj.beta=-1}.
#' @param trans logic variable, whether to only consider the trans-correlation (between genes from two different chromosomes or regions)
#' @param location vector, the (chromosome) locations for items
#'
#'
#'
#' @return A list of p-by-p differential matrix \eqn{D^{(k,l)}, k < l}.
#'
#' @export


get_diffnet_list_from_exp = function(exp_list, adj.beta = -1, trans = FALSE, location = NULL){

  # net_list, list of networks
  k = length(exp_list)

  # collect all k(k-1)/2 pairwise diffnet
  diffnet_list = vector(mode = "list", length = k*(k-1)/2)
  name_vec = rep(NA, k*(k-1)/2)
  count = 1
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      diffnet_list[[count]] = get_diffnet_from_exp(exp_list[[j]], exp_list[[i]], adj.beta, trans, location)
      name_vec[count] = paste0(i,"-",j)
      count = count + 1
    }
  }
  names(diffnet_list) = name_vec

  return(diffnet_list)
}

# Brick function 2: get pairwise correlation differential network from expression ------------------
# get difference matrix ----------

#' The differential matrix
#'
#' @description Given observations from two populations X and Y,
#' compute the differential matrix
#' \deqn{D = N(Y) - N(X)}
#' where N() is the covariance matrix, or the weighted adjacency matrices defined as
#' \deqn{N_{ij} = |corr(i, j)|^beta}
#' for some constant beta > 0, 1 <= i, j <= p.
#' Let N represent the regular correlation matrix when beta=0, and covariance matrix when beta<0.
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
#'@param trans logic variable, whether to only consider the trans-correlation (between genes from two different chromosomes or regions); see "details"
#' @param location vector, the (chromosome) locations for items
#'
#' @return The p-by-p differential matrix \eqn{D = N(Y) - N(X)}.
#'
#' @export

get_diffnet_from_exp = function(X, Y, adj.beta = -1, trans = FALSE, location = NULL){

  if (adj.beta < 0) {
    Dmat <- cov(Y) - cov(X)
  } else if (adj.beta == 0) {
    Dmat <- cor(Y) - cor(X)
  } else {
    Dmat <- abs(cor(Y))^adj.beta - abs(cor(X))^adj.beta
  }

  if(trans){
    # location: the locations of each gene, only positive numbers refer to a valid location
    locs = unique(location)[unique(location) >= 0]
    for (i in locs) {
      ind_vec = location == i
      # delete the entries when two loci within the same location
      Dmat[ind_vec, ind_vec] = 0
    }
  }

  return(Dmat)
}

# Tools --------------------------------

normalize <- function(x) {
  return(x / sqrt(sum(x^2)))
}

multi_permuteIndex <- function(n_vec){

  # n_vec = c(n1, n2, ..., nk)

  n_sum = sum(n_vec)
  k = length(n_vec)
  i.sample <- sample(1:n_sum, replace=FALSE)

  start_tmp = 1
  n_list = vector(mode = "list", length = k)
  for (j in 1:k) {
    n_list[[j]] = i.sample[start_tmp:(start_tmp+n_vec[j]-1)]
    start_tmp = start_tmp+n_vec[j]
  }

 return(n_list)
}

