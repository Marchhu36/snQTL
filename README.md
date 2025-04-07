# snQTL: A spectral framework to map QTLs affecting joint differential networks of gene co-expression

This is the R software for the paper  
> Hu, J., Weber, J. N., Fuess, L. E., Steinel, N. C., Bolnick, D. I., & Wang, M. (2025). A spectral framework to map QTLs affecting joint differential networks of gene co-expression. PLOS Computational Biology. [bioRxiv version](https://www.biorxiv.org/content/10.1101/2024.03.29.587398v1).

## Short introduction to snQTL

snQTL is a powerful statistical methods to analyze genotype $\rightarrow$ network $\rightarrow$ phenotype mechanism. snQTL aims to map quantitative trait loci (QTL) affecting the gene co-expression networks. 

Given gene expression read counts and the genotypes of genetic markers from the same set of samples, we want to test whether the gene co-expression networks change significantly across different genotypes. Statistically, snQTL tackles the hypothesis testing problem:

$$H_0: N_A = N_B = N_H,$$
where $A,B,H$ refer to different genotypes, and $N_A, N_B, N_H$ refer to the gene co-expression networks corresponding different genotypes. 

We propose spectral test statistics for above test, based on tensor/matrix decomposition. Permutation procedure is used to obtain the empirical p-value. 

snQTL takes a list of expression data divided by genotype as input. snQTL outputs the empirical p-value for above test and associated decomposition components, which indicate the gene contribution to the co-expression network changes and compose the joint differential network. 

## Installation

This package can be installed through `devtools` in R:

```
install.packages("devtools") ## if not installed
library("devtools")
install_github("Marchhu36/snQTL")
```

### Main functions

## Demo

We use synthetic data to illustrate the usage of the main function `network_QTL_test()`. Consider the following setup. 

```
library(snQTL)

# sample size for different genotypes
n1 = 50 # AA/A
n2 = 60 # BB/B
n3 = 100 # AB/H

# number of genes
p = 200 

# location information for genes
location = c(rep(1,20), rep(2, 50), rep(3, 100), rep(4, 30))
```
**Note**

- In current software, main function generates the co-expression networks from expression read counts through Pearson's correlation. Though, function `net_to_stats()` still works to calculate snQTL test statistics with any given pairwise differential networks generated from other methods. More options for co-expression networks constructions are under development. 

- The input expression read counts are usually normalized and pre-processed with quality control. We skip this part in this demo for simplicity. Our package does not contain the pre-processing functions.

- We focus on the test with trans-correlation in this demo. The correlations between two genes from the same location/chromosome are set to 0 when `trans = TRUE`. 

### Data from null hypothesis

We fist try snQTL when the synthetic data follows the null hypothesis. 

```
set.seed(0416) # random seeds for example data
exp1 = matrix(rnorm(n1*p, mean = 0, sd = 1), nrow = n1)
exp2 = matrix(rnorm(n2*p, mean = 0, sd = 1), nrow = n2)
exp3 = matrix(rnorm(n3*p, mean = 0, sd = 1), nrow = n3)

exp_list = list(exp1, exp2, exp3)
```
Then, we apply snQTL with tensor statistics and 100 permutations and leave all other arguments as default (except random seed). 

```
result = network_QTL_test(exp_list = exp_list, method = 'tensor', 
                          npermute = 100, seeds = 1:100, stats_seed = 0416,
                          trans = TRUE, location = location)
```
Check empirical p-value.
```
result$emp_p_value
## [1] 0.49
```
The p-value is quite large, which agrees with the expectation. To get a more accurate p-value, increase the number of permutations as needed.

### Data from alternative hypothesis

Now, we suppose the heterozygous genotype, H, has a special effect to the co-expression network. We add extra correlation signal when generating heterozygous expression read counts. Note that such extra signal exists between the genes from different locations. The correlation signal within the same location will be excluded when we only consider the trans-correlation (`trans = TRUE`).

```
Sigma = diag(p)
# trans-correlation
Sigma[20:50, 20:50] = Sigma[20:50, 20:50] + 0.5

set.seed(0416) # random seeds for example data
exp1 = matrix(rnorm(n1*p, mean = 0, sd = 1), nrow = n1)
exp2 = matrix(rnorm(n2*p, mean = 0, sd = 1), nrow = n2)
exp3 = MASS::mvrnorm(n3, mu = rep(0,p), Sigma = Sigma)

exp_list = list(exp1, exp2, exp3)
```
Then, we apply snQTL with the same parameter setup. 
```
result = network_QTL_test(exp_list = exp_list, method = 'tensor', 
                          npermute = 100, seeds = 1:100, stats_seed = 0416,
                          trans = TRUE, location = location)
```
Check empirical p-value.
```
result$emp_p_value
## [1] 0.01
```
The p-value is quite small, which indicates there are significant signals among the networks. This example also verifies that snQTL works for non-linear network effect: in this case, only heterozygous genotype leads to a network change while the other two homozygous genotype share similar co-expression networks. 

### Leverage, loading, and joint differential network

When using tensor statistics `method = "tensor"`, we are able to obtain the (1) gene leverage indicating the contribution of each gene to the network changes, (2) the joint differential network, and (3) the loading revealing the weights of pairwise comparisons towards the joint differential network.

We can find leverage, loading, and adjacency matrix for the joint differential network as
```
leverage = (result$res_original$decomp_result$v_hat)^2
loading = result$res_original$decomp_result$u_hat
joint_diff_network = result$res_original$decomp_result$v_hat %*% t(result$res_original$decomp_result$v_hat)
```

## Generalization

snQTL is a flexible framework to test multiple biological networks. Here, we provide a generalization when samples are separated into $K$ groups based on other interested biological factors (e.g., location, treatment) other than genotype. Then, snQTL tackles the hypothesis testing question:

$$H_0: N_1 = ... = N_K,$$

where $N_k$ refers to the network corresponding to the group $k$. 

The main change compared with the original problem is that we need to consider $K(K-1)/2$ pairwise differential networks. Specifically, we calculate a list of pairwise differential networks $D^{(k,l)} = N_l - N_k$ for all $1 \leq k < l \leq K$. Thus differential tensor $\mathcal{D}$ has dimension $p \times p \time q$, where $q = K(K-1)/2$. We calculate the test statistics based on the matrix spectral statistics for all $D^{(k,l)}$'s or the tensor spectral statistics of $\mathcal{D}$. 

Below, we provide a demo for the generalized case. 
