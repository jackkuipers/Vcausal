# Vcausal
Code to compute the variance of causal estimators for binary v-structures to accompany the article "The variance of causal effect estimators for binary v-structures", Kuipers and Moffa, Journal of Causal Inference, 10, 90-105 (2022) [https://doi.org/10.1515/jci-2021-0025 and arXiv version at https://arxiv.org/abs/2004.09181]

The R file VcausalMC.R runs Monte Carlo simulations, while Vcausalfns.R and the Maple file contain functions to compute the numerical result.

The R code now also computes the covariance of the different causal estimators.

Also included are the R files VcausalMC_fixedX.R and Vcausalfns_fixedX.R which run Monte Carlo simulations and compute the numerical results when X is fixed by block randomisation and the number of cases is predefined.
