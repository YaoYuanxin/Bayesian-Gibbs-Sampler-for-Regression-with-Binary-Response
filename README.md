# Bayesian-Regression-Gibbs-Sampers
Gibbs Samplers for Linear Regression and Logistic/Probit Regression.

This my implementation of estimating regression coefficients in a Logistic/Probit Regression problem. 

Bernoulli response (i.e values of y is 0 and 1, often occuring in classification problems) in a regression problem does not have a mathematical closed form 
solution. 
From the perspective of Bayesian statistics, a Gibbs sampler with generates posterior solutions that are comparable to other state-of the art methods.

Additionally, in each step of the Gibbs sampler, the full conditional distribution is generated using a small, randomly selected subset of the data, rather than the
larger full data. This sub-sampling step speeds up the computation while preserving the quality of the inferential resuls.
