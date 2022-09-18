
# Set seed for reproducible results
set.seed(42)

# Create 400 covariate data x which lie between (-1, 1)
N <- 20000
x <- seq(-1, 1, length.out = N)

sub_n=20

# Create n x D design matrix
D <- 2
# We learn a linear function
X <- matrix(c(rep(1, N), x), ncol = D)

# True values of regression coeffiecients theta
true_theta <- c(-2, 5)


tphi_1 = 1
tphi_2 = 5
# use mis-specified models for split variance, expect coreset to be performing better
# try and make discrepancies in signals, try to create heteroscadekstic noise 
tphi=0.5


#####corruption#######
y_o = X%*%true_theta + sqrt(tphi)*rnorm(N)
corruption = sample(N,50)
y_o[corruption] = y_o[corruption]+sqrt(tphi_2)*rnorm(50)


# Obtain the vector with probabilities of success p using the probit link
#p <- pnorm(X %*% true_theta+0.5*rnorm(N))
p <- pnorm(y_o)

# Generate binary observation data y
y <- rbinom(N, 1, p)

# Variables that we will need later
N1  <- sum(y)  # Number of successes
N0  <- N - N1  # Number of failures

# Fit the model to the data
fit <- glm(y ~ x, family = binomial(link = probit))

# MLE estimates of the regression coefficients
# (Intercept)           x 
#  -0.5490056   3.1296282
mle_theta <- fit$coefficients

# Library for sampling from Multivariate Normal distribution
require(mvtnorm)
# Library for sampling from Truncated Normal distribution
require(truncnorm)
require(pbivnorm)




#Probit-log-likelihood
Probit_Log_likelihood_n = function(b, x, y) {
  return(
    y*log(pnorm(x %*% b))
    + (1-y)*(1 -log(1 - pnorm(x %*% b))) 
  )
}

# Conjugate prior on the coefficients \theta ~ N(theta_0, Q_0)
theta_0 <- rep(0, D)
Q_0 <- diag(10, D)

J = 100

b_0 = rmvnorm(J,theta_0,Q_0)


# Conjugate prior on the coefficients \theta ~ N(theta_0, Q_0)
theta_0 <- rep(0, D)
Q_0 <- diag(10, D)

# Initialize parameters
theta <- rep(0, D)


# Number of simulations for Gibbs sampler
N_sim <- 10000 
# Burn in period
burn_in <- 5000
# Matrix storing samples of the \theta parameter
theta_chain <- matrix(0, nrow = N_sim, ncol = D)


# ---------------------------------
# Gibbs sampling algorithm
# ---------------------------------

# Compute posterior variance of theta
prec_0 <- solve(Q_0)


for (t in 2:N_sim) {
  
  
  
  
  
  sub_sample = sample(seq(1,N),sub_n)
  

  
  y_tilde = y[sub_sample]
  sub_sample_n = length(y_tilde)
  X_tilde = X[sub_sample,]
  z <- rep(0, sub_sample_n)

  
  V <- solve(prec_0 + crossprod(X_tilde, X_tilde))

  
  M1  <- sum(y_tilde)  # Number of successes
  M0  <- sub_sample_n - M1  # Number of failures
  
  
  
  
  # Update Mean of z
  mu_z <- X_tilde %*% theta
  # Draw latent variable z from its full conditional: z | \theta, y, X
  
  if(M1==0){
    z[y_tilde == 0] <- rtruncnorm(M0, mean = mu_z[y_tilde == 0], sd = 1, a = -Inf, b = 0)
  }
  else if (M0==0){
    z[y_tilde == 1] <- rtruncnorm(M1, mean = mu_z[y_tilde == 1], sd = 1, a = 0, b = Inf)
  }
  else{
    z[y_tilde == 0] <- rtruncnorm(M0, mean = mu_z[y_tilde == 0], sd = 1, a = -Inf, b = 0)
    z[y_tilde == 1] <- rtruncnorm(M1, mean = mu_z[y_tilde == 1], sd = 1, a = 0, b = Inf)
  }
  
  
  # Compute posterior mean of theta
  M <- V %*% (prec_0 %*% theta_0 + crossprod(X_tilde, z))
  # Draw variable \theta from its full conditional: \theta | z, X
  theta <- c(rmvnorm(1, M, V))
  
  # Store the \theta draws
  theta_chain[t, ] <- theta
}

# ---------------------------
# Get posterior mean of \theta
# ---------------------------
# (Intercept)           x 
#  -1.436119   3.631905
post_theta <- colMeans(theta_chain[-(1:burn_in), ])

# Plot covariates x versus observations y
plot(x, y, main = "Synthetic data")
# Show the fitted function using the posterior mean estimates
lines(x = x, y = pnorm(X %*% mle_theta), col = "red3", lwd = 2)
lines(x = x, y = pnorm(X %*% post_theta), col = "blue3", lwd = 2)
legend("bottomright", legend=c("MLE","Post. Mean"), col=c("red3","blue3"), 
       bty = 'n', lwd = 2, inset = c(0.02, 0.08), lty = 1, cex = 0.9)

