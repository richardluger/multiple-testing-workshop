# ---------------------------------------------------------------
# Opener: the multiple testing problem inside a linear regression
# ---------------------------------------------------------------
# We simulate M datasets under a pure-null classical linear regression model:
#     y_i = eps_i,   eps_i ~ N(0,1),
# with a design matrix X of k independent N(0,1) regressors that
# are unrelated to y.  All slope coefficients are therefore zero.
#
# In each dataset we compute two things at the 5% level:
#   (1) reject_any: does AT LEAST ONE of the k individual t-tests
#                   reject H0: beta_j = 0?
#   (2) reject_F:   does the joint F-test of
#                   H0: beta_1 = ... = beta_k = 0 reject?
#
# Under correct size, both Monte Carlo rejection rates should be
# close to alpha = 0.05.  We will see that only the F-test is.
# ---------------------------------------------------------------

rm(list=ls())

M <- 10000    	                                # number of Monte Carlo replications

n <- 50         	                            # sample size per dataset
alpha <- 0.05       	                        # nominal level for every test

reject_any <- reject_F <- logical(M)			# Storage for the two rejection indicators



k <- 1    # <---------------------------------- number of regressors (all null) : try k =1, 2, 3, ...



for (m in seq_len(M)) {
  
  # Generate a dataset under the null 
  X <- matrix(rnorm(n * k), n, k)             # n x k design, iid N(0,1)
  y <- rnorm(n)                               # response independent of X
  
  # Fit the linear regression 
  fit <- summary(lm(y ~ X))                   # OLS with intercept + k slopes
  

  # (1) Naive "any significant t" rule 
  # Extract the k slope t-statistics (drop the intercept row)
  t_stats <- fit$coefficients[-1, "t value"]

  # Two-sided critical value with n - k - 1 residual degrees of freedom
  crit    <- qt(1 - alpha/2, df = n - k - 1)

  # Flag the replication if any |t_j| exceeds the critical value
  reject_any[m] <- any(abs(t_stats) > crit)
  

  # (2) Joint F-test of all slopes = 0 
  f <- fit$fstatistic                         # c(value, numdf, dendf)

  # p-value from the F distribution; reject if below alpha
  reject_F[m] <- pf(f[1], f[2], f[3], lower.tail = FALSE) < alpha

}


# Empirical rejection rates under the null:
c(any_t = mean(reject_any), F = mean(reject_F))











