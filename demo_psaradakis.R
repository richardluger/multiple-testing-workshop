rm(list = ls())

# ============================================================
#  demo_psaradakis.R
#  Replication of Psaradakis (2001), Table 2 -- United States
#  "p-Value Adjustments for Multiple Tests for Nonlinearity"
#  Studies in Nonlinear Dynamics & Econometrics 4(3): 95-100
# ------------------------------------------------------------
#  Self-contained: US GDP growth rates are embedded directly.
#  Source: FRED GDPC1 (Real GDP, billions of chained 2012 USD,
#          seasonally adjusted quarterly), downloaded 2025.
#          Log first differences, 1959 Q2 -- 1997 Q2, n = 153.
#
#  Note: Psaradakis (2001) uses Hochberg (1988) as the
#  step-up alternative to Bonferroni. This demo uses Holm
#  (step-down) instead, which is covered in the workshop
#  slides. Results for the US are numerically identical
#  because McLEOD is the only test near significance.
#
#  Only external dependency: tseries (for bds.test).
#  install.packages("tseries")
# ============================================================

library(tseries)   # bds.test()
set.seed(2001)     # fixes random weights in V23 and NEURAL2


# ============================================================
# 1.  Data: log first differences of US Real GDP
#     1959 Q2 to 1997 Q2,  n = 153 observations
# ============================================================

x <- c(
  0.02228419, 0.00069702, 0.00284575, 0.02223718, -0.00539812, 0.00488735, -0.01291434, 0.00672750,
  0.01683602, 0.01902159, 0.01942677, 0.01767943, 0.00900686, 0.01222069, 0.00328861, 0.01086170,
  0.01116162, 0.02175005, 0.00653763, 0.02086788, 0.01083181, 0.01550519, 0.00308714, 0.02391092,
  0.01255556, 0.02198809, 0.02277763, 0.02404690, 0.00340889, 0.00843267, 0.00817056, 0.00881893,
  0.00061333, 0.00941629, 0.00751026, 0.02018756, 0.01656576, 0.00771535, 0.00392005, 0.01552320,
  0.00302895, 0.00657784, -0.00489405, -0.00149003, 0.00141618, 0.00917248, -0.01077356, 0.02679917,
  0.00539360, 0.00819083, 0.00234589, 0.01820543, 0.02243922, 0.00939661, 0.01660572, 0.02444811,
  0.01082859, -0.00527226, 0.00944754, -0.00863422, 0.00237415, -0.00949773, -0.00389162, -0.01225108,
  0.00712159, 0.01697554, 0.01338100, 0.02224540, 0.00730698, 0.00545915, 0.00720005, 0.01179253,
  0.01923309, 0.01787101, 0.00001993, 0.00319332, 0.03791960, 0.01000725, 0.01335172, 0.00179433,
  0.00106688, 0.00740059, 0.00249733, 0.00314077, -0.02081958, -0.00118925, 0.01847747, 0.01940473,
 -0.00743835, 0.01190389, -0.01095576, -0.01565840, 0.00455140, -0.00382999, 0.00040001, 0.01309183,
  0.02250016, 0.01979238, 0.02064720, 0.01935930, 0.01713006, 0.00959393, 0.00817423, 0.00964315,
  0.00876620, 0.01515726, 0.00740699, 0.00929415, 0.00449286, 0.00952121, 0.00534959, 0.00739107,
  0.01072407, 0.00863440, 0.01702384, 0.00515551, 0.01305474, 0.00584260, 0.01323657, 0.01011264,
  0.00760203, 0.00738006, 0.00196787, 0.01086953, 0.00362313, 0.00066560, -0.00914574, -0.00469039,
  0.00776665, 0.00504086, 0.00347939, 0.01190142, 0.01078467, 0.00983260, 0.01037455, 0.00166805,
  0.00580483, 0.00476020, 0.01350982, 0.00965685, 0.01346082, 0.00582882, 0.01139058, 0.00354153,
  0.00297879, 0.00847127, 0.00676787, 0.00746259, 0.01654187, 0.00893034, 0.01033016, 0.00643369,
  0.01651505
)

n_tot <- length(x)
cat(sprintf("US GDP growth: n = %d obs  (1959 Q2 -- 1997 Q2)\n\n", n_tot))


# ============================================================
# 2.  Lag selection by AICC  (Hurvich & Tsai 1989)
#     AICc = AIC + 2k(k+1)/(n-k-1),  k = p + 2 params
# ============================================================

aicc_ar <- function(x, p) {
  fit <- arima(x, order = c(p, 0, 0), method = "ML")
  k   <- p + 2                              # AR coefs + mean + variance
  -2 * fit$loglik + 2*k + 2*k*(k+1) / (length(x) - k - 1)
}

aicc_vals <- sapply(1:8, aicc_ar, x = x)
m         <- which.min(aicc_vals)

{
cat(sprintf("AICC-selected lag order:  m = %d\n", m))
cat(sprintf("(Psaradakis obtains m = 2 for the United States)\n\n"))
}


# ============================================================
# 3.  OLS design matrix, coefficients, residuals
# ============================================================

Xm  <- embed(x, m + 1)
y_  <- Xm[, 1]
Z_  <- cbind(1, Xm[, -1])
N   <- nrow(Z_)
p_  <- ncol(Z_)

b_  <- solve(crossprod(Z_), crossprod(Z_, y_))
e_  <- as.numeric(y_ - Z_ %*% b_)
yh_ <- as.numeric(Z_ %*% b_)

cat(sprintf("Effective OLS sample:  N = %d\n\n", N))


# ============================================================
# 4.  The seven nonlinearity tests
# ============================================================

keenan_test <- function(y, Z, e, yhat) {
  N <- length(e);  p <- ncol(Z)
  bv <- solve(crossprod(Z), crossprod(Z, yhat^2))
  v  <- yhat^2 - Z %*% bv
  eta    <- sum(e * v) / sum(v^2)
  ss_reg <- eta^2 * sum(v^2)
  ss_res <- sum(e^2) - ss_reg
  df_den <- N - 2 * p
  pf(ss_reg / (ss_res / df_den), 1, df_den, lower.tail = FALSE)
}

tsay_test <- function(y, Z) {
  N <- nrow(Z);  p <- ncol(Z);  m <- p - 1
  lags <- Z[, -1, drop = FALSE]
  idx <- which(upper.tri(matrix(0, m, m), diag = TRUE), arr.ind = TRUE)
  W   <- matrix(NA, N, nrow(idx))
  for (j in seq_len(nrow(idx)))
    W[, j] <- lags[, idx[j, 1]] * lags[, idx[j, 2]]
  ord <- order(lags[, 1])
  y_o <- y[ord];  Z_o <- Z[ord, ];  W_o <- W[ord, ]
  n0 <- 2 * p
  ep <- numeric(N - n0)
  for (t in seq(n0 + 1L, N)) {
    sl  <- seq_len(t - 1L)
    Zt  <- Z_o[sl, , drop = FALSE]
    bt  <- solve(crossprod(Zt) + diag(1e-9, p), crossprod(Zt, y_o[sl]))
    ep[t - n0] <- y_o[t] - Z_o[t, ] %*% bt
  }
  W_t  <- W_o[seq(n0 + 1L, N), , drop = FALSE]
  lm_t <- lm(ep ~ W_t)
  fs   <- summary(lm_t)$fstatistic
  unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
}

reset_test <- function(y, Z, e, yhat) {
  N <- length(e)
  Z_aug <- cbind(Z, yhat^2, yhat^3, yhat^4)
  e_aug <- residuals(lm(y ~ Z_aug - 1))
  ss0   <- sum(e^2);  ss1 <- sum(e_aug^2)
  df1   <- 3;  df2 <- N - ncol(Z_aug)
  pf((ss0 - ss1) / df1 / (ss1 / df2), df1, df2, lower.tail = FALSE)
}

mcleod_test <- function(e)
  Box.test(e^2, lag = 20, type = "Ljung-Box")$p.value

bds_test2 <- function(e)
  suppressWarnings(bds.test(e, m = 2, eps = sd(e))$p.value[1])

v23_test <- function(y, Z, e, q = 3) {
  lags <- Z[, -1, drop = FALSE]
  Gam  <- matrix(rnorm(q * ncol(lags)), q, ncol(lags))
  Gam  <- Gam / sqrt(rowSums(Gam^2))
  Psi  <- apply(Gam, 1, function(g) plogis(lags %*% g))
  lm_v <- lm(e ~ cbind(Z, Psi) - 1)
  pchisq(N * summary(lm_v)$r.squared, df = q, lower.tail = FALSE)
}

neural2_test <- function(y, Z, e) {
  lags <- Z[, -1, drop = FALSE];  m <- ncol(lags)
  W    <- matrix(runif(20 * m, -2, 2), 20, m)
  H    <- apply(W, 1, function(w) plogis(lags %*% w))
  Phi  <- prcomp(H, center = TRUE, scale. = FALSE)$x[, 2:(m + 2), drop = FALSE]
  lm_n <- lm(e ~ cbind(Z, Phi) - 1)
  pchisq(N * summary(lm_n)$r.squared, df = m + 1, lower.tail = FALSE)
}


# ============================================================
# 5.  Raw p-values
# ============================================================

cat("Running seven nonlinearity tests...\n")

pvals <- c(
  KEENAN  = keenan_test(y_, Z_, e_, yh_),
  TSAY    = tsay_test(y_, Z_),
  RESET   = reset_test(y_, Z_, e_, yh_),
  McLEOD  = mcleod_test(e_),
  BDS     = bds_test2(e_),
  V23     = v23_test(y_, Z_, e_),
  NEURAL2 = neural2_test(y_, Z_, e_)
)

cat("Done.\n\n")


# ============================================================
# 6.  Multiplicity adjustments: Bonferroni, Holm, Bootstrap
# ============================================================

k <- length(pvals)

# Bonferroni
p_bonferroni <- pmin(k * pvals, 1)

# Holm (step-down Bonferroni), recursive formula:
#   q_(1) = min(1, k * p_(1))
#   q_(i) = min(1, max((k - i + 1) * p_(i), q_(i-1)))
ord    <- order(pvals)
p_sort <- pvals[ord]
q_holm <- numeric(k)
q_holm[1] <- min(1, k * p_sort[1])
for (i in 2:k)
  q_holm[i] <- min(1, max((k - i + 1) * p_sort[i], q_holm[i - 1]))
p_holm <- q_holm[order(ord)]

# Bootstrap (residual bootstrap, B replications)
# B = 199 used here.
# Psaradakis uses B = 999 for Table 2; increase here for production use.
B         <- 199
alpha_hat <- as.numeric(b_[1])
phi_hat   <- as.numeric(b_[-1])

cat(sprintf("Running %d bootstrap replications...\n", B))
t0         <- proc.time()
min_p_boot <- numeric(B)

for (r in seq_len(B)) {
  e_star <- sample(e_ - mean(e_), N, replace = TRUE)   # recentred residuals
  xs <- numeric(N + m)
  xs[seq_len(m)] <- x[seq_len(m)]
  for (t in seq(m + 1L, N + m))
    xs[t] <- alpha_hat + sum(phi_hat * xs[seq(t - 1L, t - m)]) + e_star[t - m]
  xb <- xs[seq(m + 1L, N + m)]
  Xb  <- embed(xb, m + 1)
  yb  <- Xb[, 1];   Zb <- cbind(1, Xb[, -1])
  bb  <- solve(crossprod(Zb), crossprod(Zb, yb))
  eb  <- as.numeric(yb - Zb %*% bb)
  yhb <- as.numeric(Zb %*% bb)
  pv_b <- c(
    keenan_test(yb, Zb, eb, yhb),
    tsay_test(yb, Zb),
    reset_test(yb, Zb, eb, yhb),
    mcleod_test(eb),
    bds_test2(eb),
    v23_test(yb, Zb, eb),
    neural2_test(yb, Zb, eb)
  )
  min_p_boot[r] <- min(pv_b)
}

cat(sprintf("Done  (%.1f sec).\n\n", (proc.time() - t0)["elapsed"]))

p_bootstrap <- vapply(pvals, function(pi) mean(min_p_boot <= pi), numeric(1))


# ============================================================
# 7.  Results
# ============================================================

results <- as.data.frame(rbind(
  None       = round(pvals,        3),
  Bonferroni = round(p_bonferroni, 3),
  Holm       = round(p_holm,       3),
  Bootstrap  = round(p_bootstrap,  3)
))

cat(sprintf("=== Replication  (FRED GDPC1, current vintage, m = %d) ===\n\n", m))
print(results)


# ============================================================
# 8.  Takeaways
# ============================================================

# 1. McLEOD-Li is the only test with any signal.
#    All six others are comfortably non-significant before
#    and after adjustment.
#
# 2. Holm dominates Bonferroni uniformly: same assumptions,
#    always at least as powerful. Here they agree because
#    McLEOD is the only borderline case.
#
# 3. The bootstrap accounts for dependence among the seven
#    test statistics and is the most powerful of the three.
#    With revised data, even the bootstrap cannot push
#    McLEOD below 5% after adjustment.
#
# 4. Sources of discrepancy from Psaradakis (2001):
#    -- GDP data revisions: raw McLEOD p moves from 0.005
#       to ~0.010, enough to change the adjusted conclusion.
#    -- Bootstrap: B = 199 is noisy; increase to 999 for
#       production use.





