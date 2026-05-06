rm(list = ls())

# ---------------------------------------------------------------
#  Replication of Hothorn, Bretz & Westfall (2008), Section 6.2
#  "Prediction of total body fat"
# ---------------------------------------------------------------
#  n = 71 healthy German women (Garcia et al., 2005)
#  Response:
#    DEXfat       -- body fat mass (kg) by Dual Energy X-Ray
#                    Absorptiometry (DXA), the expensive "gold standard"
#  Predictors (p = 9):
#    age          -- age in years
#    waistcirc    -- waist circumference (cm)
#    hipcirc      -- hip circumference (cm)
#    elbowbreadth -- elbow breadth (cm), a bone-size measure
#    kneebreadth  -- knee breadth (cm), a bone-size measure
#    anthro3a     -- log(chin SF) + log(triceps SF) + log(subscapular SF)
#    anthro3b     -- log(triceps SF) + log(subscapular SF) + log(abdominal SF)
#    anthro3c     -- log-sum of three other skinfold sites (candidate)
#    anthro4      -- log-sum of four skinfold sites (candidate)
# ---------------------------------------------------------------


library(multcomp)
data("bodyfat")


# ---------------------------------------------------------------
#  (a) Fit model and collect multcomp reference results
# ---------------------------------------------------------------
lmod <- lm(DEXfat ~ ., data = bodyfat)
summary(lmod)                          # unadjusted t-tests & global F

K <- cbind(0, diag(length(coef(lmod)) - 1))
rownames(K) <- names(coef(lmod))[-1]
lmod_glht <- glht(lmod, linfct = K)

summary(lmod_glht, test = Ftest())     # joint F-test
summary(lmod_glht)                     # single-step max-t (multcomp SS)

confint(lmod_glht)                     # simultaneous 95% CIs
par(mar = c(5, 10, 4, 2))
plot(confint(lmod_glht), main = "Simultaneous 95% CIs for slopes")

# Store multcomp p-values now -- used for comparison at the end
p_SS_multcomp <- summary(lmod_glht)$test$pvalues
p_SD_multcomp <- summary(lmod_glht, test = adjusted(type = "free"))$test$pvalues


# ---------------------------------------------------------------
#  (b) Extract ingredients for MC simulation
# ---------------------------------------------------------------
n <- nrow(bodyfat)          # 71 observations
k <- length(coef(lmod))     # 10 parameters (intercept + 9 slopes)
m <- k - 1                  # 9 slope hypotheses

# Observed t-statistics for the 9 slopes (drop the intercept row)
t_obs     <- coef(summary(lmod))[-1, "t value"]
T_max_obs <- max(abs(t_obs))

{
cat("Observed t-statistics:\n")
print(round(t_obs, 3))
cat("Observed T_max =", round(T_max_obs, 3), "\n\n")
}


# ---------------------------------------------------------------
#  Precompute fixed matrices (once, outside the loop)
#
#  Under H_0 with epsilon = sigma * u, u ~ N_n(0, I_n):
#
#      t_j = [ e_j' (X'X)^{-1} X' u / sqrt([(X'X)^{-1}]_{jj}) ]
#            / sqrt( u'M_Xu / (n-k) )
#
#  Stack all m slopes:
#
#      t_vec = A %*% u / sqrt( u'M_Xu / (n-k) )
#
#  where A = D^{-1/2} (X'X)^{-1} X'  (m x n),
#  D = diag([(X'X)^{-1}]_{jj}, j = 2, ..., k).
# ---------------------------------------------------------------
X       <- model.matrix(lmod)
XtX_inv <- solve(crossprod(X))
M_X     <- diag(n) - X %*% XtX_inv %*% t(X)   # residual maker (n x n)

d <- sqrt(diag(XtX_inv)[-1])                    # sqrt of diagonal, slopes only
A <- sweep((XtX_inv %*% t(X))[-1, ], MARGIN = 1, STATS = d, FUN = "/")    # m x n

# Verify: A %*% y / sqrt(RSS/(n-k)) recovers t_obs to machine precision
y       <- bodyfat$DEXfat
t_check <- as.vector(A %*% y) / sqrt(sum(resid(lmod)^2) / (n - k))
stopifnot(max(abs(t_check - t_obs)) < 1e-10)
cat("Verification passed: formula recovers observed t-statistics.\n\n")


# ---------------------------------------------------------------
#  (c) MC simulation
#  Build step-down maxima inside the loop.
#
#  Sort observed |t_j| once, largest to smallest.
#  For each draw b, reorder |t_sim| according to the same observed
#  hypothesis ordering, then build the step-specific maxima backward:
#
#    M_SD[b, m]    = |t_{(m)}^(b)|
#    M_SD[b, step] = max(M_SD[b, step+1], |t_{(step)}^(b)|)
#
#  Column step of M_SD is the reference distribution for step step.
#  Column 1 is T_max -- the single-step reference.
# ---------------------------------------------------------------
B <- 1000
set.seed(2025)

ord      <- order(abs(t_obs), decreasing = TRUE)   # observed hypothesis ordering
t_sorted <- abs(t_obs)[ord]                        # observed |t_(1)| >= ... >= |t_(m)|

M_SD <- matrix(0, nrow = B - 1, ncol = m)

for (b in seq_len(B - 1)) {
  u <- rnorm(n)

  t_sim <- as.vector(A %*% u) /
    sqrt(drop(crossprod(u, M_X %*% u)) / (n - k))

  t_sorted_b <- abs(t_sim)[ord]

  # Backward recursion: step-specific reduced maxima
  M_SD[b, m] <- t_sorted_b[m]
  for (step in (m - 1):1) {
    M_SD[b, step] <- max(M_SD[b, step + 1], t_sorted_b[step])
  }
}

# Column 1 = T_max_sim = single-step reference distribution
T_max_sim <- M_SD[, 1]


# ---------------------------------------------------------------
#  (d) Single-step adjusted p-values
#  Every hypothesis is ranked against column 1, the global T_max.
# ---------------------------------------------------------------
p_SS_sorted <- sapply(seq_len(m), function(step) {
  ref   <- c(T_max_sim, t_sorted[step])
  R_hat <- rank(ref, ties.method = "first")[B]
  (B - R_hat + 1) / B
})

p_SS_MC        <- numeric(m)
p_SS_MC[ord]   <- p_SS_sorted
names(p_SS_MC) <- names(t_obs)


# ---------------------------------------------------------------
#  (e) Step-down: raw p-values then monotonicity enforcement
#  Hypothesis at step k is ranked against column k of M_SD.
# ---------------------------------------------------------------
p_raw_sorted <- sapply(seq_len(m), function(step) {
  ref   <- c(M_SD[, step], t_sorted[step])
  R_hat <- rank(ref, ties.method = "first")[B]
  (B - R_hat + 1) / B
})

# Enforce monotonicity
p_SD_sorted    <- numeric(m)
p_SD_sorted[1] <- p_raw_sorted[1]
for (step in 2:m) {
  p_SD_sorted[step] <- max(p_SD_sorted[step - 1], p_raw_sorted[step])
}

# Recover original order
p_SD_MC        <- numeric(m)
p_SD_MC[ord]   <- p_SD_sorted
names(p_SD_MC) <- names(t_obs)


# ---------------------------------------------------------------
#  Results table
# ---------------------------------------------------------------
result <- data.frame(
  t_obs         = round(t_obs,                                3),
  p_unadj       = round(coef(summary(lmod))[-1, "Pr(>|t|)"], 3),
  p_SS_multcomp = round(p_SS_multcomp,                        3),
  p_SS_MC       = round(p_SS_MC,                              3),
  p_SD_multcomp = round(p_SD_multcomp,                        3),
  p_SD_MC       = round(p_SD_MC,                              3)
)

{
cat("Single-step and step-down max-t adjusted p-values:\n")
print(result)
cat("\nMax abs diff MC SS vs multcomp SS:", round(max(abs(p_SS_MC - p_SS_multcomp)), 4))
cat("\nMax abs diff MC SD vs multcomp SD:", round(max(abs(p_SD_MC - p_SD_multcomp)), 4), "\n")
}


# ---------------------------------------------------------------
#  Takeaways
# ---------------------------------------------------------------

# 1. p_SS_MC and p_SS_multcomp agree to 2-3 decimal places.
#    p_SD_MC and p_SD_multcomp also agree closely.
#    Both pairs target the same joint null distribution --
#    multcomp via the Genz-Bretz integral, MC via simulation.

# 2. For the same simulated reference distribution, step-down adjusted
#    p-values are no larger than the corresponding single-step adjusted
#    p-values. Gains are small here but arise from using reduced-set
#    maxima after each rejection step.

# 3. The only structural change from SS to SD: build M_SD inside
#    the loop using the backward recursion, then rank each
#    |t_{(step)}| against its own column rather than column 1.
#    The same B-1 draws are reused -- no extra simulation needed.

# 4. kneebreadth (raw p = 0.018) does not survive either procedure.
#    waistcirc and hipcirc survive both.






