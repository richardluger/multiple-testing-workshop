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
#                    [composite in Garcia et al.'s final women's equation]
#    anthro3b     -- log(triceps SF) + log(subscapular SF) + log(abdominal SF)
#                    [composite in Garcia et al.'s final men's equation]
#    anthro3c     -- log-sum of three other skinfold sites (candidate)
#    anthro4      -- log-sum of four skinfold sites (candidate)
#  All four anthro variables are log-scale composites of skinfold
#  thickness measurements taken at different body sites.
# ---------------------------------------------------------------

## Package
library(multcomp)

## Data
data("bodyfat")
head(bodyfat)

## Model
lmod <- lm(DEXfat ~ ., data = bodyfat)
summary(lmod)                          # unadjusted t-tests & global F

# ---------------------------------------------------------------
#  Raw p-values for the 9 slopes (drop intercept row)
# ---------------------------------------------------------------
praw <- coef(summary(lmod))[-1, "Pr(>|t|)"]
k    <- length(praw)                   # k = 9


# ---------------------------------------------------------------
#  Bonferroni: multiply each p by k, cap at 1
# ---------------------------------------------------------------
p_bonf <- pmin(k * praw, 1)


# ---------------------------------------------------------------
#  Holm (step-down Bonferroni), recursive formula:
#    q_(1) = min(1, (k) * p_(1))
#    q_(i) = min(1, max((k - i + 1) * p_(i), q_(i-1)))
#  i.e. each adjusted p-value is the larger of the current
#  weighted raw p and the previous adjusted p, capped at 1.
# ---------------------------------------------------------------
ord      <- order(praw)                # indices that sort ascending
p_sorted <- praw[ord]
q_holm   <- numeric(k)
q_holm[1] <- min(1, k * p_sorted[1])
for (i in 2:k) {
  q_holm[i] <- min(1, max((k - i + 1) * p_sorted[i], q_holm[i - 1]))
}
p_holm <- q_holm[order(ord)]          # back to original order


# ---------------------------------------------------------------
#  Side-by-side comparison table
# ---------------------------------------------------------------
result <- data.frame(raw  = round(praw,   3), Bonf = round(p_bonf, 3), Holm = round(p_holm, 3))
print(result)


# ---------------------------------------------------------------
#  Takeaways
# ---------------------------------------------------------------

# 1. Clear signals survive everything.
#    hipcirc and waistcirc are significant under every procedure.
#    Adjustment does not suppress genuinely strong effects.

# 2. Clear noise goes nowhere.
#    age, elbowbreadth, anthro3c, anthro4 are never significant.
#    No method rescues a weak signal.

# 3. The action is at the margin -- and that is exactly where it matters.
#    kneebreadth looks significant raw (p = 0.018) but loses its star
#    under every FWER procedure (Bonf: 0.165, Holm: 0.128).
#    This is the variable you would have believed in had you only read
#    the unadjusted table.

# 4. Bonferroni and Holm: valid under arbitrary dependence.
#    They use only the union bound -- no assumption on the joint
#    distribution of the test statistics whatsoever.
#    Holm dominates Bonferroni uniformly: same assumptions, same strong
#    FWER control, always at least as powerful. There is no reason to
#    use Bonferroni over Holm.






