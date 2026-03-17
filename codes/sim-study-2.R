set.seed(2026)

N <- 20000
n <- 1000

x1 <- rnorm(N, 1, 1)
x2 <- rexp(N, 1)
alp <- rnorm(N)
epsilon <- rnorm(N)
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)


# ============================================================
# NEW ADDITIONS
# designed to create situations where quantile balancing
# / quantile calibration / quantile-based nuisance models
# should have an advantage over purely smooth global models
# ============================================================

# ------------------------------------------------------------
# 1. POPULATION QUANTILES OF X1 AND X2
# Added because quantile-based methods should use cutpoints
# anchored to the target population (or probability sample),
# not the nonprobability sample.
# ------------------------------------------------------------
qx1_10_90 <- quantile(x1, probs = c(0.10, 0.90))
qx1_20_80 <- quantile(x1, probs = c(0.20, 0.80))
qx2_10_90 <- quantile(x2, probs = c(0.10, 0.90))
qx2_20_80 <- quantile(x2, probs = c(0.20, 0.80))
qx2_deciles <- quantile(x2, probs = seq(0.10, 0.90, 0.10))

# Optional helper variables for threshold/rank-type effects
x1_low20  <- as.numeric(x1 <= qx1_20_80[1])
x1_high80 <- as.numeric(x1 >  qx1_20_80[2])
x2_low20  <- as.numeric(x2 <= qx2_20_80[1])
x2_high80 <- as.numeric(x2 >  qx2_20_80[2])

# Middle region indicator for x2
x2_mid40_60 <- as.numeric(x2 > quantile(x2, 0.40) & x2 <= quantile(x2, 0.60))


# ------------------------------------------------------------
# 2. NEW OUTCOME MODELS WITH THRESHOLD / PIECEWISE STRUCTURE
# Added because quantile-based methods should perform well
# when E(Y|X) depends on tail membership or threshold effects
# rather than only on smooth global polynomials.
# ------------------------------------------------------------

# Continuous outcome with threshold effect in upper tail of x2
# and lower tail of x1
y13 <- 1 + 0.5 * x1 + x2 + 2.0 * x2_high80 - 1.0 * x1_low20 + alp + epsilon

# Continuous outcome depending on rank-like oscillation in x2.
# First map x2 to its empirical population CDF rank.
Fx2 <- rank(x2, ties.method = "average") / N
y14 <- 1 + 0.5 * x1 + sin(2 * pi * Fx2) + alp + epsilon

# Binary outcome with threshold structure
y23 <- rbinom(
  N, 1,
  plogis(-1 + 0.5 * x1 + 0.5 * x2 + 2.0 * x2_high80 - 1.0 * x1_low20 + alp)
)

# Binary outcome with nonlinear rank-based pattern
y24 <- rbinom(
  N, 1,
  plogis(-0.5 + 0.5 * x1 + 1.2 * sin(2 * pi * Fx2) + alp)
)


# ------------------------------------------------------------
# 3. NEW SELECTION MECHANISMS FAVORING QUANTILE METHODS
# Added because quantile balancing should help particularly
# when selection depends on tail membership / threshold regions.
#
# IMPORTANT:
# These are designed to preserve overlap, not destroy it.
# Probabilities remain strictly between 0 and 1.
# ------------------------------------------------------------

# p3: tail-based / threshold-based nonprobability mechanism
# overrepresents high x2 and underrepresents low x1
p3 <- plogis(-2 + 2.5 * x2_high80 + 1.5 * x1_low20)

# p4: "missing middle" mechanism
# underrepresents the middle of x2, which is hard for moment-based
# adjustment to detect if mean/variance are similar
p4 <- plogis(1.5 - 3.0 * x2_mid40_60)

# p5: rank-sensitive smooth-but-local mechanism based on decile regions
# stronger selection in upper x2 deciles, but not extreme enough
# to cause zero-overlap problems
x2_gt_70 <- as.numeric(x2 > quantile(x2, 0.70))
x2_gt_90 <- as.numeric(x2 > quantile(x2, 0.90))
p5 <- plogis(-1.5 + 1.2 * x2_gt_70 + 1.0 * x2_gt_90 + 0.5 * x1)

# p6: weak-overlap stress test with bounded probabilities
# This is useful for studying instability without complete positivity failure.
eta6 <- -2.5 + 2.0 * x2_high80 - 1.5 * x2_low20 + 0.3 * x1
eps_overlap <- 0.02
p6 <- eps_overlap + (1 - 2 * eps_overlap) * plogis(eta6)


# ------------------------------------------------------------
# 4. OPTIONAL MULTIMODAL VERSION OF X FOR EXTRA SCENARIOS
# Added because quantile methods can help more when moments
# do not describe the covariate distribution well.
#
# This does NOT replace the baseline x1, x2.
# It gives an alternative DGP for robustness experiments.
# ------------------------------------------------------------

# Example alternative x1 with two modes
x1_mix <- ifelse(
  rbinom(N, 1, 0.5) == 1,
  rnorm(N, -1, 0.4),
  rnorm(N,  2, 0.6)
)

# Example alternative x2 with mixture-exponential behavior
z_mix <- rbinom(N, 1, 0.7)
x2_mix <- ifelse(z_mix == 1, rexp(N, rate = 1.0), rexp(N, rate = 0.2))

# Quantiles for multimodal / mixture scenario
qx1_mix_20_80 <- quantile(x1_mix, probs = c(0.20, 0.80))
qx2_mix_20_80 <- quantile(x2_mix, probs = c(0.20, 0.80))

# Example threshold indicators in mixture scenario
x1_mix_low20  <- as.numeric(x1_mix <= qx1_mix_20_80[1])
x2_mix_high80 <- as.numeric(x2_mix >  qx2_mix_20_80[2])

# Example additional outcome and selection under mixture covariates
y15_mix <- 1 + 0.5 * x1_mix + x2_mix + 2.0 * x2_mix_high80 - 1.0 * x1_mix_low20 +
  rnorm(N) + rnorm(N)

p7_mix <- plogis(-2 + 2.0 * x2_mix_high80 + 1.5 * x1_mix_low20)


# ------------------------------------------------------------
# 5. BIN COUNT DIAGNOSTICS FOR OVERLAP
# Added because quantile methods require that the
# nonprobability sample has at least some support in each
# population quantile region.
#
# These diagnostics should be computed after drawing the
# nonprobability sample B under each p_k.
# ------------------------------------------------------------

# Example helper function:
check_bin_counts <- function(x, probs, s_B) {
  # x   : population covariate
  # probs : quantile probabilities defining bins
  # s_B : logical or 0/1 indicator of membership in B-sample
  qs <- quantile(x, probs = probs)
  brks <- c(-Inf, qs, Inf)
  tab <- table(cut(x[s_B == 1], breaks = brks, include.lowest = TRUE))
  return(tab)
}

# Example of usage after drawing a B-sample:
# s_B3 <- rbinom(N, 1, p3)
# check_bin_counts(x2, probs = seq(0.1, 0.9, 0.1), s_B = s_B3)

# If some bins are empty or nearly empty, use fewer quantiles
# or collapse adjacent bins.


# ------------------------------------------------------------
# 6. FINAL DATA FRAME
# Added all new variables while retaining original ones
# ------------------------------------------------------------

pop_data <- data.frame(
  x1, x2,
  y11, y12, y13, y14,
  y21, y22, y23, y24,
  p1, p2, p3, p4, p5, p6,
  Fx2,
  x1_low20, x1_high80, x2_low20, x2_high80, x2_mid40_60,
  w = N / n
)

# Optional alternative population for multimodal sensitivity analysis
pop_data_mix <- data.frame(
  x1 = x1_mix,
  x2 = x2_mix,
  y15 = y15_mix,
  p7 = p7_mix,
  w = N / n
)