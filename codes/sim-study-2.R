library(survey)
library(data.table)
library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)
library(doRNG)

source("codes/functions.R")

sims <- 1000
cores <- 8

set.seed(2026)
N <- 20000
n <- 1000
x1 <- rnorm(N, 1, 1)
x2 <- rexp(N, 1)
alp <- rnorm(N)
epsilon <- rnorm(N)
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)

## helper variables for threshold/rank effects
qx1_20_80 <- quantile(x1, probs = c(0.20, 0.80))
qx2_20_80 <- quantile(x2, probs = c(0.20, 0.80))
x1_low20    <- as.numeric(x1 <= qx1_20_80[1])
x1_high80   <- as.numeric(x1 >  qx1_20_80[2])
x2_low20    <- as.numeric(x2 <= qx2_20_80[1])
x2_high80   <- as.numeric(x2 >  qx2_20_80[2])
x2_mid40_60 <- as.numeric(x2 > quantile(x2, 0.40) & x2 <= quantile(x2, 0.60))
Fx2 <- rank(x2, ties.method = "average") / N

## outcomes — threshold and rank-based
y11 <- 1 + 0.5 * x1 + x2 + 2.0 * x2_high80 - 1.0 * x1_low20 + alp + epsilon
y12 <- 1 + 0.5 * x1 + sin(2 * pi * Fx2) + alp + epsilon
y21 <- rbinom(N, 1, plogis(-1 + 0.5 * x1 + 0.5 * x2 + 2.0 * x2_high80 - 1.0 * x1_low20 + alp))
y22 <- rbinom(N, 1, plogis(-0.5 + 0.5 * x1 + 1.2 * sin(2 * pi * Fx2) + alp))

## selection mechanisms
p1 <- plogis(-2 + 2.5 * x2_high80 + 1.5 * x1_low20)
p2 <- plogis(1.5 - 3.0 * x2_mid40_60)

pop_data <- data.frame(x1, x2, y11, y12, y21, y22, p1, p2, w = N / n)

control_ipw <- control_sel(est_method = "gee",
                           nleqslv_method = "Newton",
                           nleqslv_global = "qline")

formulas <- list(
  y = ~ y11 + y12 + y21 + y22,

  lin = y11 + y12 ~ x1 + x2,
  bin = y21 + y22 ~ x1 + x2,
  ipw = ~ x1 + x2,

  ipw_q1 = ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  ipw_q2 = ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q1 = y11 + y12 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q1 = y21 + y22 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q2 = y11 + y12 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q2 = y21 + y22 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,

  ipw_d1 = ~ x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  ipw_d2 = ~ x1 + x2 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  lin_d1 = y11 + y12 ~ x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  bin_d1 = y21 + y22 ~ x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  lin_d2 = y11 + y12 ~ x1 + x2 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  bin_d2 = y21 + y22 ~ x1 + x2 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 +
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 +
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9
)

configs <- list(
  list(basis = "main",           ipw = "ipw",    lin = "lin",    bin = "bin"),
  list(basis = "quartiles",      ipw = "ipw_q1", lin = "lin_q1", bin = "bin_q1"),
  list(basis = "main+quartiles", ipw = "ipw_q2", lin = "lin_q2", bin = "bin_q2"),
  list(basis = "deciles",        ipw = "ipw_d1", lin = "lin_d1", bin = "bin_d1"),
  list(basis = "main+deciles",   ipw = "ipw_d2", lin = "lin_d2", bin = "bin_d2")
)

a <- Sys.time()
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation <- foreach(k = 1:sims,
                              .combine = rbind,
                              .packages = c("survey", "nonprobsvy", "data.table"),
                              .options.snow = opts,
                              .errorhandling = "remove"
) %dorng% {

  sample_prob <- pop_data[sample(1:N, n), ]
  sample_bd1  <- pop_data[rbinom(N, 1, pop_data$p1) == 1, ]
  sample_bd2  <- pop_data[rbinom(N, 1, pop_data$p2) == 1, ]

  sample_prob_svy <- svydesign(ids = ~1, weights = ~w, data = sample_prob)
  q1_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles1)
  q2_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles2)

  aug_data <- function(dat) {
    b1 <- make_quantile_basis(dat, q1_est)
    b2 <- make_quantile_basis(dat, q2_est)
    cbind(dat, b1$cumulative, b2$cumulative)
  }

  sample_prob_aug <- aug_data(sample_prob)
  sample_bd1_aug  <- aug_data(sample_bd1)
  sample_bd2_aug  <- aug_data(sample_bd2)

  sample_prob_aug_svy <- svydesign(ids = ~1, weights = ~w, data = sample_prob_aug)

  datasets <- list(bd1 = sample_bd1_aug, bd2 = sample_bd2_aug)
  results_list <- list()

  for (cfg in configs) {
    for (ds_name in names(datasets)) {
      ds <- datasets[[ds_name]]

      ## IPW
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]], target = formulas$y,
                     data = ds, svydesign = sample_prob_aug_svy,
                     control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "ipw"
        r$dataset   <- ds_name
        r$basis     <- cfg$basis
        r$k         <- k
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## MI linear
      res <- tryCatch({
        m <- nonprob(outcome = formulas[[cfg$lin]],
                     data = ds, svydesign = sample_prob_aug_svy)
        r <- extract(m)
        r$estimator <- "mi"
        r$dataset   <- ds_name
        r$basis     <- cfg$basis
        r$k         <- k
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## MI binomial
      res <- tryCatch({
        m <- nonprob(outcome = formulas[[cfg$bin]], family = "binomial",
                     data = ds, svydesign = sample_prob_aug_svy)
        r <- extract(m)
        r$estimator <- "mi"
        r$dataset   <- ds_name
        r$basis     <- cfg$basis
        r$k         <- k
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## DR linear
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]],
                     outcome = formulas[[cfg$lin]],
                     data = ds, svydesign = sample_prob_aug_svy,
                     control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "dr"
        r$dataset   <- ds_name
        r$basis     <- cfg$basis
        r$k         <- k
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## DR binomial
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]],
                     outcome = formulas[[cfg$bin]], family = "binomial",
                     data = ds, svydesign = sample_prob_aug_svy,
                     control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "dr"
        r$dataset   <- ds_name
        r$basis     <- cfg$basis
        r$k         <- k
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))
    }
  }

  rbindlist(Filter(Negate(is.null), results_list))
}
stopCluster(cl)
print(Sys.time() - a)

saveRDS(results_simulation, "results/sim-results-2.rds")

## ---- multimodal mixture scenario (for future use) ----
# x1_mix <- ifelse(rbinom(N, 1, 0.5) == 1, rnorm(N, -1, 0.4), rnorm(N, 2, 0.6))
# z_mix <- rbinom(N, 1, 0.7)
# x2_mix <- ifelse(z_mix == 1, rexp(N, rate = 1.0), rexp(N, rate = 0.2))
# qx1_mix_20_80 <- quantile(x1_mix, probs = c(0.20, 0.80))
# qx2_mix_20_80 <- quantile(x2_mix, probs = c(0.20, 0.80))
# x1_mix_low20  <- as.numeric(x1_mix <= qx1_mix_20_80[1])
# x2_mix_high80 <- as.numeric(x2_mix >  qx2_mix_20_80[2])
# y15_mix <- 1 + 0.5 * x1_mix + x2_mix + 2.0 * x2_mix_high80 - 1.0 * x1_mix_low20 + rnorm(N) + rnorm(N)
# p7_mix <- plogis(-2 + 2.0 * x2_mix_high80 + 1.5 * x1_mix_low20)
