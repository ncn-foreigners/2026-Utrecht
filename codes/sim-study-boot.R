library(survey)
library(data.table)
library(nonprobsvy)
library(foreach)
library(doRNG)
library(doFuture)
library(progressr)

source("codes/functions.R")

## NOTE: plan(multicore) uses forking — run from terminal, NOT RStudio
## For testing set sims/B small; for production use sims = 500, B = 50
sims  <- 8
cores <- 8
B     <- 10

set.seed(2026)
N <- 20000
n <- 1000
x1 <- rnorm(N, 1, 1)
x2 <- rexp(N, 1)
alp <- rnorm(N)
epsilon <- rnorm(N)
y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5 * (x1 - 1.5)^2 + x2^2 + alp + epsilon
y21 <- rbinom(N, 1, plogis(1 + x1 + x2 + alp))
y22 <- rbinom(N, 1, plogis(0.5 * (x1 - 1.5)^2 + x2^2 + alp))
p1 <- plogis(x2)
p2 <- plogis(-3 + (x1 - 1.5)^2 + (x2 - 2)^2)
pop_data <- data.frame(x1, x2, y11, y12, y21, y22, p1, p2, w = N / n)
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)

control_ipw <- control_sel(est_method = "gee",
                           nleqslv_method = "Newton",
                           nleqslv_global = "qline")

formulas <- list(
  y = ~ y11 + y12 + y21 + y22,

  ## models with main effects
  lin = y11 + y12 ~ x1 + x2,
  bin = y21 + y22 ~ x1 + x2,
  ipw = ~ x1 + x2,

  ## models with quartiles
  ipw_q1 = ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  ipw_q2 = ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q1 = y11 + y12 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q1 = y21 + y22 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q2 = y11 + y12 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q2 = y21 + y22 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,

  ## models with deciles
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

## config for simulation loop
configs <- list(
  list(basis = "main",           ipw = "ipw",    lin = "lin",    bin = "bin"),
  list(basis = "quartiles",      ipw = "ipw_q1", lin = "lin_q1", bin = "bin_q1"),
  list(basis = "main+quartiles", ipw = "ipw_q2", lin = "lin_q2", bin = "bin_q2"),
  list(basis = "deciles",        ipw = "ipw_d1", lin = "lin_d1", bin = "bin_d1"),
  list(basis = "main+deciles",   ipw = "ipw_d2", lin = "lin_d2", bin = "bin_d2")
)

## helper: augment a data.frame with quantile/decile indicator bases
aug_data <- function(dat, q1, q2) {
  b1 <- make_quantile_basis(dat, q1)
  b2 <- make_quantile_basis(dat, q2)
  cbind(dat, b1$cumulative, b2$cumulative)
}

## helper: fit all models for given augmented data
fit_all <- function(prob_svy, bd1, bd2) {
  results_list <- list()
  datasets <- list(bd1 = bd1, bd2 = bd2)

  for (cfg in configs) {
    for (ds_name in names(datasets)) {
      ds <- datasets[[ds_name]]

      ## IPW
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]], target = formulas$y,
                     data = ds, svydesign = prob_svy,
                     pop_size = N, control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "ipw"; r$dataset <- ds_name; r$basis <- cfg$basis
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## MI linear
      res <- tryCatch({
        m <- nonprob(outcome = formulas[[cfg$lin]],
                     pop_size = N, data = ds, svydesign = prob_svy)
        r <- extract(m)
        r$estimator <- "mi"; r$dataset <- ds_name; r$basis <- cfg$basis
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## MI binomial
      res <- tryCatch({
        m <- nonprob(outcome = formulas[[cfg$bin]], family = "binomial",
                     pop_size = N, data = ds, svydesign = prob_svy)
        r <- extract(m)
        r$estimator <- "mi"; r$dataset <- ds_name; r$basis <- cfg$basis
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## DR linear
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]],
                     outcome = formulas[[cfg$lin]],
                     data = ds, svydesign = prob_svy,
                     pop_size = N, control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "dr"; r$dataset <- ds_name; r$basis <- cfg$basis
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))

      ## DR binomial
      res <- tryCatch({
        m <- nonprob(selection = formulas[[cfg$ipw]],
                     outcome = formulas[[cfg$bin]], family = "binomial",
                     data = ds, svydesign = prob_svy,
                     pop_size = N, control_selection = control_ipw)
        r <- extract(m)
        r$estimator <- "dr"; r$dataset <- ds_name; r$basis <- cfg$basis
        r
      }, error = function(e) NULL)
      results_list <- c(results_list, list(res))
    }
  }

  rbindlist(Filter(Negate(is.null), results_list))
}

## ============================================================================
## PASS 1: draw samples and estimate quantiles (sequential — fast, no models)
## ============================================================================
cat("Pass 1: drawing samples...\n")
a <- Sys.time()

set.seed(2026 + 1)
samples_list <- vector("list", sims)
for (k in 1:sims) {
  sample_prob <- pop_data[sample(1:N, n), ]
  sample_bd1  <- pop_data[rbinom(N, 1, pop_data$p1) == 1, ]
  sample_bd2  <- pop_data[rbinom(N, 1, pop_data$p2) == 1, ]

  sample_prob_svy <- svydesign(ids = ~1, weights = ~w, data = sample_prob)
  q1_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles1)
  q2_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles2)

  samples_list[[k]] <- list(
    sample_prob = sample_prob,
    sample_bd1  = sample_bd1,
    sample_bd2  = sample_bd2,
    q1_est      = q1_est,
    q2_est      = q2_est
  )
}
cat(sprintf("Pass 1 done in %s\n", format(Sys.time() - a)))

## ============================================================================
## PASS 2: fit models — parallel over all (k, b) pairs using future (forking)
## boot = 0  : point estimate from original samples
## boot > 0  : bootstrap replicate (resample all three datasets)
## ============================================================================
kb_grid <- expand.grid(k = seq_len(sims), b = 0:B)
n_tasks <- nrow(kb_grid)

cat(sprintf("Pass 2: fitting models (%d tasks = %d reps x %d boot)...\n",
            n_tasks, sims, B + 1))

plan(multicore, workers = cores)
registerDoFuture()
handlers(global = TRUE)
handlers("txtprogressbar")

b2 <- Sys.time()

results_simulation <- with_progress({
  p <- progressor(steps = n_tasks)

  foreach(
    i = 1:n_tasks,
    .combine = rbind,
    .errorhandling = "remove"
  ) %dorng% {
    k_idx <- kb_grid$k[i]
    b_idx <- kb_grid$b[i]
    smp   <- samples_list[[k_idx]]

    if (b_idx == 0L) {
      ## ---- point estimate: use original samples directly ----
      prob_aug     <- aug_data(smp$sample_prob, smp$q1_est, smp$q2_est)
      bd1_aug      <- aug_data(smp$sample_bd1,  smp$q1_est, smp$q2_est)
      bd2_aug      <- aug_data(smp$sample_bd2,  smp$q1_est, smp$q2_est)
      prob_aug_svy <- svydesign(ids = ~1, weights = ~w, data = prob_aug)

      res <- fit_all(prob_aug_svy, bd1_aug, bd2_aug)
      res[, `:=`(k = k_idx, boot = 0L)]

    } else {
      ## ---- bootstrap replicate: resample all three datasets ----
      sp  <- smp$sample_prob
      sb1 <- smp$sample_bd1
      sb2 <- smp$sample_bd2

      prob_boot <- sp[sample.int(nrow(sp),   replace = TRUE), ]
      bd1_boot  <- sb1[sample.int(nrow(sb1), replace = TRUE), ]
      bd2_boot  <- sb2[sample.int(nrow(sb2), replace = TRUE), ]

      ## re-estimate quantiles from bootstrap probability sample
      prob_boot_svy <- svydesign(ids = ~1, weights = ~w, data = prob_boot)
      q1_boot <- tryCatch(svyquantile(formulas$ipw, prob_boot_svy, p_quantiles1),
                          error = function(e) smp$q1_est)
      q2_boot <- tryCatch(svyquantile(formulas$ipw, prob_boot_svy, p_quantiles2),
                          error = function(e) smp$q2_est)

      ## augment bootstrap samples with bootstrap quantile bases
      prob_boot_aug     <- aug_data(prob_boot, q1_boot, q2_boot)
      bd1_boot_aug      <- aug_data(bd1_boot,  q1_boot, q2_boot)
      bd2_boot_aug      <- aug_data(bd2_boot,  q1_boot, q2_boot)
      prob_boot_aug_svy <- svydesign(ids = ~1, weights = ~w, data = prob_boot_aug)

      res <- tryCatch(
        fit_all(prob_boot_aug_svy, bd1_boot_aug, bd2_boot_aug),
        error = function(e) NULL
      )
      if (is.null(res)) { p(); return(NULL) }
      res[, `:=`(k = k_idx, boot = b_idx)]
    }

    p()
    res
  }
})

plan(sequential)
cat(sprintf("Pass 2 done in %s\n", format(Sys.time() - b2)))
cat(sprintf("Total time: %s\n", format(Sys.time() - a)))

saveRDS(results_simulation, "results/sim-results-boot.rds")
