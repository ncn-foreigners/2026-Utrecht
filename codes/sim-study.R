library(survey)
library(data.table)
library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)

## TODO:
## estimators to compare: mi and dr with x1 + x2 + x1^2 + x2^2
## new estimators with quantiles dec / 
## what about vars selection?
## proper example 

source("codes/functions.R")

set.seed(2026)
N <- 20000
n <- 1000
x1 <- rnorm(N,1,1)
x2 <- rexp(N,1)
alp <- rnorm(N)
epsilon <- rnorm(N)
y11 <- 1 + x1 + x2 + alp + epsilon
y12 <- 0.5*(x1-1.5)^2 + x2^2 + alp + epsilon
y21 <- rbinom(N,1,plogis(1 + x1 + x2 + alp))
y22 <- rbinom(N,1,plogis(0.5*(x1-1.5)^2 + x2^2 + alp))
p1 <- plogis(x2)
p2 <- plogis(-3+(x1-1.5)^2+(x2-2)^2)
pop_data <- data.table(x1,x2,y11,y12,y21,y22,p1,p2, w = N/n)
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)

control_ipw <- control_sel(est_method = "gee", 
                           nleqslv_method = "Newton",
                           nleqslv_global = "pwldog")

formulas <- list(
  ## ys
  y = ~ y11 + y12 + y21 + y22,
  
  ## models with main effects
  lin = y11 + y12 ~ x1 + x2,
  bin = y21 + y22 ~ x1 + x2,
  ipw = ~ x1 + x2,
  
  
  ## models with quartiles
  ## ipw with q
  ipw_q1 = ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  ipw_q2 = ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q1 = y11 + y12 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q1 = y21 + y22 ~ x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  lin_q2 = y11 + y12 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  bin_q2 = y21 + y22 ~ x1 + x2 + x1_0.25 + x1_0.5 + x1_0.75 + x2_0.25 + x2_0.5 + x2_0.75,
  
  ## models with deciles
  ipw_q2 = ~ x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + 
    x1_0.7 + x1_0.8 + x1_0.9 + x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + 
    x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9,
  ipw_q2 = ~ x1 + x2 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + 
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

## function

set.seed(2026)
sample_prob <- pop_data[sample(1:N, n),]
sample_bd1 <- pop_data[rbinom(N,1,pop_data$p1)==1, ]
sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]

sample_prob_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob)
q1_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles1)
q2_est <- svyquantile(formulas$ipw, sample_prob_svy, p_quantiles2)
x_totals <- svytotal(formulas$ipw, sample_prob_svy)

sample_prob_basis1 <- make_quantile_basis(sample_prob, q1_est)
sample_prob_basis2 <- make_quantile_basis(sample_prob, q2_est)

sample_prob_aug <- cbind(sample_prob, 
                         sample_prob_basis1$cumulative, 
                         sample_prob_basis2$cumulative)

sample_prob_aug_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob_aug)

sample_bd1_basis1 <- make_quantile_basis(sample_bd1, q1_est)
sample_bd1_basis2 <- make_quantile_basis(sample_bd1, q2_est)
sample_bd2_basis1 <- make_quantile_basis(sample_bd2, q1_est)
sample_bd2_basis2 <- make_quantile_basis(sample_bd2, q2_est)

sample_bd1_aug <- cbind(sample_bd1, 
                        sample_bd1_basis1$cumulative, 
                        sample_bd1_basis2$cumulative)

sample_bd2_aug <- cbind(sample_bd2, 
                        sample_bd2_basis1$cumulative, 
                        sample_bd2_basis2$cumulative)

## estimators standard: mi, ipw (gee), dr (gee)
### ipw
model_bd1_ipw_all_main <- nonprob(selection = formulas$ipw, target = formulas$y, 
                             data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                             control_selection = control_ipw)

model_bd2_ipw_all_main <- nonprob(selection = formulas$ipw, target = formulas$y, 
                             data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                             control_selection = control_ipw)
### mass imputation
model_bd1_mi_lin_main <- nonprob(outcome = formulas$lin,
                            data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_lin_main <- nonprob(outcome = formulas$lin,
                            data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

model_bd1_mi_bin_main <- nonprob(outcome = formulas$bin, family = "binomial",
                            data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_bin_main <- nonprob(outcome = formulas$bin, family = "binomial",
                            data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

### doubly robust
model_bd1_dr_lin_main <- nonprob(selection = formulas$ipw,  
                            outcome = formulas$lin,
                            data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                            control_selection = control_ipw)

model_bd2_dr_lin_main <- nonprob(selection = formulas$ipw,  
                            outcome = formulas$lin,
                            data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                            control_selection = control_ipw)

model_bd1_dr_bin_main <- nonprob(selection = formulas$ipw,  
                            outcome = formulas$bin, family = "binomial",
                            data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                            control_selection = control_ipw)

model_bd2_dr_bin_main <- nonprob(selection = formulas$ipw,  
                            outcome = formulas$bin, family = "binomial",
                            data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                            control_selection = control_ipw)


## estimators with quartiles: mi, ipw (gee), dr (gee)
### ipw
model_bd1_ipw_all_q1 <- nonprob(selection = formulas$ipw_q1, target = formulas$y, 
                                  data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                                  control_selection = control_ipw)

model_bd2_ipw_all_q1 <- nonprob(selection = formulas$ipw_q1, target = formulas$y, 
                                  data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                                  control_selection = control_ipw)
### mass imputation
model_bd1_mi_lin_q1 <- nonprob(outcome = formulas$lin_q1,
                                 data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_lin_q1 <- nonprob(outcome = formulas$lin_q1,
                                 data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

model_bd1_mi_bin_q1 <- nonprob(outcome = formulas$bin_q1, family = "binomial",
                                 data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_bin_q1 <- nonprob(outcome = formulas$bin_q1, family = "binomial",
                                 data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

### doubly robust
model_bd1_dr_lin_q1 <- nonprob(selection = formulas$ipw_q1,  
                                 outcome = formulas$lin_q1,
                                 data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                                 control_selection = control_ipw)

model_bd2_dr_lin_q1 <- nonprob(selection = formulas$ipw_q1,  
                                 outcome = formulas$lin_q1,
                                 data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                                 control_selection = control_ipw)

model_bd1_dr_bin_q1 <- nonprob(selection = formulas$ipw_q1,  
                                 outcome = formulas$bin_q1, family = "binomial",
                                 data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                                 control_selection = control_ipw)

model_bd2_dr_bin_q1 <- nonprob(selection = formulas$ipw_q1,  
                                 outcome = formulas$bin_q1, family = "binomial",
                                 data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                                 control_selection = control_ipw)


## estimators with deciles: mi, ipw (gee), dr (gee)
### ipw
model_bd1_ipw_all_q2 <- nonprob(selection = formulas$ipw_q2, target = formulas$y, 
                                data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                                control_selection = control_ipw)

model_bd2_ipw_all_q2 <- nonprob(selection = formulas$ipw_q2, target = formulas$y, 
                                data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                                control_selection = control_ipw)
### mass imputation
model_bd1_mi_lin_q2 <- nonprob(outcome = formulas$lin_q2,
                               data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_lin_q2 <- nonprob(outcome = formulas$lin_q2,
                               data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

model_bd1_mi_bin_q2 <- nonprob(outcome = formulas$bin_q2, family = "binomial",
                               data = sample_bd1_aug, svydesign = sample_prob_aug_svy)

model_bd2_mi_bin_q2 <- nonprob(outcome = formulas$bin_q2, family = "binomial",
                               data = sample_bd2_aug, svydesign = sample_prob_aug_svy)

### doubly robust
model_bd1_dr_lin_q2 <- nonprob(selection = formulas$ipw_q2,  
                               outcome = formulas$lin_q2,
                               data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                               control_selection = control_ipw)

model_bd2_dr_lin_q2 <- nonprob(selection = formulas$ipw_q2,  
                               outcome = formulas$lin_q2,
                               data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                               control_selection = control_ipw)

model_bd1_dr_bin_q2 <- nonprob(selection = formulas$ipw_q2,  
                               outcome = formulas$bin_q2, family = "binomial",
                               data = sample_bd1_aug, svydesign = sample_prob_aug_svy, 
                               control_selection = control_ipw)

model_bd2_dr_bin_q2 <- nonprob(selection = formulas$ipw_q2,  
                               outcome = formulas$bin_q2, family = "binomial",
                               data = sample_bd2_aug, svydesign = sample_prob_aug_svy, 
                               control_selection = control_ipw)



a <- Sys.time()
sims <- 550
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation1 <- foreach(k=1:sims, 
                               .combine = rbind,
                               .packages = c("survey", "nonprobsvy", "jointCalib", "sampling", "data.table", "laeken"),
                               .options.snow = opts,
                               .errorhandling = "remove"
) %dopar% {
  yang_sim(k)                              
}
stopCluster(cl)
print(Sys.time() - a) 
