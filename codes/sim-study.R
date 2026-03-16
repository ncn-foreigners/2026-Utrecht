library(sampling)
library(survey)
library(data.table)
library(nonprobsvy)
library(jointCalib)
library(doSNOW)
library(progress)
library(foreach)
library(laeken)
library(ggplot2)

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
pop_data <- data.frame(x1,x2,y11,y12,y21,y22,p1,p2) |> setDT()
p_quantiles1 <- seq(0.25, 0.75, 0.25)
p_quantiles2 <- seq(0.10, 0.90, 0.10)


## function

set.seed(2026)
sample_prob <- pop_data[sample(1:N, n),]
sample_prob$w <- N/n
sample_prob_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob)
q1_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles1)
q2_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles2)
x_totals <- svytotal( ~ x1 + x2, sample_prob_svy)
sample_bd1 <- pop_data[rbinom(N,1,pop_data$p1)==1, ]
sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]


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


y_list <- list(
  mod1 = ~ y11 + y12,
  mod2 = ~ y21 + y22
)

f_p1a <- make_formulas(basis1, y_list)
f_p1b <- make_formulas(sample_bd1_basis1, y_list, x_vars = c("x1", "x2"))
f_p2a <- make_formulas(sample_bd1_basis2, y_list)
f_p2b <- make_formulas(sample_bd1_basis2, y_list, x_vars = c("x1", "x2"))

sample_bd1_aug <- cbind(sample_bd1, 
                        sample_bd1_basis1$cumulative, 
                        sample_bd1_basis2$cumulative)

sample_bd2_aug <- cbind(sample_bd2, 
                        sample_bd2_basis1$cumulative, 
                        sample_bd2_basis2$cumulative)


test <- nonprob(outcome = f_p1a$mod1,
                data =   sample_bd1_aug,
                svydesign = sample_prob_aug_svy)

test2 <- nonprob(outcome = f_p1a$mod2,
                 data =   sample_bd1_aug,
                 svydesign = sample_prob_aug_svy,
                 family = "binomial")


##

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
