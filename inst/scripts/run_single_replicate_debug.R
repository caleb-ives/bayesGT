#######################
# single_replicate.R  #
# Author: Caleb Ives  #
#######################


rm(list = ls())
library(bayesGT)
library(groupTesting)

## Settings
settings <- list(
  seed       = 123,
  N          = 5000,
  model      = "M1",
  beta_true  = NULL,
  psz        = c(5, 1),
  assay_id   = c(1, 2),
  se_t       = c(0.95, 0.98),
  sp_t       = c(0.98, 0.99),
  known_acc  = FALSE,
  se_0       = 0.9,
  sp_0       = 0.9,
  nsim       = 1,
  post_git   = 6000,
  burn       = 1000,
  alpha      = 0.05,
  g          = stats::plogis,
  keep_raw   = FALSE
)

settings$beta_true <- switch(settings$model,
                             "M1" = c(-3, 2),
                             "M2" = c(-3, 2, -1),
                             "M3" = c(-3, 2, -0.5))

## Simulate a single data set
set.seed(settings$seed)

x1 <- rnorm(settings$N)
x2 <- rbinom(settings$N, 1, 0.5)

X <- switch(settings$model,
            "M1" = cbind(1, x1),
            "M2" = cbind(1, x1, x2),
            "M3" = cbind(1, x1, x1^2))

p.t <- settings$g(X %*% settings$beta_true)

sim_out <- hier.gt.simulation(
  N       = nrow(X),
  p       = p.t,
  S       = length(settings$psz),
  psz     = settings$psz,
  Se      = settings$se_t,
  Sp      = settings$sp_t,
  assayID = settings$assay_id
)

test_dat <- list(
  X    = X,
  Z    = sim_out$gtData,
  tsts = sim_out$testsExp
)

## Inference (no parallel processing)
posterior      <- infer_posterior(test_dat, settings)

## pack what compute_results() expects for a *single* run
stub <- list(
  beta_means = matrix(posterior$beta_mean, nrow = 1),
  beta_sds   = matrix(posterior$beta_sd,   nrow = 1),
  cred_ints  = list(posterior$cred_int),
  tests      = posterior$tests
)

if (!settings$known_acc) {
  stub$se_samps <- list(posterior$se)
  stub$sp_samps <- list(posterior$sp)
}
if (settings$keep_raw) {
  stub$beta_raw <- list(posterior$beta_all)
  if (!settings$known_acc) {
    stub$se_raw <- list(posterior$se)
    stub$sp_raw <- list(posterior$sp)
  }
}

## Output
summary_out <- compute_results(stub, settings)
print_results(summary_out, settings, digits = 2)


## END OF FILE ##
