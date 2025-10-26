#################################
# run_single_replicate_debug.R  #
# Author: Caleb Ives            #
#################################

rm(list = ls())
library(bayesGT)
library(groupTesting)
library(profvis)


##########################################
# Model and Link Function Configuration
##########################################

link_tag  <- "logit"
model_tag <- "M1"

beta_val <- switch(model_tag,
                   "M1" = c(-3, 2),
                   "M2" = c(-3, 2, -1),
                   "M3" = c(-3, 2, -0.5))

glmLink <- function(fn.name = c("logit", "probit", "cloglog")) {
  fn.name <- match.arg(fn.name)
  clip <- function(x, eps = 1e-10) pmax(pmin(x, 1 - eps), eps)
  if (fn.name == "logit") {
    g <- function(u) stats::plogis(u)
  } else if (fn.name == "probit") {
    g <- function(u) stats::pnorm(u)
  } else {
    g <- function(u) 1 - exp(-exp(u))
  }
  return(list(name = fn.name, g = g))
}

link_obj <- glmLink(link_tag)
link_fun <- link_obj$g


##########################################
# Create Settings Object
##########################################

settings <- create_settings(list(
  ## Data structure
  N         = 5000,
  psz       = c(5, 1),
  assay_id  = c(1, 2),

  ## True parameters
  beta_true = beta_val,
  se_true   = c(0.95, 0.98),
  sp_true   = c(0.98, 0.99),

  ## Se/Sp estimation
  known_acc = FALSE,
  se_init   = 0.9,
  sp_init   = 0.9,

  ## Sampling configuration
  gibbs_iter = 6000,
  burn_in    = 1000,
  mh_control = list(
    type = "rw",
    args = list(proposal_sd = 0.05)
  ),
  link_fcn   = link_fun,
  link_name  = link_tag,

  ## Session parameters
  nsims      = 1,
  seed       = 123,
  alpha      = 0.05,
  keep_raw   = TRUE
))


##########################################
# Simulate Data for One Replicate
##########################################

set.seed(settings$seed)
x1 <- rnorm(settings$N)
x2 <- rbinom(settings$N, 1, 0.5)

X <- switch(model_tag,
            "M1" = cbind(1, x1),
            "M2" = cbind(1, x1, x2),
            "M3" = cbind(1, x1, x1^2))

p.t <- settings$link_fcn(X %*% settings$beta_true)

sim_out <- hier.gt.simulation(
  N       = settings$N,
  p       = p.t,
  S       = length(settings$psz),
  psz     = settings$psz,
  Se      = settings$se_true,
  Sp      = settings$sp_true,
  assayID = settings$assay_id
)

test_data <- list(X = X, Z = sim_out$gtData, tsts = sim_out$testsExp)


##########################################
# Run inference with profvis
##########################################

profvis({
  results <- infer_posterior(settings, test_data)
  print_results(results, digits = 4)
})


## END OF FILE ##
