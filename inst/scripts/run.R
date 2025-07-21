#######################
# run.R               #
# Author: Caleb Ives  #
#######################


rm(list=ls())
library(bayesGT)
library(groupTesting)


#----------------------------------------------
# Model configuration for testing
#----------------------------------------------

# Using predefined model tags ("M1", "M2", "M3")
# Based on McMahan et al. (2017, Biometrics)
# Used for comparison with published tables
model_tag <- "M2"
beta_val <- switch(model_tag,
                   "M1" = c(-3, 2),
                   "M2" = c(-3, 2, -1),
                   "M3" = c(-3, 2, -0.5)
)

glmLink <- function(fn.name=c("logit","probit","cloglog")){
  fn.name <- match.arg(fn.name)
  clip <- function(x,eps=1e-10) pmax(pmin(x,1-eps),eps)
  if (fn.name == "logit"){
    g <- function(u) stats::plogis(u)
    dg <- function(u) g(u)*(1-g(u))
    d2g <- function(u) dg(u)*(1-2*g(u))
    g.inv <- function(p) stats::qlogis(clip(p))
  }
  if (fn.name == "probit"){
    g <- function(u) stats::pnorm(u)
    dg <- function(u) stats::dnorm(u)
    d2g <- function(u) -u*stats::dnorm(u)
    g.inv <- function(p) stats::qnorm(clip(p))
  }
  if (fn.name == "cloglog"){
    g <- function(u) 1 - exp(-exp(u))
    dg <- function(u) exp(u)*exp(-exp(u))
    d2g <- function(u) exp(u-exp(u))*(1-exp(u))
    g.inv <- function(p) log(-log(1-clip(p)))
  }
  structure(list(g=g,dg=dg,d2g=d2g,gInv=g.inv),class="glmLink")
}


#----------------------------------------------
# bayesGT - Create settings object
#----------------------------------------------

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
  link_fcn = glmLink("probit")$g,

  ## Session parameters
  nsims     = 5,
  seed      = 123456,
  alpha     = 0.05,
  keep_raw  = TRUE
))


#------------------------------------------------------------
# Helper method for simulating group testing data using
# the groupTesting package written by Dr. Shamim Sarker
#------------------------------------------------------------

simulate_data <- function(settings, model_tag) {
  set.seed(settings$seed)

  lapply(seq_len(settings$nsims), function(i) {
    x1 <- rnorm(settings$N)
    x2 <- rbinom(settings$N, 1, 0.5)

    X <- switch(model_tag,
                "M1" = cbind(1, x1),
                "M2" = cbind(1, x1, x2),
                "M3" = cbind(1, x1, x1^2),
    )

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

    list(X = X, Z = sim_out$gtData, tsts = sim_out$testsExp)
  })
}
test_data <- simulate_data(settings, model_tag)


#----------------------------------------------
# bayesGT - Running replicates
#----------------------------------------------

res <- run_replicates(settings, test_data)
print_results(res, digits=4)

# # Saves results to a .csv file:
 save_results(res, "example")

# Diagnostic plots for single replicate (keep_raw
# must be true):
diag_histogram(res, name = "beta0", replicate = 2)
diag_beta_grid(res, replicate = 2)


# # This commented-out code is used to run a single
# # replicate rather than multiple via parallel computing.
# res <- infer_posterior(settings, test_data[[1]])
# res


## END OF FILE ##
