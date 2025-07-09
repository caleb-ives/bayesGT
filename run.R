#######################
# run.R               #
# Author: Caleb Ives  #
#######################

rm(list=ls())
library(bayesGT)
library(groupTesting)

###############
##  SETTINGS ##
###############
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
  nsim       = 2,
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


##########################
##  SIMULATE TEST DATA  ##
##########################
set.seed(settings$seed)
test_data <- lapply(1:settings$nsim, function(i) {
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

  list(X = X, Z = sim_out$gtData, tsts = sim_out$testsExp)
})


##########################
##  RUNNING REPLICATES  ##
##########################
res <- run_replicates(test_data, settings)

print_results(res, settings, digits=2)
#print_settings(settings)

#plot_trace(res, settings, replicate=2, parameter=1, type="se")
#plot_post_hist(res, settings, replicate=2, parameter=1, type="se")

#save_results(res, settings)


## END OF FILE ##
