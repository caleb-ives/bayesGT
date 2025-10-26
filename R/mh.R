##############################
# mh.R                       #
# Author:  Caleb Ives        #
##############################

random_walk <- function(beta0, X, y, rw.sd, link_fcn, a = NULL, R = NULL) {
  # Propose new value
  beta.star <- beta0 + rnorm(length(beta0), 0, rw.sd)

  # Acceptance calculation
  log_alpha <- log_post(beta.star, X, y, link_fcn, a, R) -
               log_post(beta0,      X, y, link_fcn, a, R)
  alpha <- min(1, exp(log_alpha))
  accept <- rbinom(1, 1, alpha)

  # Decision
  if (accept == 1) {
    return(list(param = beta.star, accept = 1))
  } else {
    return(list(param = beta0,       accept = 0))
  }
}

random_walk_adaptive <- function(beta0, X, y, rw.sd, link_fcn,
                          target_acceptance, a = NULL, R = NULL) {
  print("Incomplete")
}

weighted_least_squares <- function(beta0, X, y, rw.sd, link_fcn, a = NULL, R = NULL) {
  print("Incomplete")
}

# Log-Posterior: Helper function for random_walk
log_post <- function(b, X, y, link_fcn, a = NULL, R = NULL) {
  XtB <- X %*% b
  p <- link_fcn(XtB)
  eps <- 1e-10

  ll <- sum(y * log(p + eps) + (1 - y) * log(1 - p + eps))

  if(!is.null(a) && !is.null(R)) {
    R.inv <- solve(R)
    lp <- -0.5 * as.numeric(t(b - a) %*% R.inv %*% (b - a))
    ll <- ll + lp
  }

  return(ll)
}
