################################
# diagnostics.R                #
# Author: Caleb Ives           #
################################

######################
## USER FUNCTIONS   ##
######################

# Parse parameter name like "beta0", "se1", "sp2"
.parse_name <- function(name) {
  if (grepl('^beta', name)) {
    return(list(type = 'beta', j = as.integer(sub('^beta', '', name)) + 1))
  } else if (grepl('^se', name)) {
    return(list(type = 'se', j = as.integer(sub('^se', '', name))))
  } else if (grepl('^sp', name)) {
    return(list(type = 'sp', j = as.integer(sub('^sp', '', name))))
  }
  stop('Unknown parameter name: ', name)
}

# Traceplot with burn-in marker
diag_traceplot <- function(results, name, replicate = 1, main = NULL, ...) {
  if (!results$settings$keep_raw) {
    stop("Trace plots require settings$keep_raw = TRUE.")
  }

  pinfo <- .parse_name(name)
  rep_obj <- results$replicates[[replicate]]
  burn <- rep_obj$settings$burn_in

  chain <- switch(pinfo$type,
                  beta = {
                    if (is.null(rep_obj$raw$beta_all)) stop("No beta samples stored.")
                    rep_obj$raw$beta_all[, pinfo$j]
                  },
                  se = {
                    if (rep_obj$settings$known_acc) stop("SE not estimated when known_acc = TRUE.")
                    if (is.null(rep_obj$raw$se)) stop("No se samples stored.")
                    rep_obj$raw$se[, pinfo$j]
                  },
                  sp = {
                    if (rep_obj$settings$known_acc) stop("SP not estimated when known_acc = TRUE.")
                    if (is.null(rep_obj$raw$sp)) stop("No sp samples stored.")
                    rep_obj$raw$sp[, pinfo$j]
                  },
                  stop("Unknown parameter type.")
  )

  if (is.null(main)) {
    main <- sprintf("Traceplot: %s (Replication %d)", name, replicate)
  }

  oldpar <- par(mar = c(4, 4, 2, 1)); on.exit(par(oldpar))
  plot(chain,
       type = "l",
       col  = "steelblue",
       lwd  = 1.5,
       xlab = "Iteration",
       ylab = "Value",
       main = main,
       ...)
  abline(v = burn, col = "red", lty = 2, lwd = 1.5)
  invisible(chain)
}

# Posterior histogram (after burn-in)
diag_histogram <- function(results, name, replicate = 1, breaks = 30, main = NULL, ...) {
  if (!results$settings$keep_raw) {
    stop("Posterior plots require settings$keep_raw = TRUE.")
  }

  pinfo <- .parse_name(name)
  rep_obj <- results$replicates[[replicate]]
  burn <- rep_obj$settings$burn_in

  chain <- switch(pinfo$type,
                  beta = {
                    if (is.null(rep_obj$raw$beta_all)) stop("No beta samples stored.")
                    rep_obj$raw$beta_all[, pinfo$j]
                  },
                  se = {
                    if (rep_obj$settings$known_acc) stop("SE not estimated when known_acc = TRUE.")
                    if (is.null(rep_obj$raw$se)) stop("No se samples stored.")
                    rep_obj$raw$se[, pinfo$j]
                  },
                  sp = {
                    if (rep_obj$settings$known_acc) stop("SP not estimated when known_acc = TRUE.")
                    if (is.null(rep_obj$raw$sp)) stop("No sp samples stored.")
                    rep_obj$raw$sp[, pinfo$j]
                  },
                  stop("Unknown parameter type.")
  )

  post <- chain[-seq_len(burn)]

  if (is.null(main)) {
    main <- sprintf("Histogram: %s (Replication %d)", name, replicate)
  }

  oldpar <- par(mar = c(4, 4, 2, 1)); on.exit(par(oldpar))
  hist(post, breaks = breaks, freq = FALSE,
       col = "steelblue", border = "black",
       xlab = "Value", ylab = "Density",
       main = main,
       ...)
  invisible(post)
}

# Grid of beta traceplots
diag_beta_grid <- function(results, replicate = 1) {
  if (!results$settings$keep_raw) stop("Grid plot requires settings$keep_raw = TRUE.")

  rep_obj <- results$replicates[[replicate]]
  mat <- rep_obj$raw$beta_all
  burn <- rep_obj$settings$burn_in

  if (is.null(mat)) stop("No beta samples stored.")

  p <- ncol(mat)
  nc <- ceiling(sqrt(p))
  nr <- ceiling(p / nc)

  oldpar <- par(mfrow = c(nr, nc), mar = c(3, 3, 2, 1)); on.exit(par(oldpar))
  for (j in seq_len(p)) {
    plot(mat[, j], type = "l", col = "steelblue", lwd = 1.5,
         xlab = "", ylab = "",
         main = paste0("beta_", j - 1))
    abline(v = burn, col = "red", lty = 2, lwd = 1.5)
  }
  invisible(NULL)
}
