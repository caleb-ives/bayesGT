################################
# estimator.R                  #
# Author: Caleb Ives           #
################################


#######################
## RUN REPLICATES    ##
#######################
run_replicates <- function(settings, test_data) {

  # Generate unique yet deterministic seeds for reproducibility
  # Each replicate has a unique seed that is based on the master seed
  seeds <- (settings$seed + seq_along(test_data)) %% .Machine$integer.max
  seeds[seeds == 0] <- 1

  cl <- parallel::makeCluster(parallel::detectCores())
  on.exit(parallel::stopCluster(cl), add=TRUE)

  runtime <- system.time({
    replicates <- parallel::parLapplyLB(
      cl,
      seq_along(test_data),
      function(i) {
        rep_settings <- settings
        rep_settings$seed <- seeds[i]
        infer_posterior(rep_settings, test_data[[i]])
      }
    )
  })["elapsed"]

  summary <- replicate_metrics(replicates)

  list(
    settings   = settings,
    replicates = replicates,
    summary    = summary,
    runtime    = as.numeric(runtime)
  )
}


#########################
## INFER POSTERIOR     ##
#########################
infer_posterior <- function(settings, test_data) {
  if (!is.null(settings$seed)) set.seed(settings$seed)

  X   <- test_data$X
  Z   <- test_data$Z[order(test_data$Z[,5]), ]
  tsts <- test_data$tsts

  assay_vec <- Z[,5]
  se_sp <- if (settings$known_acc) {
    Z[,3:4]
  } else {
    cbind(settings$se_init[assay_vec],
          settings$sp_init[assay_vec])
  }

  runtime <- system.time({
    sampler_out <- bayes_sampler(
      settings  = settings,
      beta_init = numeric(length(settings$beta_true)),
      Z         = Z,
      X         = X,
      se_sp     = se_sp
    )
  })["elapsed"]

  post <- sampler_out$param[-seq_len(settings$burn_in), , drop = FALSE]

  beta_mean <- colMeans(post)
  beta_sd   <- apply(post, 2, sd)
  ci_mat    <- apply(post, 2, quantile,
                     probs = c(settings$alpha/2, 1 - settings$alpha/2))
  cred_int  <- lapply(seq_len(ncol(ci_mat)), function(j) ci_mat[, j])

  acceptance <- sampler_out$accept_rate

  estimates <- list(
    beta_mean = beta_mean,
    beta_sd   = beta_sd,
    cred_int  = cred_int,
    tests     = tsts
  )

  raw <- NULL
  if (settings$keep_raw) {
    raw <- list(
      beta_all = sampler_out$param,
      se       = if (!settings$known_acc) sampler_out$se else NULL,
      sp       = if (!settings$known_acc) sampler_out$sp else NULL
    )
  } else if (!settings$known_acc) {
    raw <- list(
      se = sampler_out$se[-seq_len(settings$burn_in), , drop = FALSE],
      sp = sampler_out$sp[-seq_len(settings$burn_in), , drop = FALSE]
    )
  }

  return(list(
    settings   = settings,
    estimates  = estimates,
    acceptance = acceptance,
    raw        = raw,
    runtime    = as.numeric(runtime)
  ))
}


########################
## REPLICATE METRICS  ##
########################
replicate_metrics <- function(replicates) {
  # Support for passing in a single replicate
  if ("settings" %in% names(replicates) &&
      "estimates" %in% names(replicates) &&
      "acceptance" %in% names(replicates) &&
      "runtime"   %in% names(replicates)) {
    replicates <- list(replicates)
  }

  ## More support for passing in a single replicate
  if (length(replicates) == 1 && is.null(replicates[[1]]$summary)) {
    return(list(
      param_summary = data.frame(
        Parameter = paste0("beta", seq_along(replicates[[1]]$estimates$beta_mean) - 1L),
        TrueValue  = replicates[[1]]$settings$beta_true,
        EstMean    = replicates[[1]]$estimates$beta_mean,
        Bias       = replicates[[1]]$estimates$beta_mean - replicates[[1]]$settings$beta_true,
        CP         = mapply(function(ci, true)
          as.numeric(ci[1] <= true && true <= ci[2]),
          replicates[[1]]$estimates$cred_int,
          replicates[[1]]$settings$beta_true),
        SSD        = NA,
        ESE        = replicates[[1]]$estimates$beta_sd,
        stringsAsFactors = FALSE
      ),
      acc_summary   = NULL,
      acceptance    = replicates[[1]]$acceptance,
      avg_tests     = replicates[[1]]$estimates$tests,
      pct_reduction = 100 * (1 - replicates[[1]]$estimates$tests / replicates[[1]]$settings$N),
      runtimes      = c(replicates[[1]]$runtime)
    ))
  }

  ## Main function: ##

  settings   <- replicates[[1]]$settings
  beta_true  <- settings$beta_true
  alpha      <- settings$alpha

  beta_means <- do.call(rbind, lapply(replicates, function(r) r$estimates$beta_mean))
  beta_sds   <- do.call(rbind, lapply(replicates, function(r) r$estimates$beta_sd))
  tests      <- vapply(replicates, function(r) r$estimates$tests, numeric(1))
  runtimes   <- vapply(replicates, function(r) r$runtime,        numeric(1))
  cred_list  <- lapply(replicates, function(r) r$estimates$cred_int)
  acceptance <- vapply(replicates, function(r) r$acceptance, numeric(1))

  avg_tests     <- mean(tests)
  pct_reduction <- 100 * (1 - avg_tests / settings$N)

  # Beta parameter performance summary
  param_summary <- data.frame(
    Parameter = paste0("beta", seq_along(beta_true) - 1L),
    TrueValue = beta_true,
    EstMean   = colMeans(beta_means),
    Bias      = colMeans(beta_means) - beta_true,
    CP        = sapply(seq_along(beta_true), function(j) {
      ci_j <- t(vapply(cred_list,
                       function(x) x[[j]],
                       numeric(2)))
      mean(ci_j[,1] <= beta_true[j] & ci_j[,2] >= beta_true[j])
    }),
    SSD = apply(beta_means, 2, sd),
    ESE = colMeans(beta_sds),
    stringsAsFactors = FALSE
  )

  # Accuracy performance
  acc_summary <- NULL
  if (!settings$known_acc) {
    se_chains <- lapply(replicates, function(r) r$raw$se[-seq_len(settings$burn_in), , drop=FALSE])
    sp_chains <- lapply(replicates, function(r) r$raw$sp[-seq_len(settings$burn_in), , drop=FALSE])

    se_means <- do.call(rbind, lapply(se_chains, colMeans))
    sp_means <- do.call(rbind, lapply(sp_chains, colMeans))

    se_post_mean <- colMeans(se_means)
    sp_post_mean <- colMeans(sp_means)
    se_ese       <- colMeans(do.call(rbind, lapply(se_chains, function(m) apply(m, 2, sd))))
    sp_ese       <- colMeans(do.call(rbind, lapply(sp_chains, function(m) apply(m, 2, sd))))

    lb <- alpha/2; ub <- 1 - lb
    se_cp <- sapply(seq_along(se_post_mean), function(j) {
      ci_j <- t(vapply(se_chains, function(m) quantile(m[,j], c(lb, ub)), numeric(2)))
      mean(ci_j[,1] <= settings$se_true[j] & ci_j[,2] >= settings$se_true[j])
    })
    sp_cp <- sapply(seq_along(sp_post_mean), function(j) {
      ci_j <- t(vapply(sp_chains, function(m) quantile(m[,j], c(lb, ub)), numeric(2)))
      mean(ci_j[,1] <= settings$sp_true[j] & ci_j[,2] >= settings$sp_true[j])
    })

    acc_summary <- data.frame(
      Assay    = seq_along(se_post_mean),
      True_Se  = settings$se_true,
      Est_Se   = se_post_mean,
      Bias_Se  = se_post_mean - settings$se_true,
      CP_Se    = se_cp,
      SSD_Se   = apply(se_means, 2, sd),
      ESE_Se   = se_ese,
      True_Sp  = settings$sp_true,
      Est_Sp   = sp_post_mean,
      Bias_Sp  = sp_post_mean - settings$sp_true,
      CP_Sp    = sp_cp,
      SSD_Sp   = apply(sp_means, 2, sd),
      ESE_Sp   = sp_ese,
      stringsAsFactors = FALSE
    )
  }

  return(list(
    param_summary = param_summary,
    acc_summary   = acc_summary,
    acceptance    = acceptance,
    avg_tests     = avg_tests,
    pct_reduction = pct_reduction,
    runtimes      = runtimes
  ))
}


## END OF FILE ##
