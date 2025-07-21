#########################################
# bayes_sampler.R                       #
# Author:  Caleb Ives                   #
# Purpose: Runs one replicate of        #
# the Bayesian Gibbs sampler            #
#########################################


bayes_sampler <- function(settings, beta_init, Z, X, se_sp){
  N          <- settings$N
  S          <- length(settings$psz)
  link_fcn   <- settings$link_fcn
  rw.sd      <- settings$mh_control$args$proposal_sd
  post_git   <- settings$gibbs_iter
  known_acc  <- settings$known_acc
  se_0       <- settings$se_init
  sp_0       <- settings$sp_init

  mh <- switch(
    settings$mh_control$type,
    rw          = random_walk,
    rw_adaptive = random_walk_adaptive,
    wls         = weighted_least_squares
  )

  X  <- as.matrix(X)
  Yt <- rbinom(N, 1, link_fcn(X %*% beta_init))

  # Sort and extract assay info
  Z             <- Z[order(Z[,5]), ]
  assay_vec     <- Z[,5]
  unique_assays <- sort(unique(assay_vec))
  num_assays    <- length(unique_assays)

  # Construct Y‑til matrix
  tmp <- Z[, -(1:5)]
  ytm <- matrix(-9L, N, S)
  for (d in seq_len(N)) {
    idx <- which(tmp == d, arr.ind = TRUE)[, "row"]
    if (length(idx) > 0) ytm[d, seq_along(idx)] <- sort(idx)
  }
  Ytmat <- cbind(
    Yt,
    rowSums(ytm > 0),
    ytm
  )
  Ycol <- ncol(Ytmat)

  # Prepare Fortran inputs
  Z_obs     <- Z[, -(3:5)]
  Zrow      <- nrow(Z_obs)
  Zcol      <- ncol(Z_obs)
  Ztil_out  <- integer(Zrow * Zcol)
  se_out    <- integer(2L * num_assays)
  sp_out    <- integer(2L * num_assays)
  U_all     <- matrix(runif(N * post_git), nrow = N, ncol = post_git)
  assay_map <- match(assay_vec, unique_assays)

  # Pre-allocate
  beta_sv <- matrix(NA_real_, nrow = post_git, ncol = length(beta_init))
  accept  <- integer(post_git)

  if (!known_acc) {
    se_sv <- matrix(NA_real_, nrow = post_git, ncol = num_assays)
    sp_sv <- matrix(NA_real_, nrow = post_git, ncol = num_assays)
    se    <- rep(se_0, length.out = num_assays)
    sp    <- rep(sp_0, length.out = num_assays)
  } else {
    se <- se_sp[,1]
    sp <- se_sp[,2]
  }


  # Begin Gibbs sampling:
  for (s in 1:post_git) {
    # --- Sample Ytil
    res <- .Fortran("sample",
                    p         = as.double(link_fcn(X %*% beta_init)),
                    Ytmat     = as.integer(Ytmat),
                    Z_mat     = as.integer(Z_obs),
                    N         = as.integer(N),
                    Ycols     = as.integer(Ycol),
                    Zrows     = as.integer(Zrow),
                    Zcols     = as.integer(Zcol),
                    U         = as.double(U_all[, s]),
                    Ztil      = as.integer(Ztil_out),
                    assay_vec = as.integer(assay_vec),
                    L         = as.integer(num_assays),
                    se_in     = as.double(se),
                    sp_in     = as.double(sp),
                    se_counts = as.integer(se_out),
                    sp_counts = as.integer(sp_out)
    )
    Ytmat     <- matrix(res$Ytmat, N, Ycol)
    se_counts <- matrix(res$se_counts, nrow = 2)
    sp_counts <- matrix(res$sp_counts, nrow = 2)

    # --- Sample beta
    mh_out        <- mh(beta_init, X, Ytmat[,1], rw.sd, link_fcn)
    beta_init     <- mh_out$param
    beta_sv[s, ]  <- beta_init
    accept[s]     <- mh_out$accept

    # Gibbs sampling for Se/Sp
    if (!known_acc) {
      for (a in seq_len(num_assays)) {
        tp <- se_counts[1, a]; fn <- se_counts[2, a]
        tn <- sp_counts[1, a]; fp <- sp_counts[2, a]
        se[a] <- rbeta(1, 1 + tp, 1 + fn)
        sp[a] <- rbeta(1, 1 + tn, 1 + fp)
      }
      se_sv[s, ] <- se
      sp_sv[s, ] <- sp
      se_sp[,1]  <- se[assay_map]
      se_sp[,2]  <- sp[assay_map]
    }

  }

  output <- list(
    param       = beta_sv,
    convergence = 0,
    accept      = accept,
    accept_rate = mean(accept)
  )
  if (!known_acc) {
    colnames(se_sv) <- paste0("Se_a", unique_assays)
    colnames(sp_sv) <- paste0("Sp_a", unique_assays)
    output$se <- se_sv
    output$sp <- sp_sv
  }
  return(output)
}


## END OF FILE ##
