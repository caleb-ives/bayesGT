#######################
# settings.R          #
# Author: Caleb Ives  #
#######################


# Uses default settings except for those which are
# overridden by the overrides list
create_settings <- function(overrides = list()) {
  defaults <- list(

    ## === Data structure ===
    N          = 5000,
    psz        = c(5, 1),
    assay_id   = c(1, 1),

    ## === True parameters ===
    beta_true  = c(-3, 2, -1),
    se_true    = c(0.95, 0.95),
    sp_true    = c(0.98, 0.98),

    ## === Se/Sp estimation ===
    known_acc  = FALSE,
    se_init    = 0.9,
    sp_init    = 0.9,

    ## === Sampling configuration ===
    gibbs_iter = 6000,
    burn_in    = 1000,
    mh_control = list(
      type = "rw",
      args = list(
        proposal_sd = 0.05
      )
    ),
    link_fcn   = stats::plogis,

    ## === Session parameters ===
    nsims      = 3,
    seed       = 123,
    alpha      = 0.05,
    keep_raw   = TRUE
  )

  settings <- modifyList(defaults, overrides)
  settings <- validate_settings(settings)
  return(settings)
}

validate_settings <- function(settings) {
  if (!is.numeric(settings$N) || settings$N <= 0)
    stop("Invalid setting for N")
  if (!is.numeric(settings$nsims) || settings$nsims <= 0)
    stop("Invalid setting for nsims")
  if (!is.numeric(settings$seed))
    stop("Invalid setting for seed")
  if (!is.numeric(settings$alpha) || settings$alpha <= 0 || settings$alpha >= 1)
    stop("Invalid setting for alpha")
  if (!is.numeric(settings$beta_true) || length(settings$beta_true) < 1)
    stop("Invalid setting for beta_true")
  if (!is.numeric(settings$burn_in) || settings$burn_in < 0)
    stop("Invalid setting for burn_in")
  if (!is.numeric(settings$gibbs_iter) || settings$gibbs_iter <= 0)
    stop("Invalid setting for gibbs_iter")
  if (settings$burn_in >= settings$gibbs_iter)
    stop("burn_in must be less than gibbs_iter")

  if (length(settings$psz) != length(settings$assay_id))
    stop("Length of psz must match length of assay_id")

  n_stages   <- length(settings$psz)
  assay_ids  <- settings$assay_id
  n_assays   <- max(assay_ids)

  expand_stage_param <- function(param, name) {
    if (length(param) == 1) {
      return(rep(param, n_stages))
    } else if (length(param) == n_stages) {
      return(param)
    } else {
      stop(paste("Invalid setting for", name))
    }
  }

  expand_assay_param <- function(param, name) {
    if (length(param) == 1) {
      return(rep(param, n_assays))
    } else if (length(param) == n_assays) {
      return(param)
    } else {
      stop(paste("Invalid setting for", name))
    }
  }

  # assay‑specific initial values
  settings$se_init <- expand_assay_param(settings$se_init, "se_init")
  settings$sp_init <- expand_assay_param(settings$sp_init, "sp_init")

  # stage‑specific true Se/Sp for simulation
  settings$se_true <- expand_stage_param(settings$se_true, "se_true")
  settings$sp_true <- expand_stage_param(settings$sp_true, "sp_true")

  if (!is.list(settings$mh_control) || !"type" %in% names(settings$mh_control))
    stop("Invalid setting for mh_control")
  if (!settings$mh_control$type %in% c("rw", "rw_adaptive", "wls"))
    stop("Unsupported mh_control$type")

  if (!is.function(settings$link_fcn))
    stop("Invalid setting for link_fcn")

  return(settings)
}


## END OF FILE ##
