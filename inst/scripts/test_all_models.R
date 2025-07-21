#######################
# test_all_models.R   #
# Author: Caleb Ives  #
#######################


rm(list = ls())
setwd("C:/Users/caleb/Desktop/AllModelsTest")
library(bayesGT)
library(groupTesting)


###################
#  Link function  #
###################
glmLink <- function(fn.name=c("logit","probit","cloglog")) {
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


#################
#   Constants   #
#################
BASE_SEED  <- 123
BASE_N     <- 5000
BASE_NSIM  <- 500
BASE_PSZ   <- c(5, 1)
BASE_BURN  <- 1000
BASE_ITER  <- 6000
BASE_ALPHA <- 0.05
BASE_SE0   <- 0.9
BASE_SP0   <- 0.9
BASE_RWSD  <- 0.05
LINKS      <- c("logit", "probit", "cloglog")

model_defs <- list(
  M1 = c(-3, 2),
  M2 = c(-3, 2, -1),
  M3 = c(-3, 2, -0.5)
)

configurations <- list(
  list(id = 1, known_acc = TRUE,  assay_id = c(1, 1), se_true = c(0.95, 0.95), sp_true = c(0.98, 0.98)),
  list(id = 2, known_acc = FALSE, assay_id = c(1, 1), se_true = c(0.95, 0.95), sp_true = c(0.98, 0.98)),
  list(id = 3, known_acc = FALSE, assay_id = c(1, 2), se_true = c(0.95, 0.98), sp_true = c(0.98, 0.99))
)

dir.create("results", showWarnings = FALSE)
all_results <- list()
counter <- 1
total <- length(LINKS) * length(configurations) * length(model_defs)


#################
#   Main Loop   #
#################
for (link_type in LINKS) {
  link <- glmLink(link_type)

  for (cfg in configurations) {
    for (model_name in names(model_defs)) {
      cat(sprintf("▶ [%02d/%02d] LINK: %s | TABLE: %d | MODEL: %s\n",
                  counter, total, toupper(link$name), cfg$id, model_name))

      beta_val <- model_defs[[model_name]]

      settings <- create_settings(list(
        ## Data structure
        N         = BASE_N,
        psz       = BASE_PSZ,
        assay_id  = cfg$assay_id,

        ## True parameters
        beta_true = beta_val,
        se_true   = cfg$se_true,
        sp_true   = cfg$sp_true,

        ## Se/Sp estimation
        known_acc = cfg$known_acc,
        se_init   = BASE_SE0,
        sp_init   = BASE_SP0,

        ## Sampling configuration
        gibbs_iter = BASE_ITER,
        burn_in    = BASE_BURN,
        mh_control = list(
          type = "rw",
          args = list(proposal_sd = BASE_RWSD)
        ),
        link_fcn   = link$g,
        link_name  = link$name,

        ## Session parameters
        nsims    = BASE_NSIM,
        seed     = BASE_SEED,
        alpha    = BASE_ALPHA,
        keep_raw = FALSE
      ))

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
      test_data <- simulate_data(settings, model_tag = model_name)
      results <- run_replicates(settings, test_data)
      print_results(results, digits = 3)

      fname <- sprintf("table%d_model%s_link%s", cfg$id, model_name, link$name)
      save_results(results, paste0("results/", fname))
      cat("→ Saved to: results/", fname, ".csv\n\n", sep = "")

      summary_df <- results$summary$param_summary
      summary_df$link  <- link$name
      summary_df$model <- model_name
      summary_df$table <- cfg$id
      all_results[[length(all_results) + 1]] <- summary_df

      counter <- counter + 1
    }
  }
}


################################
#   Combining/Saving Results   #
################################
combined_df <- do.call(rbind, all_results)
write.csv(combined_df, "results/all_results.csv", row.names = FALSE)
cat("All results saved to: results/all_results.csv\n")


## END OF FILE ##
