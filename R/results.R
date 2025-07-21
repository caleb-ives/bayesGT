################################
# results.R                    #
# Author: Caleb Ives           #
################################


#######################
## PRINT RESULTS     ##
#######################
print_results <- function(results, digits = 4) {
  settings <- results$settings

  if (is.null(results$summary)) {
    results$summary <- replicate_metrics(results)
  }
  summary  <- results$summary

  master_seed <- results$seed

  nsim       <- length(results$replicates)
  N          <- settings$N
  gibbs_iter <- settings$gibbs_iter
  burn_in    <- settings$burn_in
  alpha      <- settings$alpha
  mh_type    <- settings$mh_control$type
  mh_args    <- settings$mh_control$args
  keep_raw   <- settings$keep_raw
  known_acc  <- settings$known_acc

  total_runtime <- results$runtime
  avg_runtime   <- mean(summary$runtimes)

  ## Configuration table
  cat("==========================================\n")
  cat("Configuration Overview\n")
  cat("==========================================\n")

  ## Session parameters
  cat(sprintf("Master Seed         : %d\n", master_seed))
  cat(sprintf("Replicates          : %d\n", settings$nsims))
  cat(sprintf("Credible level      : %.0f%%\n", 100 * (1 - settings$alpha)))
  cat(sprintf("Keep Raw Samples    : %s\n", ifelse(settings$keep_raw, "Yes", "No")))

  ## Data structure
  cat(sprintf("Sample size (N)     : %d\n", settings$N))
  cat(sprintf("Pool sizes (psz)    : %s\n", paste(settings$psz, collapse = ", ")))
  cat(sprintf("Assay IDs           : %s\n", paste(settings$assay_id, collapse = ", ")))

  ## True parameters
  cat(sprintf("Number of betas     : %d\n", length(settings$beta_true)))
  cat(sprintf("True betas          : %s\n", paste(settings$beta_true, collapse = ", ")))
  cat(sprintf("True Se             : %s\n", paste(settings$se_true, collapse = ", ")))
  cat(sprintf("True Sp             : %s\n", paste(settings$sp_true, collapse = ", ")))

  ## Se/Sp estimation
  cat(sprintf("Known Accuracy      : %s\n", ifelse(settings$known_acc, "Yes", "No")))
  if (!settings$known_acc) {
    cat(sprintf("Initial Se          : %s\n", paste(settings$se_init, collapse = ", ")))
    cat(sprintf("Initial Sp          : %s\n", paste(settings$sp_init, collapse = ", ")))
  }

  ## Sampling configuration
  cat(sprintf("Gibbs samples       : %d\n", settings$gibbs_iter))
  cat(sprintf("Burn-in             : %d\n", settings$burn_in))
  cat(sprintf("Link function       : %s\n", settings$link_name))
  cat(sprintf("MH Type             : %s\n", settings$mh_control$type))
  cat("MH Settings         : ")
  print(settings$mh_control$args)

  cat("\n")

  ## Parameter summary table
  cat("==========================================\n")
  cat("Posterior Summary\n")
  cat("==========================================\n")

  df <- summary$param_summary
  df$TrueValue <- round(df$TrueValue, digits)
  df$EstMean   <- round(df$EstMean,   digits)
  df$CP        <- round(df$CP,        2)
  df$SSD       <- round(df$SSD,       digits)
  df$ESE       <- round(df$ESE,       digits)

  print_df <- data.frame(
    Parameter = df$Parameter,
    True      = df$TrueValue,
    Estimate  = df$EstMean,
    CP95      = df$CP,
    SSD       = df$SSD,
    ESE       = df$ESE,
    stringsAsFactors = FALSE
  )

  widths <- mapply(function(col, name) {
    # Convert NA to empty string for nchar to return 0, then take max with name length
    # Ensure a minimum width of 1 to avoid invalid 'width' argument
    max(nchar(as.character(replace(col, is.na(col), ""))), nchar(name), 1)
  }, as.list(print_df), names(print_df))

  header   <- paste(
    mapply(function(name, w) format(name, width = w, justify = "centre"),
           names(print_df), widths),
    collapse = " | "
  )
  underline <- paste(
    mapply(function(w) paste(rep("-", w), collapse = ""), widths),
    collapse = "-|-"
  )

  cat(header, "\n")
  cat(underline, "\n")

  apply(print_df, 1, function(row) {
    row <- as.list(row)
    line <- paste(
      mapply(function(cell, w) format(cell, width = w, justify = "centre"),
             row, widths),
      collapse = " | "
    )
    cat(line, "\n")
  })
  cat("\n")
  cat("\n")

  ## Runtime summary
  cat("==========================================\n")
  cat("Runtime Summary\n")
  cat("==========================================\n")
  cat(sprintf("Total runtime     : %.2f seconds\n", total_runtime))
  cat(sprintf("Avg per replicate : %.2f seconds\n", avg_runtime))
  cat(sprintf("Acceptance rate   : %.1f%%\n", 100 * mean(summary$acceptance)))
  cat("==========================================\n")
}


#######################
## PRINT RESULTS     ##
#######################
save_results <- function(results, filename = NULL) {
  settings <- results$settings
  summary  <- results$summary

  # Determine output filename
  if (is.null(filename) || !nzchar(filename)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("bayesGT_results_", timestamp, ".csv")
  }
  if (!grepl("\\.csv$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".csv")
  }

  # Construct metadata table
  meta <- data.frame(
    Key = c(
      "Runtime_sec", "Avg_Tests", "Pct_Reduction",
      "Seed", "Sample_Size", "Simulations",
      "Pool_Sizes", "Assay_IDs",
      "Known_Accuracy", "Initial_Se", "Initial_Sp",
      "True_Se", "True_Sp",
      "Gibbs_Samples", "Burn_in", "Credible_Level",
      "MH_Type", "Proposal_SD", "Keep_Raw"
    ),
    Value = c(
      summary$runtime,
      summary$avg_tests,
      sprintf("%.1f", summary$pct_reduction),
      settings$seed,
      settings$N,
      settings$nsims,
      paste(settings$psz, collapse = ";"),
      paste(settings$assay_id, collapse = ";"),
      ifelse(settings$known_acc, "Yes", "No"),
      if (!settings$known_acc) paste(settings$se_init, collapse = ";") else NA,
      if (!settings$known_acc) paste(settings$sp_init, collapse = ";") else NA,
      paste(settings$se_true, collapse = ";"),
      paste(settings$sp_true, collapse = ";"),
      settings$gibbs_iter,
      settings$burn_in,
      sprintf("%.1f%%", 100 * (1 - settings$alpha)),
      settings$mh_control$type,
      settings$mh_control$args$proposal_sd,
      ifelse(settings$keep_raw, "Yes", "No")
    ),
    stringsAsFactors = FALSE
  )

  # Format posterior parameter summary
  param_df <- summary$param_summary
  out_df <- data.frame(
    Parameter = param_df$Parameter,
    TrueValue = param_df$TrueValue,
    Estimate  = param_df$EstMean,
    Bias      = param_df$Bias,
    CP95      = param_df$CP,
    SSD       = param_df$SSD,
    ESE       = param_df$ESE,
    stringsAsFactors = FALSE
  )

  # Add Se/Sp posterior summary if available
  if (!is.null(summary$acc_summary)) {
    acc <- summary$acc_summary
    se_df <- data.frame(
      Parameter = paste0("Se", acc$Assay),
      TrueValue = acc$True_Se,
      Estimate  = acc$Est_Se,
      Bias      = acc$Bias_Se,
      CP95      = acc$CP_Se,
      SSD       = acc$SSD_Se,
      ESE       = acc$ESE_Se,
      stringsAsFactors = FALSE
    )
    sp_df <- data.frame(
      Parameter = paste0("Sp", acc$Assay),
      TrueValue = acc$True_Sp,
      Estimate  = acc$Est_Sp,
      Bias      = acc$Bias_Sp,
      CP95      = acc$CP_Sp,
      SSD       = acc$SSD_Sp,
      ESE       = acc$ESE_Sp,
      stringsAsFactors = FALSE
    )
    out_df <- rbind(out_df, se_df, sp_df)
  }

  # Write both tables to file
  con <- file(filename, open = "wt")
  on.exit(close(con), add = TRUE)

  write.table(meta, file = con,
              sep = ",", row.names = FALSE,
              col.names = TRUE, quote = TRUE)
  writeLines("", con)

  write.table(out_df, file = con,
              sep = ",", row.names = FALSE,
              col.names = TRUE, quote = TRUE)

  message("Results written to: ", normalizePath(filename))
}

## END OF FILE ##
