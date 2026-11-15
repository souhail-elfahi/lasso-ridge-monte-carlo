# Monte Carlo experiments -----------------------------------------------------

# Squared estimation error on the 4 relevant slopes (all equal to 1).
squared_error <- function(beta_hat) sum((beta_hat[RELEVANT] - BETA_TRUE[-1])^2)

# Main experiment: selection and estimation error for one (C, n) configuration.
# Returns one row of averages over n_rep replications.
run_configuration <- function(C, n, n_rep = 400) {
  relevant_kept <- irrelevant_kept <- exact_recovery <- ridge_keeps_all <- 0
  mse <- c(lasso = 0, ridge = 0, ols_full = 0, oracle = 0)
  kept_count <- numeric(N_CANDIDATES)   # how often LASSO keeps each regressor

  for (r in 1:n_rep) {
    s   <- simulate_sample(n, C)
    fit <- fit_estimators(s$X, s$y)
    sel <- selected(fit$lasso)

    relevant_kept   <- relevant_kept   + length(intersect(sel, RELEVANT))
    irrelevant_kept <- irrelevant_kept + length(intersect(sel, IRRELEVANT))
    exact_recovery  <- exact_recovery  + setequal(sel, RELEVANT)
    ridge_keeps_all <- ridge_keeps_all + setequal(selected(fit$ridge), 1:N_CANDIDATES)
    kept_count[sel] <- kept_count[sel] + 1

    mse <- mse + sapply(fit, squared_error)
  }

  out <- data.frame(
    relevant_kept   = relevant_kept   / n_rep,   # target: 4
    irrelevant_kept = irrelevant_kept / n_rep,   # target: 0
    exact_recovery  = exact_recovery  / n_rep,   # share of runs with the true model
    ridge_keeps_all = ridge_keeps_all / n_rep,
    mse_lasso       = mse[["lasso"]]    / n_rep,
    mse_ridge       = mse[["ridge"]]    / n_rep,
    mse_ols_full    = mse[["ols_full"]] / n_rep,
    mse_oracle      = mse[["oracle"]]   / n_rep
  )
  freq <- as.data.frame(t(kept_count / n_rep))
  names(freq) <- paste0("kept_x", 2:10)
  cbind(out, freq)
}

# LASSO selection along a fixed grid of penalties, plus the penalty that
# cross-validation picks in each replication.
run_lambda_path <- function(C, n, lambda_grid, n_rep = 300) {
  relevant_kept <- irrelevant_kept <- numeric(length(lambda_grid))
  lambda_cv <- numeric(n_rep)

  for (r in 1:n_rep) {
    s <- simulate_sample(n, C)
    B <- as.matrix(glmnet(s$X, s$y, alpha = 1, lambda = lambda_grid)$beta)
    for (j in seq_along(lambda_grid)) {
      sel <- selected(B[, j])
      relevant_kept[j]   <- relevant_kept[j]   + length(intersect(sel, RELEVANT))
      irrelevant_kept[j] <- irrelevant_kept[j] + length(intersect(sel, IRRELEVANT))
    }
    lambda_cv[r] <- cv.glmnet(s$X, s$y, alpha = 1, lambda = lambda_grid,
                              nfolds = 5)$lambda.min
  }

  data.frame(lambda          = lambda_grid,
             relevant_kept   = relevant_kept   / n_rep,
             irrelevant_kept = irrelevant_kept / n_rep,
             lambda_cv_mean  = mean(lambda_cv))
}

# Between-blocks case: how often LASSO keeps the irrelevant regressor that is
# paired (correlated) with a relevant one. Averaged over the 4 pairs.
run_paired_retention <- function(C, n, n_rep = 400) {
  kept <- total <- 0
  for (r in 1:n_rep) {
    s   <- simulate_sample(n, C)
    sel <- selected(fit_lasso_cv(s$X, s$y))
    for (j in 1:4) { total <- total + 1; if ((4 + j) %in% sel) kept <- kept + 1 }
  }
  kept / total
}

# General case: how often LASSO drops at least one relevant regressor.
run_exclusion <- function(C, n, n_rep = 400) {
  dropped <- 0
  for (r in 1:n_rep) {
    s   <- simulate_sample(n, C)
    sel <- selected(fit_lasso_cv(s$X, s$y))
    dropped <- dropped + !all(RELEVANT %in% sel)
  }
  dropped / n_rep
}
