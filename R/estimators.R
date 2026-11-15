# Estimators ------------------------------------------------------------------
#
# Four estimators of the slopes are compared:
#   lasso    - LASSO, penalty chosen by 5-fold cross-validation (lambda.min)
#   ridge    - Ridge, penalty chosen by 5-fold cross-validation (lambda.min)
#   ols_full - OLS on all 9 candidates (the "no penalty" reference)
#   oracle   - OLS on the 4 relevant regressors only (best achievable result,
#              unavailable in practice because it assumes the true model is known)

# Slopes of the LASSO at the cross-validated penalty (intercept dropped).
fit_lasso_cv <- function(X, y, n_folds = 5) {
  as.numeric(coef(cv.glmnet(X, y, alpha = 1, nfolds = n_folds), s = "lambda.min"))[-1]
}

# Slopes of Ridge at the cross-validated penalty (intercept dropped).
fit_ridge_cv <- function(X, y, n_folds = 5) {
  as.numeric(coef(cv.glmnet(X, y, alpha = 0, nfolds = n_folds), s = "lambda.min"))[-1]
}

# Fit the four estimators on one sample. The order of the calls matters for
# reproducibility: cross-validation draws random folds, so changing the order
# would change the random numbers each estimator receives.
fit_estimators <- function(X, y) {
  lasso    <- fit_lasso_cv(X, y)
  ridge    <- fit_ridge_cv(X, y)
  ols_full <- as.numeric(coef(lm(y ~ X)))[-1]
  oracle   <- as.numeric(coef(lm(y ~ X[, RELEVANT])))[-1]
  list(lasso = lasso, ridge = ridge, ols_full = ols_full, oracle = oracle)
}

# Indices of the regressors an estimator keeps (non-zero coefficient).
selected <- function(beta_hat, tol = 1e-8) which(abs(beta_hat) > tol)
