# Data-generating process -----------------------------------------------------
#
# True model:   y = 1 + x2 + x3 + x4 + x5 + u,   u ~ N(0, 1)
# Candidates:   a constant plus 9 random regressors, of which
#               4 are relevant (x2..x5) and 5 are irrelevant (x6..x10).
#
# The regressors are built as X = Z C', with Z a matrix of independent N(0, 1)
# factors. Their covariance matrix is therefore Sigma = C C', which is positive
# definite by construction: we choose C and read the correlation off Sigma,
# instead of fixing Sigma and hoping it is valid.

N_CANDIDATES <- 9                 # random candidate regressors (constant excluded)
RELEVANT     <- 1:4               # columns of the relevant regressors
IRRELEVANT   <- 5:9               # columns of the irrelevant regressors
BETA_TRUE    <- c(1, 1, 1, 1, 1)  # constant + the 4 relevant slopes

# Case 1 - independence: every regressor is its own factor, Sigma = I.
C_independent <- function() diag(N_CANDIDATES)

# Case 2 - correlation between blocks only: each relevant regressor shares one
# factor (weight b) with one irrelevant regressor. Correlation is b^2 within a
# pair and 0 everywhere else; the fifth irrelevant regressor stays unpaired.
C_between_blocks <- function(b) {
  C <- matrix(0, N_CANDIDATES, 4 + N_CANDIDATES)
  for (j in 1:4) { C[j, j] <- b; C[4 + j, j] <- b }
  for (i in 1:N_CANDIDATES) C[i, 4 + i] <- sqrt(max(1 - sum(C[i, 1:4]^2), 0))
  C
}

# Case 3 - general correlation: all regressors load on one common factor
# (weight a), which gives an equicorrelation matrix with correlation a^2.
C_general <- function(a) {
  C <- matrix(0, N_CANDIDATES, N_CANDIDATES + 1); C[, 1] <- a
  for (i in 1:N_CANDIDATES) C[i, i + 1] <- sqrt(1 - a^2)
  C
}

# Build C from a case label and the target correlation rho.
build_C <- function(case, rho = 0) {
  switch(case,
         independent    = C_independent(),
         between_blocks = C_between_blocks(sqrt(rho)),
         general        = C_general(sqrt(rho)),
         stop("Unknown case: ", case))
}

# Draw n observations of the regressors: a constant followed by Z C'.
generate_X <- function(n, C) {
  m <- ncol(C)
  Z <- matrix(rnorm(n * m), n, m)
  cbind(1, Z %*% t(C))
}

# Draw one full sample. Returns the 9 random regressors (no constant) and y.
simulate_sample <- function(n, C) {
  X <- generate_X(n, C)
  y <- X[, 1:5] %*% BETA_TRUE + rnorm(n)
  list(X = X[, -1], y = y)
}
