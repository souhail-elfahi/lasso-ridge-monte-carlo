# Run the full Monte Carlo study and save the results as CSV files.
#
# Usage (from the project root):   Rscript scripts/01_run_simulation.R
# or, in RStudio: open lasso-ridge-monte-carlo.Rproj, then
#                 source("scripts/01_run_simulation.R")
#
# Runtime: about 10 to 15 minutes on a laptop.
#
# The experiments run in a fixed order with a single seed. Do not reorder them:
# the results reproduce the tables of the thesis only if every random draw
# happens in the same sequence.

if (!file.exists("R/dgp.R")) stop("Run this script from the project root.")

library(glmnet)
source("R/dgp.R")
source("R/estimators.R")
source("R/experiments.R")

set.seed(20260628)
dir.create("results", showWarnings = FALSE)

rho_levels   <- c(0.3, 0.6, 0.9)
sample_sizes <- c(50, 500)

# 1. Main grid: selection and estimation error -------------------------------

run_one <- function(case, rho, n) {
  message(sprintf("  %-14s rho = %.1f  n = %d", case, rho, n))
  cbind(case = case, rho = rho, n = n, run_configuration(build_C(case, rho), n))
}

message("1/4  Main grid (14 configurations x 400 replications)")
grid <- data.frame()
for (n in sample_sizes)
  grid <- rbind(grid, run_one("independent", 0, n))
for (n in sample_sizes) for (rho in rho_levels)
  grid <- rbind(grid, run_one("between_blocks", rho, n))
for (n in sample_sizes) for (rho in rho_levels)
  grid <- rbind(grid, run_one("general", rho, n))

kept_columns <- grep("^kept_x", names(grid), value = TRUE)
write.csv(grid[, setdiff(names(grid), kept_columns)],
          "results/main_grid.csv", row.names = FALSE)
write.csv(grid[, c("case", "rho", "n", kept_columns)],
          "results/selection_frequency.csv", row.names = FALSE)

# 2. LASSO selection along the penalty path ----------------------------------

message("2/4  Penalty path (general case, rho = 0.6, n = 200)")
lambda_grid <- exp(seq(log(2), log(0.005), length.out = 40))
path <- run_lambda_path(build_C("general", 0.6), n = 200, lambda_grid)
write.csv(path, "results/lambda_path.csv", row.names = FALSE)

# 3. Retention of the paired irrelevant regressor (between-blocks case) ------

message("3/4  Paired retention (between-blocks case)")
paired <- data.frame()
for (rho in c(0.6, 0.9)) for (n in sample_sizes)
  paired <- rbind(paired, data.frame(
    rho = rho, n = n,
    paired_kept = run_paired_retention(build_C("between_blocks", rho), n)))
write.csv(paired, "results/paired_retention.csv", row.names = FALSE)

# 4. Exclusion of a relevant regressor (general case) ------------------------

message("4/4  Exclusion of a relevant regressor (general case)")
exclusion <- data.frame()
for (rho in rho_levels) for (n in sample_sizes)
  exclusion <- rbind(exclusion, data.frame(
    rho = rho, n = n,
    relevant_dropped = run_exclusion(build_C("general", rho), n)))
write.csv(exclusion, "results/exclusion.csv", row.names = FALSE)

# Record the software versions used for this run.
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

message("Done. Results saved in results/")
