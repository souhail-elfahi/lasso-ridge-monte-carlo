# Draw the figures from the saved results (no simulation is run here).
#
# Usage (from the project root):   Rscript scripts/02_make_figures.R
# Requires results/*.csv, produced by scripts/01_run_simulation.R.

if (!file.exists("results/main_grid.csv"))
  stop("Results not found. Run scripts/01_run_simulation.R from the project root first.")

dir.create("figures", showWarnings = FALSE)

grid <- read.csv("results/main_grid.csv")
path <- read.csv("results/lambda_path.csv")

# Shared style ----------------------------------------------------------------

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"; YELLOW <- "#eda100"
INK <- "#0b0b0b"; INK_SOFT <- "#52514e"; MUTED <- "#898781"
GRID_COL <- "#e1e0d9"; AXIS_COL <- "#c3c2b7"; SURFACE <- "#fcfcfb"

open_png <- function(file, width = 1500, height = 950) {
  png(file, width = width, height = height, res = 170)
  par(bg = SURFACE, fg = AXIS_COL, col.axis = MUTED, col.lab = INK_SOFT,
      col.main = INK, family = "sans", las = 1, mgp = c(2.6, 0.7, 0),
      mar = c(4.2, 4.6, 4.4, 1.2), cex.axis = 0.85, cex.lab = 0.95)
}

# Title on the first line, plain-language subtitle underneath, both left-aligned.
headline <- function(title, subtitle) {
  mtext(title, side = 3, line = 2.3, adj = 0, font = 2, cex = 1.05, col = INK)
  mtext(subtitle, side = 3, line = 1.0, adj = 0, cex = 0.82, col = INK_SOFT)
}

h_grid <- function(at) abline(h = at, col = GRID_COL, lwd = 1)

# Figure 1 - the penalty path -------------------------------------------------
# The key result: a penalty that selects almost perfectly exists, but
# cross-validation does not choose it.

x     <- log(path$lambda)
x_cv  <- log(path$lambda_cv_mean[1])
at_cv <- approx(x, path$irrelevant_kept, xout = x_cv)$y   # irrelevant kept at the CV penalty

open_png("figures/lambda_path.png")
plot(x, path$irrelevant_kept, type = "n", ylim = c(0, 5.2), axes = FALSE,
     xlab = expression(log(lambda) ~ "  (penalty: weak on the left, strong on the right)"),
     ylab = "Average number of regressors kept")
h_grid(1:5)
abline(h = 0, col = AXIS_COL)
axis(1, lwd = 0, lwd.ticks = 1, col.ticks = AXIS_COL)
axis(2, at = 0:5, col = NA, col.ticks = NA)

# Penalties that keep the 4 relevant regressors and almost no irrelevant one.
good <- path$relevant_kept >= 3.99 & path$irrelevant_kept <= 0.25
if (any(good)) {
  rect(min(x[good]), 0, max(x[good]), 5.2, col = "#2a78d61a", border = NA)
  text(mean(range(x[good])), 4.72, "Near-perfect\nselection", cex = 0.78, col = INK_SOFT)
}

segments(x_cv, 0, x_cv, 5.2, col = INK_SOFT, lwd = 1)
lines(x, path$irrelevant_kept, col = ORANGE, lwd = 2.5)
lines(x, path$relevant_kept,   col = BLUE,   lwd = 2.5)
points(x_cv, at_cv, pch = 21, bg = ORANGE, col = SURFACE, lwd = 2, cex = 1.5)

text(x_cv, 5.0, "Penalty chosen by\ncross-validation", pos = 2, cex = 0.78, col = INK_SOFT)
text(x_cv, at_cv, sprintf("%.1f irrelevant regressors\nstill in the model", at_cv),
     pos = 4, offset = 0.8, cex = 0.78, col = INK)
text(-2.4, 4.25, "Relevant (target: 4)", pos = 4, offset = 0, cex = 0.8, col = INK)
text(-1.55, 1.6, "Irrelevant (target: 0)", pos = 4, offset = 0, cex = 0.8, col = INK)
legend("right", inset = c(0.01, 0), legend = c("Relevant", "Irrelevant"),
       col = c(BLUE, ORANGE), lwd = 2.5, bty = "n", cex = 0.8, text.col = INK_SOFT)
headline("Cross-validation picks a penalty too weak to select the true model",
         "LASSO along the penalty path. General correlation 0.6, n = 200, 300 replications.")
dev.off()

# Figure 2 - selection in the three correlation cases --------------------------

sel <- rbind(grid[grid$case == "independent"    & grid$n == 500, ],
             grid[grid$case == "between_blocks" & grid$n == 500 & grid$rho == 0.6, ],
             grid[grid$case == "general"        & grid$n == 500 & grid$rho == 0.6, ])
heights <- rbind(sel$relevant_kept, sel$irrelevant_kept)

open_png("figures/selection_by_case.png")
bp <- barplot(heights, beside = TRUE, plot = FALSE, space = c(0.08, 1.1))
plot(NA, xlim = c(min(bp) - 0.9, max(bp) + 0.9), ylim = c(0, 5.3), axes = FALSE,
     xlab = "", ylab = "Average number of regressors kept", xaxs = "i")
h_grid(1:4)
barplot(heights, beside = TRUE, add = TRUE, axes = FALSE, border = NA,
        col = c(BLUE, ORANGE), space = c(0.08, 1.1))
abline(h = 0, col = AXIS_COL)
axis(2, at = 0:4, col = NA, col.ticks = NA)
axis(1, at = colMeans(bp), tick = FALSE, line = -0.4, col.axis = INK_SOFT,
     labels = c("Independent", "Between blocks (0.6)", "General (0.6)"))
text(bp, heights + 0.2, sprintf("%.2f", heights), cex = 0.8, col = INK)
legend("topright", legend = c("Relevant kept (target: 4)", "Irrelevant kept (target: 0)"),
       fill = c(BLUE, ORANGE), border = NA, bty = "n", cex = 0.8, text.col = INK_SOFT)
headline("LASSO keeps the relevant regressors but does not remove the irrelevant ones",
         "Cross-validated LASSO, n = 500, 400 replications per case.")
dev.off()

# Figure 3 - estimation error against the oracle -------------------------------

err <- grid[grid$case == "between_blocks" & grid$n == 50, ]
err <- err[order(err$rho), ]
heights <- t(as.matrix(err[, c("mse_lasso", "mse_ridge", "mse_ols_full", "mse_oracle")]))
estimator_cols <- c(BLUE, ORANGE, AQUA, YELLOW)

open_png("figures/estimation_error.png")
bp <- barplot(heights, beside = TRUE, plot = FALSE, space = c(0.08, 1.1))
plot(NA, xlim = c(min(bp) - 0.9, max(bp) + 0.9), ylim = c(0, 0.72), axes = FALSE,
     xlab = "Correlation within each relevant-irrelevant pair",
     ylab = "Mean squared error on the 4 relevant slopes", xaxs = "i")
h_grid(seq(0.1, 0.6, 0.1))
barplot(heights, beside = TRUE, add = TRUE, axes = FALSE, axisnames = FALSE,
        border = NA, col = estimator_cols, space = c(0.08, 1.1))
abline(h = 0, col = AXIS_COL)
axis(2, at = seq(0, 0.6, 0.1), col = NA, col.ticks = NA)
axis(1, at = colMeans(bp), tick = FALSE, line = -0.4, col.axis = INK_SOFT,
     labels = sprintf("%.1f", err$rho))
last <- ncol(bp)   # label only the strongest correlation, where the gap opens
text(bp[, last], heights[, last] + 0.025, sprintf("%.2f", heights[, last]), cex = 0.8, col = INK)
legend("topleft", legend = c("LASSO", "Ridge", "OLS, all candidates", "Oracle (true model known)"),
       fill = estimator_cols, border = NA, bty = "n", cex = 0.8, text.col = INK_SOFT, ncol = 2)
headline("The gap to the oracle widens as correlation increases",
         "Between-blocks correlation, n = 50, 400 replications per configuration.")
dev.off()

message("Figures saved in figures/")
