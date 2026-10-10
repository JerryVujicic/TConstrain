# Run from the repository root after run_validation.R.
results <- "validation/results"
cal <- read.csv(file.path(results, "calibration.csv"))
png(file.path(results, "calibration.png"), width = 1800, height = 960, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.6, 2.2, 1.5), family = "sans", las = 1)
plot(cal$mean_probability, cal$observed_fraction, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Mean posterior edge probability", ylab = "Observed edge frequency",
     pch = 21, bg = "#176C87", col = "white", cex = 1 + sqrt(cal$n_edges) / 20,
     main = "Calibration under the model", cex.main = 1)
abline(0, 1, lty = 2, col = "#777777")
legend("topleft", "400 independent tree simulations", bty = "n", cex = 0.85)
negative <- read.csv(file.path(results, "negative_controls.csv"))
b <- barplot(negative$long_edge_probability, ylim = c(0, 1.12),
             names.arg = c("Baseline", "Correct rule", "Incorrect rule"),
             col = c("#9CA7AF", "#176C87", "#A75339"), border = NA,
             ylab = "Posterior probability of A to D",
             main = "A strong incorrect prior can dominate", cex.main = 1, cex.names = 0.85)
text(b, negative$long_edge_probability + 0.05, sprintf("%.1f%%", 100 * negative$long_edge_probability), cex = 0.9)
dev.off()
r <- read.csv(file.path(results, "rebound_posterior.csv"))
mat <- sapply(c("geometry", "soft", "soft_time"), function(m) r$posterior_probability[r$model == m])
png(file.path(results, "rebound.png"), width = 1800, height = 940, res = 200)
par(mar = c(4.7, 8, 4.5, 1), family = "sans", las = 1)
cols <- c("#9CA7AF", "#176C87", "#BC7A36")
barplot(t(mat), beside = TRUE, horiz = TRUE, xlim = c(0, 1.02),
        names.arg = paste(r$from[1:6], r$to[1:6], sep = " to "),
        col = cols, border = NA, xlab = "Posterior probability of a direct edge",
        main = "", cex.main = 1)
title("Rebound geometry and the effect of supervision", line = 3, cex.main = 1)
legend("top", inset = c(0, -0.11), xpd = NA, horiz = TRUE,
       c("Geometry", "Soft rules", "Soft rules and time"), fill = cols, border = NA, bty = "n", cex = 0.8)
dev.off()
