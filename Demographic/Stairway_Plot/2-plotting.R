# Open TIFF device
tiff("stairway2.plot_DEJU_YEJU.tiff", height = 6.25, width = 10, units = "in", res = 300)

# Setup plot
plot(YEJU$year / 1000, YEJU$Ne_median / 1000, log = c("xy"), type = "n",
     xlab = "Time (1k years ago)", ylab = "Effective Population Size (1k individuals)",
     xlim = c(1, 1000), ylim = c(1, 1000))

# Add shaded area for DEJU confidence interval
polygon(c(DEJU$year / 1000, rev(DEJU$year / 1000)),
        c(DEJU$Ne_2.5. / 1000, rev(DEJU$Ne_97.5. / 1000)),
        col = rgb(135 / 255, 206 / 255, 235 / 255, 0.15), border = NA)

# Add lines for DEJU
lines(DEJU$year / 1000, DEJU$Ne_median / 1000, type = "s", col = "skyblue", lwd = 3, lty = 1)

# Add shaded area for YEJU confidence interval
polygon(c(YEJU$year / 1000, rev(YEJU$year / 1000)),
        c(YEJU$Ne_2.5. / 1000, rev(YEJU$Ne_97.5. / 1000)),
        col = rgb(255 / 255, 215 / 255, 0 / 255, 0.15), border = NA)

# Add lines for YEJU
lines(YEJU$year / 1000, YEJU$Ne_median / 1000, type = "s", col = "gold", lwd = 3, lty = 1)

# Add vertical shaded area for the last glacial maximum
rect(20, 0.1, 25, 1500, col = rgb(235 / 255, 105 / 255, 100 / 255, 0.05), border = NA)

# Add vertical shaded area for the penultimate glacial period
rect(135, 0.1, 194, 1500, col = rgb(169 / 255, 169 / 255, 169 / 255, 0.1), border = NA)

# Add legend in bottom right corner with adjusted parameters
legend("bottomright", legend = c("DEJU Median", "YEJU Median", "DEJU 95% CI", "YEJU 95% CI", "LGM", "PGP"),
       col = c("skyblue", "gold", rgb(135 / 255, 206 / 255, 235 / 255, 0.15), rgb(255 / 255, 215 / 255, 0 / 255, 0.15), rgb(235 / 255, 105 / 255, 100 / 255, 0.05), rgb(169 / 255, 169 / 255, 169 / 255, 0.1)),
       lty = c(1, 1, NA, NA, NA, NA),
       lwd = c(3, 3, NA, NA, NA, NA),
       pch = c(NA, NA, 15, 15, 15, 15),
       pt.cex = 2,
       cex = 0.8,
       trace = TRUE, y.intersp = 1.5)  # Adjust the font size here

# Close TIFF device
dev.off()
