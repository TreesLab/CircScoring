#!/usr/bin/env Rscript

# Recreate the CircScoring model-generalizability figure from the bootstrap
# result table. Only base R is required.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 02_plot_model_generalizability.R <results.csv> <output_directory>")
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

green <- "#148A2A"
orange <- "#FF6B00"
ink <- "#2A2A2A"
benchmarks <- c("P1N1", "P2N2", "P2'N2'", "P1N1", "P2N2", "P2'N2'")
scores <- c(rep("CS-R", 3), rep("CS-C", 3))
y_positions <- 6:1
y_axis_labels <- expression(P1*N1, P2*N2, P2*minute*N2*minute,
                            P1*N1, P2*N2, P2*minute*N2*minute)

open_device <- function(path, type) {
  if (type == "svg") {
    svg(path, width = 12 / 2.54, height = 7.59 / 2.54,
        bg = "transparent", family = "Arial")
  } else {
    png(path, width = 12 / 2.54, height = 7.59 / 2.54, units = "in",
        res = 600, bg = "transparent", type = "cairo", family = "Arial")
  }
}

draw_panel <- function(metric, show_y_labels) {
  plot(NA, xlim = c(0.80, 1.00), ylim = c(0.35, 6.65), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  axis(1, at = seq(0.80, 1.00, 0.05), labels = sprintf("%.2f", seq(0.80, 1.00, 0.05)),
       cex.axis = 10 / 12, lwd = 0.8, lwd.ticks = 0.8, col = ink, col.axis = ink,
       family = "Arial")
  axis(2, at = y_positions, labels = if (show_y_labels) y_axis_labels else FALSE,
       las = 1, cex.axis = 12 / 12, lwd = 0.8, lwd.ticks = 0.8,
       col = ink, col.axis = ink, family = "Arial")
  box(bty = "l", lwd = 0.8, col = ink)
  mtext(metric, side = 1, line = 2.2, cex = 12 / 12, family = "Arial", col = ink)

  for (i in seq_along(benchmarks)) {
    for (model in c("XGBoost", "Elastic Net")) {
      row <- d[d$Benchmark == benchmarks[[i]] & d$Score == scores[[i]] &
                 d$Model == model & d$Metric == metric, ]
      stopifnot(nrow(row) == 1L)
      offset <- if (model == "XGBoost") 0.15 else -0.15
      color <- if (model == "XGBoost") green else orange
      y <- y_positions[[i]] + offset
      segments(row$CI_lower, y, row$CI_upper, y, col = color, lwd = 1.1)
      segments(row$CI_lower, y - 0.10, row$CI_lower, y + 0.10, col = color, lwd = 1.1)
      segments(row$CI_upper, y - 0.10, row$CI_upper, y + 0.10, col = color, lwd = 1.1)
      points(row$Estimate, y, pch = 16, cex = 0.95, col = color)
    }
  }
}

make_plot <- function(type) {
  path <- file.path(output_dir, paste0("CircScoring_model_generalizability.", type))
  open_device(path, type)
  par(mfrow = c(1, 2), mar = c(3.2, 4.4, 0.8, 0.8), oma = c(0, 0, 4.0, 0),
      family = "Arial", fg = ink)
  draw_panel("AUROC", TRUE)
  draw_panel("AUPRC", FALSE)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend("topleft", inset = c(0.14, 0.015), bty = "n", cex = 12 / 12,
         legend = c("XGBoost", "Elastic Net"), col = c(green, orange),
         pch = 16, lty = 1, lwd = 1.1, pt.cex = 0.95,
         x.intersp = 0.7, y.intersp = 0.9, text.col = ink)
  dev.off()
}

for (type in c("svg", "png")) make_plot(type)
message("Plots written to: ", normalizePath(output_dir, mustWork = TRUE))
