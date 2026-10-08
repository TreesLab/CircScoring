#!/usr/bin/env Rscript

# Recreate the model-performance stability violin/box plots from the retained
# 1,000-partition result table. Only base R is required.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 02_plot_model_performance_stability.R <1000_partitions.csv> <output_directory>")
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(input_csv, check.names = FALSE)

green <- "#2E8B57"
orange <- "#EF7D17"
ink <- "#111111"
models <- c("XGBoost_CS-R", "XGBoost_CS-C", "ElasticNet_CS-R", "ElasticNet_CS-C")
positions <- c(1.0, 1.35, 2.1, 2.45, 3.5, 3.85, 4.6, 4.95)

open_device <- function(path, type) {
  if (type == "svg") {
    svg(path, width = 15 / 2.54, height = 7 / 2.54, bg = "transparent", family = "Arial")
  } else {
    png(path, width = 15 / 2.54, height = 7 / 2.54, units = "in", res = 600,
        bg = "transparent", type = "cairo", family = "Arial")
  }
}

draw_violin <- function(x, at, color, filled, width = 0.14) {
  den <- density(x, na.rm = TRUE, n = 512, adjust = 1)
  scale <- width / max(den$y)
  polygon(c(at - den$y * scale, rev(at + den$y * scale)),
          c(den$x, rev(den$x)),
          col = if (filled) color else "white", border = color, lwd = 1.2)
  boxplot(x, at = at, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.10,
          border = ink, col = "white", staplewex = 0.55, whisklty = 1, lwd = 0.8)
}

plot_metric <- function(metric, ylim, ylab) {
  plot(NA, xlim = c(0.65, 5.3), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  axis(2, las = 1, cex.axis = 10 / 12, family = "Arial", col = ink, col.axis = ink)
  axis(1, at = c(1.175, 2.275, 3.675, 4.775), labels = FALSE, tck = 0)
  box(bty = "l", col = ink)
  mtext(ylab, side = 2, line = 2.5, cex = 12 / 12, family = "Arial")
  for (i in seq_along(models)) {
    model <- models[[i]]
    color <- if (grepl("XGBoost", model)) green else orange
    x80 <- d[[paste0(model, "_", metric, "_80")]]
    x20 <- d[[paste0(model, "_", metric, "_20")]]
    draw_violin(x80, positions[2 * i - 1], color, TRUE)
    draw_violin(x20, positions[2 * i], color, FALSE)
  }
}

make_main <- function(type) {
  ext <- type
  path <- file.path(output_dir, paste0("model_performance_stability_AUROC_AUPRC.", ext))
  open_device(path, type)
  par(mfrow = c(1, 2), mar = c(1.7, 3.5, 2.5, 0.6), oma = c(0, 0, 2.2, 0),
      family = "Arial", fg = ink)
  plot_metric("AUROC", c(0.80, 1.00), "AUROC")
  title("AUROC", cex.main = 12 / 12, family = "Arial")
  plot_metric("AUPRC", c(0.80, 1.00), "AUPRC")
  title("AUPRC", cex.main = 12 / 12, family = "Arial")
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend("top", inset = 0.015, horiz = TRUE, bty = "n", cex = 9 / 12,
         legend = c("XGBoost 80%", "XGBoost 20%", "Elastic Net 80%", "Elastic Net 20%"),
         fill = c(green, "white", orange, "white"),
         border = c(green, green, orange, orange), x.intersp = 0.5, text.width = 0.14)
  dev.off()
}

plot_difference_metric <- function(metric, ylab) {
  vals <- lapply(models, function(m) d[[paste0(m, "_", metric, "_diff_80minus20")]])
  y_lim <- range(unlist(vals), finite = TRUE)
  pad <- diff(y_lim) * 0.07
  plot(NA, xlim = c(0.5, 4.5), ylim = y_lim + c(-pad, pad), axes = FALSE, xlab = "", ylab = "")
  abline(h = 0, lty = 2, col = "#666666", lwd = 0.8)
  axis(2, las = 1, cex.axis = 10 / 12, family = "Arial", col = ink, col.axis = ink)
  axis(1, at = 1:4, labels = FALSE, tck = 0)
  box(bty = "l", col = ink)
  mtext(ylab, side = 2, line = 2.5, cex = 12 / 12, family = "Arial")
  for (i in seq_along(models)) {
    color <- if (grepl("XGBoost", models[[i]])) green else orange
    draw_violin(vals[[i]], i, color, TRUE, width = 0.22)
  }
}

make_difference <- function(type) {
  path <- file.path(output_dir, paste0("model_performance_stability_difference_80minus20.", type))
  open_device(path, type)
  par(mfrow = c(1, 2), mar = c(1.7, 3.5, 2.5, 0.6), family = "Arial", fg = ink)
  plot_difference_metric("AUROC", expression(Delta * "AUROC (80% - 20%)"))
  title(expression(Delta * "AUROC"), cex.main = 12 / 12, family = "Arial")
  plot_difference_metric("AUPRC", expression(Delta * "AUPRC (80% - 20%)"))
  title(expression(Delta * "AUPRC"), cex.main = 12 / 12, family = "Arial")
  dev.off()
}

for (type in c("svg", "png")) {
  make_main(type)
  make_difference(type)
}

message("Plots written to: ", normalizePath(output_dir, mustWork = TRUE))
