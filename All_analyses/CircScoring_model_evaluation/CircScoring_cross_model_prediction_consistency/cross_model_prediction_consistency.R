#!/usr/bin/env Rscript

# CircScoring cross-model prediction consistency analysis.
# Compares finalized XGBoost and Elastic Net scores without model refitting.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript cross_model_prediction_consistency.R <fixed_scores.tsv> <output_directory>")
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Make character tie-breaking independent of the operating-system locale.
invisible(Sys.setlocale("LC_COLLATE", "C"))
options(digits = 17)

d <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

benchmarks <- list(
  P1N1 = c("P1", "N1"),
  P2N2 = c("P2", "N2"),
  `P2'N2'` = c("P2'", "N2'")
)

score_pairs <- list(
  `CS-R` = c("XGBoost_CS-R_score", "ElasticNet_CS-R_score"),
  `CS-C` = c("XGBoost_CS-C_score", "ElasticNet_CS-C_score")
)

expected_n <- c(P1N1 = 1301L, P2N2 = 1289L, `P2'N2'` = 365L)
required_columns <- unique(c(
  "circRNA_id",
  unlist(benchmarks, use.names = FALSE),
  unlist(score_pairs, use.names = FALSE)
))
missing_columns <- setdiff(required_columns, names(d))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

deterministic_deciles <- function(score, circRNA_id) {
  if (any(!is.finite(score))) stop("Non-finite score encountered.")
  n <- length(score)
  sorted_index <- order(score, circRNA_id, method = "radix")
  decile <- integer(n)
  # Exact retained rule: deterministic cut points follow ceiling(10 * rank / n).
  decile[sorted_index] <- ceiling(10 * seq_len(n) / n)
  decile
}

quadratic_weighted_kappa <- function(x, y, categories = 1:10) {
  x <- factor(x, levels = categories)
  y <- factor(y, levels = categories)
  observed <- table(x, y)
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  index <- seq_along(categories)
  weights <- outer(index, index, function(i, j) (i - j)^2 / (length(categories) - 1)^2)
  1 - sum(weights * observed) / sum(weights * expected)
}

summary_rows <- list()
assignment_rows <- list()
matrix_rows <- list()
summary_index <- 1L
assignment_index <- 1L
matrix_index <- 1L

for (benchmark_name in names(benchmarks)) {
  group_columns <- benchmarks[[benchmark_name]]
  positive <- suppressWarnings(as.numeric(d[[group_columns[[1]]]]))
  negative <- suppressWarnings(as.numeric(d[[group_columns[[2]]]]))
  keep <- positive == 1 | negative == 1
  if (any(positive[keep] == 1 & negative[keep] == 1)) {
    stop(group_columns[[1]], " and ", group_columns[[2]], " must be mutually exclusive.")
  }
  subset <- d[keep, , drop = FALSE]
  n <- nrow(subset)
  if (n != expected_n[[benchmark_name]]) {
    stop("Unexpected sample size for ", benchmark_name, ": observed ", n,
         ", expected ", expected_n[[benchmark_name]], ".")
  }

  for (score_name in names(score_pairs)) {
    columns <- score_pairs[[score_name]]
    xgboost_score <- as.numeric(subset[[columns[[1]]]])
    elastic_score <- as.numeric(subset[[columns[[2]]]])

    correlation <- suppressWarnings(cor.test(
      xgboost_score,
      elastic_score,
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    ))

    xgboost_decile <- deterministic_deciles(xgboost_score, subset$circRNA_id)
    elastic_decile <- deterministic_deciles(elastic_score, subset$circRNA_id)
    kappa <- quadratic_weighted_kappa(xgboost_decile, elastic_decile)

    summary_rows[[summary_index]] <- data.frame(
      benchmark = benchmark_name,
      score = score_name,
      n = n,
      spearman_rho = unname(correlation$estimate),
      spearman_p = correlation$p.value,
      quadratic_weighted_kappa = kappa,
      x_axis = "XGBoost decile",
      y_axis = "Elastic Net decile",
      decile_method = paste0(
        "ascending score; ties resolved by circRNA_id ascending; ",
        "equal-count deterministic deciles"
      ),
      stringsAsFactors = FALSE
    )
    summary_index <- summary_index + 1L

    assignment_rows[[assignment_index]] <- data.frame(
      circRNA_id = subset$circRNA_id,
      benchmark = benchmark_name,
      score = score_name,
      XGBoost_score = xgboost_score,
      ElasticNet_score = elastic_score,
      XGBoost_decile = xgboost_decile,
      ElasticNet_decile = elastic_decile,
      stringsAsFactors = FALSE
    )
    assignment_index <- assignment_index + 1L

    matrix_count <- table(
      factor(xgboost_decile, levels = 1:10),
      factor(elastic_decile, levels = 1:10)
    )
    matrix_rows[[matrix_index]] <- data.frame(
      benchmark = benchmark_name,
      score = score_name,
      XGBoost_decile = rep(1:10, times = 10),
      ElasticNet_decile = rep(1:10, each = 10),
      count = as.vector(matrix_count),
      percent_of_BSJs = 100 * as.vector(matrix_count) / n,
      stringsAsFactors = FALSE
    )
    matrix_index <- matrix_index + 1L
  }
}

summary_result <- do.call(rbind, summary_rows)
assignment_result <- do.call(rbind, assignment_rows)
matrix_result <- do.call(rbind, matrix_rows)

write.csv(
  summary_result,
  file.path(output_dir, "Cross_model_consistency_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.table(
  assignment_result,
  file.path(output_dir, "Cross_model_decile_assignments.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  matrix_result,
  file.path(output_dir, "Cross_model_10x10_matrices.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Analysis completed. Outputs written to: ", normalizePath(output_dir, mustWork = TRUE))
