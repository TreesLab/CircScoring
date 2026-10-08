#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, digits = 17)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
} else {
  normalizePath("CircScoring_model_development.R", mustWork = TRUE)
}
project_dir <- dirname(script_path)

get_option_value <- function(arguments, name, default = NULL) {
  location <- match(name, arguments)
  if (is.na(location)) return(default)
  if (location == length(arguments)) stop("Missing value after ", name)
  arguments[[location + 1L]]
}

arguments <- commandArgs(trailingOnly = TRUE)
input_file <- get_option_value(
  arguments,
  "--input",
  file.path(project_dir, "data", "P1N1_training_data.tsv")
)
output_dir <- get_option_value(
  arguments,
  "--output",
  file.path(project_dir, "results")
)
seed <- as.integer(get_option_value(arguments, "--seed", "123"))
threads <- as.integer(get_option_value(arguments, "--threads", "1"))
quick_test <- "--quick-test" %in% arguments
refit_only <- "--refit-only" %in% arguments

if (!is.finite(seed)) stop("--seed must be an integer")
if (!is.finite(threads) || threads < 1L) stop("--threads must be >= 1")
if (quick_test && refit_only) stop("Use either --quick-test or --refit-only, not both")

required_packages <- c("xgboost", "glmnet")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run install_exact_packages.R first."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_dir <- file.path(output_dir, "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

feature_map <- data.frame(
  feature_id = paste0("f", 1:18),
  feature_name = c(
    "donor_site_at_the_annotated_boundary",
    "acceptor_site_at_the_annotated_boundary",
    "donor_acceptor_sites_at_the_same_transcript_isoform",
    "has_AS_event(donor)",
    "has_AS_event(acceptor)",
    "phyloP(acceptor_in)",
    "phyloP(acceptor_out)",
    "phyloP(donor_in)",
    "phyloP(donor_out)",
    "phastCons(acceptor_in)",
    "phastCons(acceptor_out)",
    "phastCons(donor_in)",
    "phastCons(donor_out)",
    "MAXENT(donor)",
    "MAXENT(acceptor)",
    "7 DB",
    "FL-circAS or circFL_seq",
    "mouse_conserved"
  ),
  mode = c(rep("CS-C", 15), rep("CS-R", 3)),
  check.names = FALSE
)
mode_features <- split(feature_map$feature_name, feature_map$mode)

expected_best <- data.frame(
  method = c("XGBoost", "XGBoost", "Elastic Net", "Elastic Net"),
  mode = c("CS-R", "CS-C", "CS-R", "CS-C"),
  max_depth = c(4, 2, NA, NA),
  learning_rate = c(0.03, 0.05, NA, NA),
  min_child_weight = c(1, 5, NA, NA),
  subsample = c(0.7, 0.7, NA, NA),
  colsample_bytree = c(0.8, 1.0, NA, NA),
  gamma = c(0.0, 0.5, NA, NA),
  nrounds = c(100, 100, NA, NA),
  alpha = c(NA, NA, 0.3, 0.1),
  lambda = c(NA, NA, 0.18420699693267201, 0.0001),
  check.names = FALSE
)

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

read_training_data <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  data <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "")
  )
  required <- c("circRNA_id", "P1", "N1", feature_map$feature_name)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing input column(s): ", paste(missing_columns, collapse = ", "))
  }
  data <- data[, required, drop = FALSE]
  if (anyNA(data$circRNA_id) || any(!nzchar(trimws(data$circRNA_id)))) {
    stop("circRNA_id cannot be missing or blank")
  }
  if (anyDuplicated(data$circRNA_id)) stop("P1/N1 circRNA_id values must be unique")
  data$P1 <- as.numeric(data$P1)
  data$N1 <- as.numeric(data$N1)
  valid_class <- (data$P1 == 1 & data$N1 == 0) |
    (data$P1 == 0 & data$N1 == 1)
  if (anyNA(valid_class) || any(!valid_class)) {
    stop("Every row must be exclusively P1-positive or N1-negative")
  }
  for (feature_name in feature_map$feature_name) {
    data[[feature_name]] <- as.numeric(data[[feature_name]])
    if (anyNA(data[[feature_name]]) || any(!is.finite(data[[feature_name]]))) {
      stop("Predictor contains missing or non-finite values: ", feature_name)
    }
  }
  label <- as.integer(data$P1 == 1)
  if (nrow(data) != 1301L || sum(label == 1L) != 797L || sum(label == 0L) != 504L) {
    stop(
      "Unexpected development population: rows=", nrow(data),
      ", P1=", sum(label == 1L), ", N1=", sum(label == 0L)
    )
  }
  list(data = data, label = label)
}

make_stratified_folds <- function(label, k = 5L, random_seed = 123L) {
  set.seed(random_seed)
  fold <- integer(length(label))
  for (class_value in sort(unique(label))) {
    indices <- which(label == class_value)
    indices <- sample(indices, length(indices), replace = FALSE)
    fold[indices] <- rep(seq_len(k), length.out = length(indices))
  }
  fold
}

calculate_auroc <- function(label, probability) {
  if (length(label) != length(probability) || anyNA(probability)) {
    stop("Invalid AUROC inputs")
  }
  positive_n <- sum(label == 1L)
  negative_n <- sum(label == 0L)
  if (!positive_n || !negative_n) return(NA_real_)
  probability_rank <- rank(probability, ties.method = "average")
  (
    sum(probability_rank[label == 1L]) - positive_n * (positive_n + 1) / 2
  ) / (positive_n * negative_n)
}

calculate_auprc <- function(label, probability) {
  if (length(label) != length(probability) || anyNA(probability)) {
    stop("Invalid AUPRC inputs")
  }
  positive_n <- sum(label == 1L)
  if (!positive_n) return(NA_real_)
  order_index <- order(probability, decreasing = TRUE)
  ordered_probability <- probability[order_index]
  ordered_label <- label[order_index]
  group_end <- c(which(diff(ordered_probability) != 0), length(ordered_probability))
  cumulative_tp <- cumsum(ordered_label == 1L)[group_end]
  cumulative_fp <- cumsum(ordered_label == 0L)[group_end]
  recall <- cumulative_tp / positive_n
  precision <- cumulative_tp / (cumulative_tp + cumulative_fp)
  sum(diff(c(0, recall)) * precision)
}

metric_summary <- function(label, probability) {
  c(AUROC = calculate_auroc(label, probability), AUPRC = calculate_auprc(label, probability))
}

as_predictor_matrix <- function(data, features) {
  matrix_data <- data.matrix(data[, features, drop = FALSE])
  storage.mode(matrix_data) <- "double"
  matrix_data
}

fit_xgboost <- function(x, y, parameters, random_seed, nthread) {
  xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = as.integer(parameters$max_depth),
      eta = as.numeric(parameters$learning_rate),
      min_child_weight = as.numeric(parameters$min_child_weight),
      subsample = as.numeric(parameters$subsample),
      colsample_bytree = as.numeric(parameters$colsample_bytree),
      gamma = as.numeric(parameters$gamma),
      seed = as.integer(random_seed),
      nthread = as.integer(nthread),
      verbosity = 0
    ),
    data = xgboost::xgb.DMatrix(x, label = y, missing = NA_real_),
    nrounds = as.integer(parameters$nrounds),
    verbose = 0
  )
}

predict_xgboost <- function(model, x) {
  as.numeric(stats::predict(model, xgboost::xgb.DMatrix(x, missing = NA_real_)))
}

fit_elastic_net <- function(x, y, parameters) {
  glmnet::glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = as.numeric(parameters$alpha),
    lambda = as.numeric(parameters$lambda),
    intercept = TRUE,
    standardize = TRUE,
    thresh = 1e-7
  )
}

predict_elastic_net <- function(model, x, lambda) {
  as.numeric(stats::predict(model, newx = x, s = lambda, type = "response"))
}

summarize_fold_metrics <- function(fold_metrics) {
  c(
    CV_AUROC_mean = mean(fold_metrics[, "AUROC"]),
    CV_AUROC_sd = stats::sd(fold_metrics[, "AUROC"]),
    CV_AUPRC_mean = mean(fold_metrics[, "AUPRC"]),
    CV_AUPRC_sd = stats::sd(fold_metrics[, "AUPRC"])
  )
}

tune_xgboost <- function(data, label, fold, grid, random_seed, nthread) {
  rows <- vector("list", nrow(grid) * 2L)
  output_index <- 0L
  for (mode in c("CS-R", "CS-C")) {
    x <- as_predictor_matrix(data, mode_features[[mode]])
    for (candidate_index in seq_len(nrow(grid))) {
      parameters <- grid[candidate_index, , drop = FALSE]
      fold_metrics <- matrix(
        NA_real_, nrow = max(fold), ncol = 2,
        dimnames = list(NULL, c("AUROC", "AUPRC"))
      )
      start_time <- proc.time()[[3]]
      for (fold_number in seq_len(max(fold))) {
        validation <- fold == fold_number
        model <- fit_xgboost(
          x[!validation, , drop = FALSE],
          label[!validation],
          parameters,
          random_seed,
          nthread
        )
        probability <- predict_xgboost(model, x[validation, , drop = FALSE])
        fold_metrics[fold_number, ] <- metric_summary(label[validation], probability)
      }
      output_index <- output_index + 1L
      rows[[output_index]] <- data.frame(
        method = "XGBoost",
        mode = mode,
        candidate_id = candidate_index,
        parameters,
        as.list(summarize_fold_metrics(fold_metrics)),
        elapsed_seconds = proc.time()[[3]] - start_time,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (candidate_index %% 25L == 0L || candidate_index == nrow(grid)) {
        message("XGBoost ", mode, ": ", candidate_index, "/", nrow(grid))
      }
    }
  }
  result <- do.call(rbind, rows[seq_len(output_index)])
  result$rank <- ave(
    seq_len(nrow(result)), result$mode,
    FUN = function(index) {
      ranking <- order(
        -result$CV_AUROC_mean[index],
        -result$CV_AUPRC_mean[index],
        result$candidate_id[index]
      )
      rank_value <- integer(length(index))
      rank_value[ranking] <- seq_along(ranking)
      rank_value
    }
  )
  result
}

tune_elastic_net <- function(data, label, fold, grid) {
  rows <- vector("list", nrow(grid) * 2L)
  output_index <- 0L
  for (mode in c("CS-R", "CS-C")) {
    x <- as_predictor_matrix(data, mode_features[[mode]])
    for (candidate_index in seq_len(nrow(grid))) {
      parameters <- grid[candidate_index, , drop = FALSE]
      fold_metrics <- matrix(
        NA_real_, nrow = max(fold), ncol = 2,
        dimnames = list(NULL, c("AUROC", "AUPRC"))
      )
      start_time <- proc.time()[[3]]
      for (fold_number in seq_len(max(fold))) {
        validation <- fold == fold_number
        model <- fit_elastic_net(
          x[!validation, , drop = FALSE],
          label[!validation],
          parameters
        )
        probability <- predict_elastic_net(
          model, x[validation, , drop = FALSE], parameters$lambda
        )
        fold_metrics[fold_number, ] <- metric_summary(label[validation], probability)
      }
      output_index <- output_index + 1L
      rows[[output_index]] <- data.frame(
        method = "Elastic Net",
        mode = mode,
        candidate_id = candidate_index,
        parameters,
        as.list(summarize_fold_metrics(fold_metrics)),
        elapsed_seconds = proc.time()[[3]] - start_time,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      if (candidate_index %% 25L == 0L || candidate_index == nrow(grid)) {
        message("Elastic Net ", mode, ": ", candidate_index, "/", nrow(grid))
      }
    }
  }
  result <- do.call(rbind, rows[seq_len(output_index)])
  result$rank <- ave(
    seq_len(nrow(result)), result$mode,
    FUN = function(index) {
      ranking <- order(
        -result$CV_AUROC_mean[index],
        -result$CV_AUPRC_mean[index],
        result$candidate_id[index]
      )
      rank_value <- integer(length(index))
      rank_value[ranking] <- seq_along(ranking)
      rank_value
    }
  )
  result
}

selected_from_tuning <- function(xgboost_tuning, elastic_tuning) {
  selected <- list()
  for (mode in c("CS-R", "CS-C")) {
    xgb_row <- xgboost_tuning[xgboost_tuning$mode == mode & xgboost_tuning$rank == 1, ]
    en_row <- elastic_tuning[elastic_tuning$mode == mode & elastic_tuning$rank == 1, ]
    selected[[paste("XGBoost", mode)]] <- xgb_row[1, , drop = FALSE]
    selected[[paste("Elastic Net", mode)]] <- en_row[1, , drop = FALSE]
  }
  selected
}

selected_from_expected <- function() {
  selected <- list()
  for (row_index in seq_len(nrow(expected_best))) {
    row <- expected_best[row_index, , drop = FALSE]
    selected[[paste(row$method, row$mode)]] <- row
  }
  selected
}

make_parameter_summary <- function(selected) {
  rows <- lapply(selected, function(row) {
    wanted <- c(
      "method", "mode", "max_depth", "learning_rate", "min_child_weight",
      "subsample", "colsample_bytree", "gamma", "nrounds", "alpha", "lambda",
      "CV_AUROC_mean", "CV_AUROC_sd", "CV_AUPRC_mean", "CV_AUPRC_sd"
    )
    for (name in setdiff(wanted, names(row))) row[[name]] <- NA
    row[, wanted, drop = FALSE]
  })
  do.call(rbind, rows)
}

generate_selected_oof <- function(data, label, fold, selected, random_seed, nthread) {
  output <- data.frame(
    circRNA_id = data$circRNA_id,
    label = label,
    fold = fold,
    stringsAsFactors = FALSE
  )
  performance <- list()
  for (method in c("XGBoost", "Elastic Net")) {
    for (mode in c("CS-R", "CS-C")) {
      parameters <- selected[[paste(method, mode)]]
      x <- as_predictor_matrix(data, mode_features[[mode]])
      probability <- rep(NA_real_, nrow(data))
      for (fold_number in seq_len(max(fold))) {
        validation <- fold == fold_number
        if (method == "XGBoost") {
          model <- fit_xgboost(
            x[!validation, , drop = FALSE], label[!validation], parameters,
            random_seed, nthread
          )
          probability[validation] <- predict_xgboost(
            model, x[validation, , drop = FALSE]
          )
        } else {
          model <- fit_elastic_net(
            x[!validation, , drop = FALSE], label[!validation], parameters
          )
          probability[validation] <- predict_elastic_net(
            model, x[validation, , drop = FALSE], parameters$lambda
          )
        }
      }
      column_name <- paste(method, mode, "OOF_probability", sep = "_")
      output[[column_name]] <- probability
      metrics <- metric_summary(label, probability)
      performance[[paste(method, mode)]] <- data.frame(
        method = method,
        mode = mode,
        OOF_AUROC = metrics[["AUROC"]],
        OOF_AUPRC = metrics[["AUPRC"]],
        stringsAsFactors = FALSE
      )
    }
  }
  list(predictions = output, performance = do.call(rbind, performance))
}

refit_all_models <- function(data, label, selected, random_seed, nthread) {
  fitted_models <- list()
  fitted_scores <- data.frame(
    circRNA_id = data$circRNA_id,
    label = label,
    stringsAsFactors = FALSE
  )
  xgb_importance <- list()
  elastic_coefficients <- list()

  for (method in c("XGBoost", "Elastic Net")) {
    for (mode in c("CS-R", "CS-C")) {
      key <- paste(method, mode)
      parameters <- selected[[key]]
      features <- mode_features[[mode]]
      x <- as_predictor_matrix(data, features)
      feature_ids <- feature_map$feature_id[match(features, feature_map$feature_name)]
      if (method == "XGBoost") {
        model <- fit_xgboost(x, label, parameters, random_seed, nthread)
        probability <- predict_xgboost(model, x)
        xgboost::xgb.save(
          model,
          file.path(model_dir, paste0("XGBoost_", mode, ".ubj"))
        )
        saveRDS(model, file.path(model_dir, paste0("XGBoost_", mode, ".rds")))
        importance <- xgboost::xgb.importance(model = model)
        all_importance <- data.frame(
          Feature = features,
          Gain = 0,
          Cover = 0,
          Frequency = 0,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        if (nrow(importance)) {
          matched <- match(importance$Feature, all_importance$Feature)
          all_importance$Gain[matched] <- importance$Gain
          all_importance$Cover[matched] <- importance$Cover
          all_importance$Frequency[matched] <- importance$Frequency
        }
        all_importance$method <- method
        all_importance$mode <- mode
        all_importance$feature_id <- feature_ids
        all_importance$feature_name <- features
        xgb_importance[[mode]] <- all_importance[, c(
          "method", "mode", "feature_id", "feature_name", "Gain", "Cover", "Frequency"
        )]
      } else {
        model <- fit_elastic_net(x, label, parameters)
        probability <- predict_elastic_net(model, x, parameters$lambda)
        saveRDS(model, file.path(model_dir, paste0("ElasticNet_", mode, ".rds")))
        coefficient_matrix <- as.matrix(stats::coef(model, s = parameters$lambda))
        coefficient_names <- rownames(coefficient_matrix)
        ids <- c(NA_character_, feature_ids)
        names(ids) <- c("(Intercept)", features)
        elastic_coefficients[[mode]] <- data.frame(
          method = method,
          mode = mode,
          feature_id = unname(ids[coefficient_names]),
          feature_name = coefficient_names,
          effect_type = ifelse(coefficient_names == "(Intercept)", "intercept", "coefficient"),
          coefficient = as.numeric(coefficient_matrix[, 1]),
          alpha = as.numeric(parameters$alpha),
          lambda = as.numeric(parameters$lambda),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
      fitted_models[[key]] <- model
      fitted_scores[[paste(method, mode, "score", sep = "_")]] <- probability
    }
  }
  list(
    models = fitted_models,
    fitted_scores = fitted_scores,
    xgboost_importance = do.call(rbind, xgb_importance),
    elastic_coefficients = do.call(rbind, elastic_coefficients)
  )
}

compare_with_expected <- function(best_parameters, elastic_coefficients) {
  checks <- list()
  check_index <- 0L
  for (row_index in seq_len(nrow(expected_best))) {
    expected <- expected_best[row_index, , drop = FALSE]
    actual <- best_parameters[
      best_parameters$method == expected$method & best_parameters$mode == expected$mode,
      , drop = FALSE
    ]
    parameter_names <- if (expected$method == "XGBoost") {
      c("max_depth", "learning_rate", "min_child_weight", "subsample", "colsample_bytree", "gamma", "nrounds")
    } else {
      c("alpha", "lambda")
    }
    for (parameter_name in parameter_names) {
      check_index <- check_index + 1L
      difference <- abs(as.numeric(actual[[parameter_name]][1]) - as.numeric(expected[[parameter_name]][1]))
      checks[[check_index]] <- data.frame(
        check = paste(expected$method, expected$mode, parameter_name),
        expected = as.numeric(expected[[parameter_name]][1]),
        actual = as.numeric(actual[[parameter_name]][1]),
        absolute_difference = difference,
        tolerance = 1e-12,
        status = ifelse(is.finite(difference) && difference <= 1e-12, "PASS", "FAIL"),
        stringsAsFactors = FALSE
      )
    }
  }

  expected_coef_path <- file.path(project_dir, "expected", "ElasticNet_coefficients.tsv")
  expected_coef <- utils::read.delim(
    expected_coef_path, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = "NA"
  )
  for (row_index in seq_len(nrow(expected_coef))) {
    expected <- expected_coef[row_index, , drop = FALSE]
    actual <- elastic_coefficients[
      elastic_coefficients$mode == expected$mode &
        elastic_coefficients$feature_name == expected$feature_name,
      , drop = FALSE
    ]
    check_index <- check_index + 1L
    actual_value <- if (nrow(actual)) actual$coefficient[[1]] else NA_real_
    difference <- abs(actual_value - expected$coefficient[[1]])
    checks[[check_index]] <- data.frame(
      check = paste("Elastic Net coefficient", expected$mode, expected$feature_name),
      expected = expected$coefficient[[1]],
      actual = actual_value,
      absolute_difference = difference,
      tolerance = 1e-8,
      status = ifelse(is.finite(difference) && difference <= 1e-8, "PASS", "FAIL"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, checks)
}

message("Reading development data: ", normalizePath(input_file, mustWork = TRUE))
training <- read_training_data(input_file)
data <- training$data
label <- training$label
archived_fold_file <- file.path(
  project_dir, "expected", "fold_assignments_seed123.tsv"
)
if (seed == 123L && file.exists(archived_fold_file)) {
  archived_folds <- utils::read.delim(
    archived_fold_file,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  fold_match <- match(data$circRNA_id, archived_folds$circRNA_id)
  if (anyNA(fold_match) || anyDuplicated(archived_folds$circRNA_id)) {
    stop("Archived fold assignments do not match the training BSJ identifiers")
  }
  if (any(archived_folds$label[fold_match] != label)) {
    stop("Archived fold-assignment labels do not match the training data")
  }
  fold <- as.integer(archived_folds$fold[fold_match])
} else {
  fold <- make_stratified_folds(label, k = 5L, random_seed = seed)
}

fold_table <- data.frame(
  circRNA_id = data$circRNA_id,
  label = label,
  fold = fold,
  stringsAsFactors = FALSE
)
write_tsv(fold_table, file.path(output_dir, "fold_assignments.tsv"))

fold_counts <- do.call(rbind, lapply(seq_len(5), function(fold_number) {
  validation <- fold == fold_number
  data.frame(
    fold = fold_number,
    training_rows = sum(!validation),
    training_P1 = sum(label[!validation] == 1L),
    training_N1 = sum(label[!validation] == 0L),
    validation_rows = sum(validation),
    validation_P1 = sum(label[validation] == 1L),
    validation_N1 = sum(label[validation] == 0L)
  )
}))
write_tsv(fold_counts, file.path(output_dir, "fold_counts.tsv"))

if (refit_only) {
  message("Refit-only mode: using archived selected hyperparameters")
  selected <- selected_from_expected()
  xgboost_tuning <- NULL
  elastic_tuning <- NULL
} else {
  xgboost_grid <- expand.grid(
    max_depth = c(2, 3, 4, 5),
    learning_rate = c(0.01, 0.03, 0.05, 0.1),
    min_child_weight = c(1, 3, 5),
    subsample = c(0.7, 0.8, 0.9),
    colsample_bytree = c(0.7, 0.8, 1.0),
    gamma = c(0, 0.5),
    nrounds = c(100, 200, 400),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  elastic_grid <- expand.grid(
    alpha = seq(0, 1, by = 0.1),
    lambda = 10^seq(-4, 1, length.out = 50),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (quick_test) {
    message("Quick-test mode: validating one archived candidate per model")
    xgboost_grid <- unique(expected_best[expected_best$method == "XGBoost", c(
      "max_depth", "learning_rate", "min_child_weight", "subsample",
      "colsample_bytree", "gamma", "nrounds"
    )])
    elastic_grid <- unique(expected_best[expected_best$method == "Elastic Net", c(
      "alpha", "lambda"
    )])
  }

  message("Tuning ", nrow(xgboost_grid), " XGBoost candidates per dimension")
  xgboost_tuning <- tune_xgboost(
    data, label, fold, xgboost_grid, random_seed = seed, nthread = threads
  )
  write_tsv(xgboost_tuning, file.path(output_dir, "XGBoost_tuning.tsv"))

  message("Tuning ", nrow(elastic_grid), " Elastic Net candidates per dimension")
  elastic_tuning <- tune_elastic_net(data, label, fold, elastic_grid)
  write_tsv(elastic_tuning, file.path(output_dir, "ElasticNet_tuning.tsv"))
  selected <- selected_from_tuning(xgboost_tuning, elastic_tuning)
}

best_parameters <- make_parameter_summary(selected)
write_tsv(best_parameters, file.path(output_dir, "best_parameters.tsv"))

message("Generating selected-model out-of-fold predictions")
oof <- generate_selected_oof(
  data, label, fold, selected, random_seed = seed, nthread = threads
)
write_tsv(oof$predictions, file.path(output_dir, "out_of_fold_predictions.tsv"))
write_tsv(oof$performance, file.path(output_dir, "out_of_fold_performance.tsv"))

message("Refitting four selected models on all development BSJs")
final <- refit_all_models(
  data, label, selected, random_seed = seed, nthread = threads
)
write_tsv(
  final$xgboost_importance,
  file.path(output_dir, "XGBoost_feature_importance.tsv")
)
write_tsv(
  final$elastic_coefficients,
  file.path(output_dir, "ElasticNet_coefficients.tsv")
)
write_tsv(
  final$fitted_scores,
  file.path(output_dir, "training_fitted_scores.tsv")
)

manifest <- data.frame(
  item = c(
    "analysis_population", "P1_positive", "N1_negative", "seed", "CV_folds",
    "CV_stratified", "selection_primary", "selection_secondary",
    "XGBoost_candidates_per_dimension", "ElasticNet_candidates_per_dimension",
    "R_version", "xgboost_version", "glmnet_version"
  ),
  value = c(
    nrow(data), sum(label == 1L), sum(label == 0L), seed, 5, TRUE,
    "Mean fold AUROC (descending)", "Mean fold AUPRC (descending)",
    if (refit_only) "not rerun" else nrow(unique(xgboost_tuning[, c(
      "max_depth", "learning_rate", "min_child_weight", "subsample",
      "colsample_bytree", "gamma", "nrounds"
    )])),
    if (refit_only) "not rerun" else nrow(unique(elastic_tuning[, c("alpha", "lambda")])),
    R.version.string,
    as.character(utils::packageVersion("xgboost")),
    as.character(utils::packageVersion("glmnet"))
  ),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(output_dir, "model_manifest.tsv"))

checks <- compare_with_expected(best_parameters, final$elastic_coefficients)
write_tsv(checks, file.path(output_dir, "reproducibility_checks.tsv"))

session_file <- file(file.path(output_dir, "sessionInfo.txt"), open = "wt")
sink(session_file)
print(utils::sessionInfo())
sink()
close(session_file)

failed_checks <- sum(checks$status != "PASS")
message("Completed. Results saved under: ", normalizePath(output_dir, mustWork = TRUE))
message("Reproducibility checks: ", nrow(checks) - failed_checks, " PASS; ", failed_checks, " FAIL")
if (failed_checks > 0L) {
  warning(
    "One or more archived-result checks failed. Review reproducibility_checks.tsv. ",
    "This can reflect a quick-test run, a different model-selection result, or ",
    "software/platform differences."
  )
}
