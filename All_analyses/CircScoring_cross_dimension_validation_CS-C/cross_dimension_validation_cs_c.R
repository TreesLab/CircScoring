#!/usr/bin/env Rscript

# CircScoring: cross-dimension validation of CS-C.
#
# Separate univariable logistic regressions quantify the association of each
# finalized CS-C score with orthogonal circRNA-reliability evidence. No model
# refitting, resampling, multiple-testing correction, or plotting is performed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    "Usage: Rscript cross_dimension_validation_cs_c.R ",
    "<analysis_input.tsv[.gz]> <output_directory> [--relaxed-qc]"
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
relaxed_qc <- length(args) == 3L && identical(args[[3]], "--relaxed-qc")
if (length(args) == 3L && !relaxed_qc) stop("Unknown option: ", args[[3]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(digits = 17)

# The original circIG export may contain a first row of column groups followed
# by the actual column-name row. Detect that layout without modifying the data.
first_lines <- readLines(input_file, n = 2L, warn = FALSE)
if (length(first_lines) < 1L) stop("Input file is empty.")
delimiter <- if (grepl("\\t", first_lines[[1]], fixed = FALSE)) "\t" else ","
split_line <- function(x) strsplit(x, delimiter, fixed = TRUE)[[1]]
header_1 <- split_line(first_lines[[1]])
header_2 <- if (length(first_lines) >= 2L) split_line(first_lines[[2]]) else character()
score_names <- c("XGBoost_CS-C_score", "ElasticNet_CS-C_score")
skip_rows <- if (all(score_names %in% header_1)) {
  0L
} else if (all(score_names %in% header_2)) {
  1L
} else {
  stop("Could not identify the header row containing the two CS-C score columns.")
}

d <- read.table(
  input_file,
  header = TRUE,
  sep = delimiter,
  skip = skip_rows,
  quote = "\"",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)

resolve_column <- function(data_names, alternatives, description) {
  hit <- alternatives[alternatives %in% data_names]
  if (length(hit) == 0L) {
    stop(
      "Missing ", description, " column. Accepted names: ",
      paste(alternatives, collapse = ", ")
    )
  }
  hit[[1]]
}

full_length_col <- resolve_column(
  names(d),
  c("FL-circAS or circFL_seq", "FL-circAS_or_circFL_seq"),
  "full-length support"
)
mouse_col <- resolve_column(
  names(d),
  c("mouse conserved", "mouse_conserved", "mouse conservation", "mouse_conservation"),
  "mouse conservation"
)
db_col <- resolve_column(
  names(d),
  c("#7_DB", "#7_DBs", "7_DB", "7_DBs"),
  "seven-database recurrence"
)

required_scores <- c("XGBoost_CS-C_score", "ElasticNet_CS-C_score")
missing_scores <- setdiff(required_scores, names(d))
if (length(missing_scores) > 0L) {
  stop("Missing score columns: ", paste(missing_scores, collapse = ", "))
}

numeric_blank_zero <- function(x, label) {
  x[is.na(x) | trimws(as.character(x)) == ""] <- 0
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value))) stop("Non-numeric value encountered in ", label, ".")
  value
}

score_data <- lapply(required_scores, function(column) {
  value <- suppressWarnings(as.numeric(d[[column]]))
  if (any(!is.finite(value))) stop("Missing or non-finite values in ", column, ".")
  value
})
names(score_data) <- c("XGBoost", "Elastic Net")

full_length <- numeric_blank_zero(d[[full_length_col]], full_length_col)
mouse <- numeric_blank_zero(d[[mouse_col]], mouse_col)
db_recurrence <- numeric_blank_zero(d[[db_col]], db_col)

outcomes <- list(
  "Full-length support" = as.integer(full_length >= 1),
  "Mouse conservation" = as.integer(mouse >= 1)
)
for (threshold in 1:7) {
  outcomes[[paste0("DB recurrence >=", threshold)]] <-
    as.integer(db_recurrence >= threshold)
}

# These controls identify the fixed analysis-eligible circIG population used in
# the retained analysis. They prevent silent use of a different atlas release.
expected_n <- 4455197L
expected_events <- c(
  "Full-length support" = 800687L,
  "Mouse conservation" = 94840L,
  "DB recurrence >=1" = 1291639L,
  "DB recurrence >=2" = 348111L,
  "DB recurrence >=3" = 214436L,
  "DB recurrence >=4" = 146985L,
  "DB recurrence >=5" = 101435L,
  "DB recurrence >=6" = 45429L,
  "DB recurrence >=7" = 10287L
)
expected_score_mean <- c(
  "XGBoost" = 0.17149202949464576,
  "Elastic Net" = 0.14720568763471145
)
expected_score_sd <- c(
  "XGBoost" = 0.2810835459295679,
  "Elastic Net" = 0.2885071754762009
)

observed_events <- vapply(outcomes, sum, numeric(1))
if (!relaxed_qc) {
  if (nrow(d) != expected_n) {
    stop("Unexpected N: observed ", nrow(d), ", expected ", expected_n, ".")
  }
  if (!identical(as.numeric(observed_events), as.numeric(expected_events[names(outcomes)]))) {
    stop(
      "Evidence-event counts do not match the fixed analysis population.\nObserved: ",
      paste(names(observed_events), observed_events, sep = "=", collapse = "; "),
      "\nExpected: ",
      paste(names(expected_events), expected_events, sep = "=", collapse = "; ")
    )
  }
  observed_score_mean <- vapply(score_data, mean, numeric(1))
  observed_score_sd <- vapply(score_data, sd, numeric(1))
  if (any(abs(observed_score_mean - expected_score_mean) > 5e-12) ||
      any(abs(observed_score_sd - expected_score_sd) > 5e-12)) {
    stop(
      "CS-C score moments do not match the fixed analysis population. ",
      "Check the atlas release and score columns."
    )
  }
}

fit_one_logistic_model <- function(y, raw_score, model_name, predictor_name) {
  score_mean <- mean(raw_score)
  score_sd <- sd(raw_score) # sample SD (denominator N - 1)
  score_z <- (raw_score - score_mean) / score_sd

  design <- cbind(`(Intercept)` = 1, score_z = score_z)
  fit <- glm.fit(
    x = design,
    y = y,
    family = binomial(link = "logit"),
    control = glm.control(epsilon = 1e-8, maxit = 100L, trace = FALSE)
  )
  fit$call <- match.call()
  class(fit) <- c("glm", "lm")
  coefficient_table <- coef(summary(fit))
  beta <- unname(coefficient_table["score_z", "Estimate"])
  standard_error <- unname(coefficient_table["score_z", "Std. Error"])
  z_value <- beta / standard_error
  p_value <- 2 * pnorm(-abs(z_value))
  log10_p <- (log(2) + pnorm(-abs(z_value), log.p = TRUE)) / log(10)
  critical_value <- qnorm(0.975)

  data.frame(
    Model = model_name,
    Predictor = predictor_name,
    N = length(y),
    N_event = sum(y),
    N_nonevent = length(y) - sum(y),
    Event_rate = mean(y),
    Predictor_mean = score_mean,
    Predictor_SD = score_sd,
    Log_OR = beta,
    SE_Log_OR = standard_error,
    Z = z_value,
    P_value = p_value,
    log10_P = log10_p,
    OR_per_1SD = exp(beta),
    CI95_low = exp(beta - critical_value * standard_error),
    CI95_high = exp(beta + critical_value * standard_error),
    Converged = isTRUE(fit$converged),
    Iterations = fit$iter,
    stringsAsFactors = FALSE
  )
}

result_rows <- list()
result_index <- 1L
for (outcome_name in names(outcomes)) {
  y <- outcomes[[outcome_name]]
  if (outcome_name == "Full-length support") {
    evidence <- "Full-length support"
    threshold <- "present"
  } else if (outcome_name == "Mouse conservation") {
    evidence <- "Mouse conservation"
    threshold <- "present"
  } else {
    evidence <- "DB recurrence"
    threshold <- sub("DB recurrence >=", "", outcome_name, fixed = TRUE)
  }

  for (model_name in names(score_data)) {
    predictor_name <- if (model_name == "XGBoost") {
      "XGBoost_CS-C_score"
    } else {
      "ElasticNet_CS-C_score"
    }
    fitted <- fit_one_logistic_model(
      y = y,
      raw_score = score_data[[model_name]],
      model_name = model_name,
      predictor_name = predictor_name
    )
    fitted <- cbind(
      Outcome = gsub(">=", "≥", outcome_name, fixed = TRUE),
      Evidence = evidence,
      Threshold = threshold,
      fitted,
      stringsAsFactors = FALSE
    )
    result_rows[[result_index]] <- fitted
    result_index <- result_index + 1L
  }
}

results <- do.call(rbind, result_rows)
event_counts <- data.frame(
  Outcome = gsub(">=", "≥", names(outcomes), fixed = TRUE),
  Evidence = c("Full-length support", "Mouse conservation", rep("DB recurrence", 7)),
  Threshold = c("present", "present", as.character(1:7)),
  N_total = nrow(d),
  N_event = as.numeric(observed_events),
  N_nonevent = nrow(d) - as.numeric(observed_events),
  Event_rate = as.numeric(observed_events) / nrow(d),
  stringsAsFactors = FALSE
)

write.table(
  results,
  file.path(output_dir, "Cross_dimension_validation_CS-C_logistic_regression.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)
write.table(
  event_counts,
  file.path(output_dir, "Cross_dimension_validation_CS-C_event_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

message(
  "Analysis completed. Results written to: ",
  normalizePath(output_dir, mustWork = TRUE)
)
