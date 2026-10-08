options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

required <- c(xgboost = "3.2.1.1", glmnet = "4.1-10")
for (package_name in names(required)) {
  target_version <- required[[package_name]]
  installed_version <- if (requireNamespace(package_name, quietly = TRUE)) {
    as.character(utils::packageVersion(package_name))
  } else {
    NA_character_
  }
  if (is.na(installed_version) || installed_version != target_version) {
    remotes::install_version(
      package_name,
      version = target_version,
      repos = getOption("repos"),
      upgrade = "never"
    )
  }
}

cat("R:", R.version.string, "\n")
for (package_name in names(required)) {
  cat(package_name, as.character(utils::packageVersion(package_name)), "\n")
}

