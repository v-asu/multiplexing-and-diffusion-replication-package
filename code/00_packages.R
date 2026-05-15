if (file.exists("code/00_packages.R")) {
  package_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
} else if (file.exists("00_packages.R") && basename(getwd()) == "code") {
  package_root <- normalizePath("..", winslash = "/", mustWork = TRUE)
} else if (file.exists("replication_package/code/00_packages.R")) {
  package_root <- normalizePath("replication_package", winslash = "/", mustWork = TRUE)
} else {
  stop("Could not locate replication_package root.")
}

if (normalizePath(getwd(), winslash = "/", mustWork = TRUE) != package_root) {
  setwd(package_root)
}

if (!isTRUE(getOption("replication_package.packages_loaded"))) {
  packages <- c(
    "tidyverse", "haven", "patchwork", "knitr", "R.matlab", "igraph",
    "Matrix", "expm", "kableExtra", "estimatr", "modelsummary",
    "glmnet", "ggrepel", "data.table", "ggdist",
    "ggsci", "factoextra", "gtools", "ICtest", "fixest",
    "marginaleffects", "viridis", "xtable", "reshape2", "broom"
  )

  missing_packages <- packages[!vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]

  if (length(missing_packages) > 0) {
    message(
      "Installing missing packages from CRAN: ",
      paste(missing_packages, collapse = ", ")
    )

    install.packages(
      missing_packages,
      repos = "https://cloud.r-project.org"
    )

    still_missing <- missing_packages[!vapply(
      missing_packages,
      requireNamespace,
      logical(1)
    )]

    if (length(still_missing) > 0) {
      stop(
        "Failed to install required packages: ",
        paste(still_missing, collapse = ", ")
      )
    }
  }

  invisible(lapply(
    packages,
    function(pkg) suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  ))

  suppressPackageStartupMessages(library(igraph))

  options(
    modelsummary_factory_latex = "kableExtra",
    ggrepel.max.overlaps = Inf,
    replication_package.packages_loaded = TRUE
  )
}
