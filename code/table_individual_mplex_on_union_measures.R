if (file.exists("code/00_packages.R")) {
  source("code/00_packages.R")
} else if (file.exists("00_packages.R")) {
  source("00_packages.R")
} else if (file.exists("replication_package/code/00_packages.R")) {
  source("replication_package/code/00_packages.R")
} else {
  stop("Could not locate replication_package bootstrap.")
}
source("code/00_functions.R")

## Estimate the retained household-level regressions of multiplexing on
## clustering, complex path length, and local bridge share.

input_path <- "data/processed/rct_villages/rct_hh_union_network_mechanisms.csv"
table_path <- "tables/tab_si_11_multiplexing_union_mechanisms.tex"

fit_models <- function(df) {
  list(
    "(1)" = lm_robust(
      multiplexing_i ~ clustering_i + factor(village),
      data = df,
      clusters = village
    ),
    "(2)" = lm_robust(
      multiplexing_i ~ plc_i + factor(village),
      data = df,
      clusters = village
    ),
    "(3)" = lm_robust(
      multiplexing_i ~ lb_i + factor(village),
      data = df,
      clusters = village
    )
  )
}

save_table <- function(model_list) {
  options(modelsummary_factory_latex = "kableExtra")

  modelsummary(
    model_list,
    coef_omit = "(Intercept)|factor\\(village\\)",
    coef_rename = c(
      "clustering_i" = "Clustering",
      "plc_i" = "Complex Path Length",
      "lb_i" = "Local Bridge Share"
    ),
    gof_omit = "AIC|BIC|Log.Lik|F|RMSE|Std.Errors|Adj",
    statistic = c("std.error", "[{p.value}]"),
    title = "Multiplexing and network structure",
    output = "latex"
  ) %>%
    add_header_above(c(" " = 1, "Individual Multiplexing" = length(model_list))) %>%
    footnote(general = "Multiplexing index is calculated at the household level using the union layer. We include village fixed effects and standard errors clustered at the village level.") %>%
    save_kable(table_path)
}

main <- function() {
  df <- read.csv(input_path) %>%
    as_tibble() %>%
    filter(
      !is.na(multiplexing_i),
      !is.na(clustering_i),
      !is.na(plc_i),
      !is.na(lb_i)
    )

  models <- fit_models(df)

  save_table(models)

  message("Saved household-level union-measure regression table.")
}

main()
