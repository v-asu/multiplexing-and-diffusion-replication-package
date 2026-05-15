# -------------------------------------------------------------------------
# Setup: load package dependencies and shared functions.
# -------------------------------------------------------------------------

if (file.exists("code/00_packages.R")) {
  source("code/00_packages.R")
} else if (file.exists("replication_package/code/00_packages.R")) {
  source("replication_package/code/00_packages.R")
} else {
  stop("Could not locate replication_package bootstrap.")
}

source("code/00_functions.R")

output_dirs <- c(
  "data/processed/mf_villages",
  "data/processed/rct_villages",
  "figures",
  "tables"
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

cleanup_env <- function() {
  keep <- c(
    "cleanup_env",
    "run_in_current_session",
    "run_isolated",
    "run_scripts",
    "run_python_script",
    "script"
  )
  rm(list = setdiff(ls(envir = .GlobalEnv), keep), envir = .GlobalEnv)
  invisible(gc())
}

run_in_current_session <- function(scripts) {
  for (script in scripts) {
    source(script)
  }
}

run_isolated <- function(script) {
  cleanup_env()
  source(script)
}

run_scripts <- function(scripts) {
  for (script in scripts) {
    run_isolated(script)
  }
}

run_python_script <- function(script) {
  python <- Sys.getenv("REPLICATION_PACKAGE_PYTHON", unset = "")
  mpl_config_dir <- Sys.getenv("MPLCONFIGDIR", unset = "")

  if (!nzchar(python)) {
    python <- Sys.which("python3")
  }

  if (!nzchar(python)) {
    python <- Sys.which("python")
  }

  if (!nzchar(python)) {
    stop("Could not locate a Python interpreter for ", script, ".")
  }

  if (!nzchar(mpl_config_dir)) {
    mpl_config_dir <- file.path(tempdir(), "matplotlib")
    dir.create(mpl_config_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cache_home <- Sys.getenv("XDG_CACHE_HOME", unset = "")

  if (!nzchar(cache_home)) {
    cache_home <- file.path(tempdir(), "cache")
    dir.create(cache_home, recursive = TRUE, showWarnings = FALSE)
  }

  status <- system2(
    python,
    script,
    env = c(
      paste0("MPLCONFIGDIR=", mpl_config_dir),
      paste0("XDG_CACHE_HOME=", cache_home)
    )
  )

  if (!identical(status, 0L)) {
    stop("Python script failed: ", script)
  }
}

# -------------------------------------------------------------------------
# Data preparation: construct processed inputs used by later scripts.
# -------------------------------------------------------------------------

run_in_current_session(c(
  "code/01_data_prep.R"
))

# -------------------------------------------------------------------------
# Main-paper and SI outputs built directly from the prepared pair-level data.
# -------------------------------------------------------------------------

run_in_current_session(c(
  "code/02_pca.R"
))

# -------------------------------------------------------------------------
# Supporting-information outputs built directly from the prepared pair-level data.
# -------------------------------------------------------------------------

run_in_current_session(c(
  "code/layer_correlation.R"
))

# -------------------------------------------------------------------------
# Data preparation: construct diffusion and backbone inputs for regressions.
# -------------------------------------------------------------------------

run_in_current_session(c(
  "code/03_backbone_with_jati.R",
  "code/03b_backbone_no_jati.R",
  "code/04_diffusion_centrality.R"
))

# -------------------------------------------------------------------------
# Main-paper outputs: figures and tables used in the main text.
# Some of these scripts also write SI outputs.
# -------------------------------------------------------------------------

run_scripts(c(
  "code/05_main_regs.R",
  "code/06_mpex_vs_diffusion.R",
  "code/07_gender_MF_wave2.R",
  "code/08_mpex_vs_degree.R"
))

run_python_script("code/09_plot_sims.py")

# -------------------------------------------------------------------------
# Supporting-information outputs: additional figures and tables used in the SI.
# -------------------------------------------------------------------------

run_scripts(c(
  "code/06b_mpex_vs_diffusion_continuous.R",
  "code/06c_mpex_vs_diffusion_all_layers.R",
  "code/hh_diffusion_centrality.R",
  "code/complex_path_length_build.R",
  "code/figure_village_mpex_vs_lb_union.R",
  "code/figure_village_mpex_vs_complex_path_union.R",
  "code/table_individual_mplex_on_union_measures.R",
  "code/network_layer_stats.R",
  "code/network_stats.R",
  "code/gender_MF_wave1.R",
  "code/gender_stratified_mpex.R"
))
