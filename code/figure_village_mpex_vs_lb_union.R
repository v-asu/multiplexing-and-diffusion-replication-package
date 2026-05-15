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

## Plot village multiplexing against threshold-averaged local bridge share.

input_path <- "data/processed/rct_villages/rct_village_union_network_mechanisms.csv"
output_path <- "figures/fig_si_08_multiplexing_local_bridges.pdf"

plot_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

main <- function() {
  plot_df <- read.csv(input_path) %>%
    as_tibble() %>%
    filter(!is.na(M_i), !is.na(lb))

  plot_lb <- ggplot(plot_df, aes(x = lb, y = M_i)) +
    geom_point(size = 3, alpha = 0.7, color = "#333333") +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray90", linewidth = 0.8) +
    labs(
      x = "Local Bridge Share (Union Layer)",
      y = "Multiplexing Index",
      title = " "
    ) +
    plot_theme

  ggsave(output_path, plot_lb, width = 7, height = 5, dpi = 300)
  message("Saved village multiplexing vs local bridge share figure.")
}

main()
