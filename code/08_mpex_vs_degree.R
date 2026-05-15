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

## This file looks at how multiplexing index varies across degree in the union 

################################# Data #######################################

## Network graphs
load("data/processed/rct_villages/rct_network_adjacency_layers.RData")

if (!exists("netlist_rct")) {
  netlist_rct = c("social", "kerorice", "advice", "decision", "jati")
}


############################### Processing ####################################

total_links = vector(mode = "list", length = 71)
total_people = vector(mode = "list", length = 71)
M_index = vector(mode = "list", length = 71)

main_layers = netlist_rct[!netlist_rct %in% c("jati", "distance_w")]
L = length(main_layers)


for (i in c(1:26, 28:71)) {

  total_links[[i]] = Reduce(`+`, map(main_layers, ~rowSums(get(.x)[[i]]))) / L

  total_people[[i]] = map_dbl(1:nrow(get(main_layers[1])[[i]]), \(x) {
    sum_matrix = Reduce(`+`, map(main_layers, ~get(.x)[[i]]))
    sum(sum_matrix[x, ] > 0)
  })

  M_index[[i]] = total_links[[i]]/total_people[[i]]

}

M_index = compact(M_index) ## remove the empty list cells
total_people = compact(total_people) ## remove the empty list cells

df = data.frame(degree = unlist(total_people), m_i = unlist(M_index),
                village = rep(c(1:26, 28:71), times = unlist(map(total_people, length))))

df_collapsed = df %>%
    group_by(degree) %>%
    summarize(m_i = mean(m_i, na.rm=T))

mpex_vs_degree_plot = ggplot(df_collapsed, aes(x = degree, y = m_i)) +
    geom_point() +
    lims(y = c(0,1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "Degree in union layer", y = "Multiplexing Index",
    title = " ")

ggsave("figures/fig_main_03c_multiplexing_by_degree.pdf", mpex_vs_degree_plot, width = 6, height = 4)
