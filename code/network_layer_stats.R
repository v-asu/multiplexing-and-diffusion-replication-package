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

# Network Layer Comparison for RCT Data
# Creates clean LaTeX table with network statistics

# Load data
load("data/processed/rct_villages/rct_network_adjacency_layers.RData")
jati = readRDS("data/processed/rct_villages/rct_hh_pair_jati_layer.rds")

# Define villages to analyze
village_ids = c(1:26, 28:71)

# Define layers to compare
layers = list(
    advice = advice,
    kerorice = kerorice,
    social = social,
    decision = decision,
    jati = jati,
    union = union_link,
    intersect = intersect_link
)

# Function to compute network statistics
compute_stats <- function(adj_mat) {
    g = graph_from_adjacency_matrix(adj_mat, mode = "max", diag = FALSE)

    deg = degree(g)

    # Components
    comp = igraph::components(g)

    # Number of unique triangles
    adj_mat = as.matrix(as_adjacency_matrix(g))
    n_triangles = sum(diag(adj_mat %*% adj_mat %*% adj_mat)) / 6

    # Only compute path-related metrics if graph is connected (or has giant component)
    avg_path_length = NA
    diameter = NA
    if (comp$no == 1) {
        # Connected graph - compute path metrics
        avg_path_length = mean_distance(g, directed = FALSE)
        diameter = diameter(g, directed = FALSE)
    } else if (comp$no > 1) {
        # Use largest component for path metrics
        g_giant = induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
        avg_path_length = mean_distance(g_giant, directed = FALSE)
        diameter = diameter(g_giant, directed = FALSE)
    }

    # Centralization measures (updated for igraph 2.0.0)
    # Degree centralization: how star-like is the network
    deg_cent = centr_degree(g)$centralization

    # Betweenness centralization: how concentrated is betweenness
    betweenness_cent = centr_betw(g)$centralization

    # Closeness centralization - use giant component if disconnected
    if (comp$no == 1) {
        closeness_cent = centr_clo(g)$centralization
    } else {
        g_giant = induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
        closeness_cent = centr_clo(g_giant)$centralization
    }

    # Bridge detection: count of bridges (edges whose removal disconnects graph)
    bridges = bridges(g)
    n_bridges = length(bridges)
    bridge_prop = n_bridges / ecount(g)  # proportion of edges that are bridges

    # Edge betweenness statistics (mean and max) 
    edge_btw = edge_betweenness(g, directed = FALSE)
    mean_edge_btw = mean(edge_btw)
    max_edge_btw = max(edge_btw)

    data.frame(
        density = edge_density(g, loops = FALSE),
        mean_degree = mean(deg),
        sd_degree = sd(deg),
        n_triangles = n_triangles,
        clustering = transitivity(g, type = "global"),
        n_components = comp$no,
        giant_comp_prop = max(comp$csize) / vcount(g),
        avg_path_length = avg_path_length,
        diameter = diameter,
        degree_centralization = deg_cent,
        betweenness_centralization = betweenness_cent,
        closeness_centralization = closeness_cent,
        n_bridges = n_bridges,
        bridge_proportion = bridge_prop,
        mean_edge_betweenness = mean_edge_btw,
        max_edge_betweenness = max_edge_btw
    )
}

# Compute stats for all layers and villages
results = data.frame()

for (layer_name in names(layers)) {
    layer_data = layers[[layer_name]]

    for (village_id in village_ids) {
        adj_mat = layer_data[[village_id]]

        if (is.null(adj_mat) || !is.matrix(adj_mat) || nrow(adj_mat) == 0) {
            next
        }

        stats = compute_stats(adj_mat)
        stats$layer = layer_name
        stats$village = village_id
        results = bind_rows(results, stats)
    }
}

# Summary statistics by layer
summary_stats = results %>%
    group_by(layer) %>%
    summarise(
        Density = mean(density),
        Mean_degree = mean(mean_degree),
        SD_degree = mean(sd_degree),
        Triangles = mean(n_triangles),
        Clustering = mean(clustering, na.rm = TRUE),
        Components = mean(n_components),
        Giant_comp = mean(giant_comp_prop),
        Avg_path_length = mean(avg_path_length, na.rm = TRUE),
        Diameter = mean(diameter, na.rm = TRUE),
        Degree_centralization = mean(degree_centralization, na.rm = TRUE),
        Betweenness_centralization = mean(betweenness_centralization, na.rm = TRUE),
        Closeness_centralization = mean(closeness_centralization, na.rm = TRUE),
        N_bridges = mean(n_bridges),
        Bridge_proportion = mean(bridge_proportion, na.rm = TRUE),
        Mean_edge_betweenness = mean(mean_edge_betweenness, na.rm = TRUE),
        Max_edge_betweenness = mean(max_edge_betweenness, na.rm = TRUE),
        .groups = "drop"
    )

# Format for table with clear names (using makecell for Overleaf compatibility)
summary_stats = summary_stats %>%
    mutate(
        Density = sprintf("%.4f", Density),
        Mean_degree = sprintf("%.2f", Mean_degree),
        SD_degree = sprintf("%.2f", SD_degree),
        Triangles = sprintf("%.0f", Triangles),
        Clustering = sprintf("%.3f", Clustering),
        Components = sprintf("%.1f", Components),
        Giant_comp = sprintf("%.3f", Giant_comp),
        Avg_path_length = sprintf("%.2f", Avg_path_length),
        Diameter = sprintf("%.1f", Diameter),
        Degree_centralization = sprintf("%.4f", Degree_centralization),
        Betweenness_centralization = sprintf("%.4f", Betweenness_centralization),
        Closeness_centralization = sprintf("%.4f", Closeness_centralization),
        N_bridges = sprintf("%.1f", N_bridges),
        Mean_edge_betweenness = sprintf("%.1f", Mean_edge_betweenness),
        Max_edge_betweenness = sprintf("%.1f", Max_edge_betweenness)
    ) %>%
    select(layer, Density, Mean_degree, SD_degree,
           Triangles, Clustering, Components, Giant_comp,
           Avg_path_length, Diameter, Degree_centralization,
           Betweenness_centralization, Closeness_centralization,
           N_bridges, Bridge_proportion,
            Mean_edge_betweenness, Max_edge_betweenness) %>%
    rename(
        Layer = layer,
        density = Density,
        degree = Mean_degree,
        SD_degree = SD_degree,
        triangles = Triangles,
        clustering = Clustering,
        components = Components,
        GC_prop = Giant_comp,
        avg_path_len = Avg_path_length,
        diameter = Diameter,
        deg_centr = Degree_centralization,
        betw_centr = Betweenness_centralization,
        close_centr = Closeness_centralization,
        n_bridges = N_bridges,
        mean_eb = Mean_edge_betweenness,
        max_eb = Max_edge_betweenness
    )

# Convert to numeric for LaTeX
table_df_out = summary_stats %>%
    mutate(across(-Layer, as.numeric))

# Split into two tables
table1_cols = c("Layer", "density", "degree", "SD_degree",
                "triangles", "clustering", "components", "GC_prop")
table1 = table_df_out %>% select(all_of(table1_cols))

table2_cols = c("Layer", "avg_path_len", "diameter", "deg_centr", "betw_centr",
                "close_centr", "n_bridges", "Bridge_proportion", "mean_eb", "max_eb")
table2 = table_df_out %>% select(all_of(table2_cols))

table1 %>% datasummary_df(output = "tables/tab_si_12_network_layer_basic_metrics.tex",
                            title = "Network Statistics by Layer (Basic Metrics)")

table2 %>% datasummary_df(output = "tables/tab_si_13_network_layer_path_metrics.tex",
                            title = "Network Statistics by Layer (Path and Centralization Metrics)")
