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

## Build canonical union-analysis datasets for the retained outputs:
## 1. household regressions of multiplexing on clustering / PLC / LB
## 2. village multiplexing vs complex path length
## 3. village multiplexing vs local bridge share

village_ids <- c(1:26, 28:39, 41:71)
thresholds <- 2:6
output_dir <- "data/processed/rct_villages"

sanitize_adjacency <- function(adj) {
  mat <- as.matrix(adj)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- ((mat + t(mat)) > 0) * 1L
  diag(mat) <- 0L
  storage.mode(mat) <- "integer"
  mat
}

build_union_adjacency <- function(social_adj, kerorice_adj, advice_adj, decision_adj) {
  sanitize_adjacency((social_adj + kerorice_adj + advice_adj + decision_adj) > 0)
}

largest_component_nodes <- function(adj) {
  graph <- graph_from_adjacency_matrix(adj, mode = "max", diag = FALSE)
  membership <- igraph::components(graph)$membership
  component_sizes <- tabulate(membership)
  which(membership == which.max(component_sizes))
}

closed_neighborhoods <- function(adj) {
  lapply(
    seq_len(nrow(adj)),
    function(node) sort.int(c(node, which(adj[node, ] > 0L)))
  )
}

closed_neighborhood_matrix <- function(neighborhoods, n_nodes) {
  mat <- matrix(FALSE, nrow = n_nodes, ncol = n_nodes)

  for (node in seq_len(n_nodes)) {
    mat[node, neighborhoods[[node]]] <- TRUE
  }

  mat
}

bridge_width_matrix <- function(adj, neighborhood_matrix) {
  n_nodes <- nrow(adj)
  widths <- matrix(0L, nrow = n_nodes, ncol = n_nodes)

  for (source in seq_len(n_nodes)) {
    source_neighborhood <- neighborhood_matrix[source, ]
    tied_to_source <- colSums(adj[source_neighborhood, , drop = FALSE]) > 0L

    for (target in seq_len(n_nodes)) {
      if (source == target) {
        next
      }

      target_neighborhood <- neighborhood_matrix[target, ]
      overlap_count <- sum(source_neighborhood & target_neighborhood)
      disjoint_target <- target_neighborhood & !source_neighborhood
      reinforcement_count <- sum(disjoint_target & tied_to_source)

      widths[source, target] <- as.integer(overlap_count + reinforcement_count)
    }
  }

  widths
}

local_bridge_share <- function(width_row, threshold_value, source) {
  bridge_targets <- which(width_row >= 1L)
  bridge_targets <- bridge_targets[bridge_targets != source]

  if (length(bridge_targets) == 0L) {
    return(0)
  }

  mean(width_row[bridge_targets] >= threshold_value)
}

complex_diffusion <- function(adj, initial_active, threshold_value) {
  active <- rep(FALSE, nrow(adj))
  active[initial_active] <- TRUE

  repeat {
    exposures <- rowSums(adj[, active, drop = FALSE])
    new_active <- (!active) & (exposures >= threshold_value)

    if (!any(new_active)) {
      break
    }

    active[new_active] <- TRUE
  }

  active
}

bfs_edge_distances <- function(adj, allowed_nodes, source) {
  n_nodes <- nrow(adj)
  allowed <- rep(FALSE, n_nodes)
  allowed[allowed_nodes] <- TRUE

  distances <- rep.int(Inf, n_nodes)
  visited <- rep(FALSE, n_nodes)
  queue <- integer(n_nodes)
  head <- 1L
  tail <- 1L

  distances[source] <- 0
  visited[source] <- TRUE
  queue[tail] <- source

  while (head <= tail) {
    current <- queue[head]
    head <- head + 1L

    neighbors <- which(adj[current, ] > 0L & allowed & !visited)
    if (length(neighbors) == 0L) {
      next
    }

    for (neighbor in neighbors) {
      visited[neighbor] <- TRUE
      distances[neighbor] <- distances[current] + 1L
      tail <- tail + 1L
      queue[tail] <- neighbor
    }
  }

  distances
}

seed_complex_path_metrics <- function(adj, neighborhoods, seed, threshold_value) {
  n_nodes <- nrow(adj)
  seed_neighborhood <- neighborhoods[[seed]]
  eligible_targets <- setdiff(seq_len(n_nodes), seed_neighborhood)

  if (length(eligible_targets) == 0L) {
    return(list(plc_i = 0, reachable_targets = 0L))
  }

  active <- complex_diffusion(adj, seed_neighborhood, threshold_value)
  active_nodes <- which(active)
  distances <- bfs_edge_distances(adj, active_nodes, seed)
  plc_ij <- rep.int(0, length(eligible_targets))

  reachable <- active[eligible_targets] & is.finite(distances[eligible_targets])
  plc_ij[reachable] <- distances[eligible_targets][reachable] + 1L

  list(
    plc_i = mean(plc_ij),
    reachable_targets = sum(plc_ij > 0L)
  )
}

compute_multiplexing_i <- function(social_adj, kerorice_adj, advice_adj, decision_adj) {
  combined_adj <- social_adj + kerorice_adj + advice_adj + decision_adj
  total_links <- rowSums(combined_adj) / 4
  total_people <- rowSums(combined_adj > 0)
  multiplexing_i <- total_links / total_people
  multiplexing_i[total_people == 0] <- NA_real_
  multiplexing_i
}

compute_village_outputs <- function(village_id, network_data) {
  social_adj <- sanitize_adjacency(network_data$social[[village_id]])
  kerorice_adj <- sanitize_adjacency(network_data$kerorice[[village_id]])
  advice_adj <- sanitize_adjacency(network_data$advice[[village_id]])
  decision_adj <- sanitize_adjacency(network_data$decision[[village_id]])
  union_adj <- build_union_adjacency(
    social_adj = social_adj,
    kerorice_adj = kerorice_adj,
    advice_adj = advice_adj,
    decision_adj = decision_adj
  )

  multiplexing_full <- compute_multiplexing_i(
    social_adj = social_adj,
    kerorice_adj = kerorice_adj,
    advice_adj = advice_adj,
    decision_adj = decision_adj
  )

  graph_union <- graph_from_adjacency_matrix(union_adj, mode = "max", diag = FALSE)
  clustering_full <- transitivity(graph_union, type = "localundirected", isolates = "zero")

  keep_nodes <- largest_component_nodes(union_adj)
  union_gc <- union_adj[keep_nodes, keep_nodes, drop = FALSE]
  neighborhoods <- closed_neighborhoods(union_gc)
  neighborhood_matrix <- closed_neighborhood_matrix(neighborhoods, nrow(union_gc))
  widths <- bridge_width_matrix(union_gc, neighborhood_matrix)

  threshold_results <- lapply(
    thresholds,
    function(threshold_value) {
      seed_metrics <- lapply(
        seq_len(nrow(union_gc)),
        function(seed) seed_complex_path_metrics(union_gc, neighborhoods, seed, threshold_value)
      )

      tibble(
        threshold = threshold_value,
        seed_node = seq_len(nrow(union_gc)),
        lb_i = vapply(
          seq_len(nrow(union_gc)),
          function(seed) local_bridge_share(widths[seed, ], threshold_value, seed),
          numeric(1)
        ),
        plc_i = vapply(seed_metrics, `[[`, numeric(1), "plc_i"),
        reachable_targets = vapply(seed_metrics, `[[`, integer(1), "reachable_targets")
      )
    }
  ) %>%
    bind_rows()

  household_df <- threshold_results %>%
    group_by(seed_node) %>%
    summarise(
      plc_i = mean(plc_i, na.rm = TRUE),
      lb_i = mean(lb_i, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      village = village_id,
      original_node_id = keep_nodes,
      multiplexing_i = multiplexing_full[keep_nodes],
      clustering_i = clustering_full[keep_nodes]
    ) %>%
    select(village, seed_node, original_node_id, multiplexing_i, plc_i, lb_i, clustering_i)

  village_df <- threshold_results %>%
    summarise(
      plc = mean(plc_i, na.rm = TRUE),
      lb = mean(lb_i, na.rm = TRUE)
    ) %>%
    mutate(
      village = village_id,
      M_i = mean(multiplexing_full, na.rm = TRUE)
    ) %>%
    select(village, M_i, plc, lb)

  list(household = household_df, village = village_df)
}

main <- function() {

  load("data/processed/rct_villages/rct_network_adjacency_layers.RData")
  network_data <- list(
    social = social,
    kerorice = kerorice,
    advice = advice,
    decision = decision,
    union_link = union_link
  )

  results <- lapply(village_ids, compute_village_outputs, network_data = network_data)

  household_df <- bind_rows(map(results, "household")) %>%
    arrange(village, seed_node)
  village_df <- bind_rows(map(results, "village")) %>%
    arrange(village)

  stopifnot(
    nrow(village_df) == length(village_ids),
    all(!is.na(village_df$M_i)),
    all(!is.na(village_df$plc)),
    all(!is.na(village_df$lb)),
    all(!is.na(household_df$multiplexing_i)),
    all(!is.na(household_df$plc_i)),
    all(!is.na(household_df$lb_i)),
    all(!is.na(household_df$clustering_i))
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    village_df,
    file.path(output_dir, "rct_village_union_network_mechanisms.csv"),
    row.names = FALSE
  )

  write.csv(
    household_df,
    file.path(output_dir, "rct_hh_union_network_mechanisms.csv"),
    row.names = FALSE
  )

}

main()
