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

## Backbone layer with collapsed layers + Jati layer

netlist_rct = union(netlist_rct, "jati")

#################### Data ####################################
df_edges_rct = readRDS("data/processed/rct_villages/rct_hh_pair_edges.rds")

##################### Eigenvalues ####################################

## Keeping the relevant columns
X_rct = df_edges_rct %>%
  select(all_of(netlist_rct))

X_rct_act = as.matrix(X_rct)
mean_X_rct = apply(X_rct_act, 2, mean)
X_rct_act = apply(X_rct_act, 2, scale, scale = TRUE)

## calculating the covariance matrix and its eigenvalues/vectors
X_rct_cov = cov(X_rct_act)
X_rct_eig = eigen(X_rct_cov)

if (all(X_rct_eig$vectors[,1] < 0)) {
  X_rct_loadings = X_rct_eig$vectors * -1
} else{
  X_rct_loadings = X_rct_eig$vectors
}


rownames(X_rct_loadings) = netlist_rct
colnames(X_rct_loadings) = paste0("PC", 1:ncol(X_rct_loadings))

laddle_out = ICtest::PCAladle(as.matrix(X_rct_act))
k = laddle_out$k

####################### Backbone Calculation ###########################

pca_scores_rct = as.matrix(X_rct) %*% X_rct_loadings

## create the score automatically such that it uses the first laddle_out$k PCs
weights = X_rct_eig$values[1:k] / sum(X_rct_eig$values[1:k])
S_rct = pca_scores_rct[, 1:k] %*% weights

df_S_rct = cbind(df_edges_rct %>% select(village, id1, id2), S_rct)

backbone_rct_cont = list()
for (i in c(1:26, 28:71)) {

  B_edge = df_S_rct %>% filter(village == i) %>%
    select(id1, id2, S_rct)

  backbone_rct_cont[[i]] = graph_from_data_frame(B_edge, directed = FALSE) %>%
    as_adjacency_matrix(sparse = FALSE, attr = "S_rct")

  diag(backbone_rct_cont[[i]]) = 0
}

saveRDS(backbone_rct_cont, "data/processed/rct_villages/rct_hh_pair_backbone_with_jati.rds")
