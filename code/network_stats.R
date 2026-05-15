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

## This file calculates the network statistics for both RCT and Microfinance villages

net_list_rct = c("social", "kerorice", "advice", "decision", "union_link", "intersect_link", "jati") ## 196 HH on avg


netlist_mf = c("social", "kerorice", "advice", "decision", "money", "templecompany", "medic", "union_link", "intersect_link", "jati")

## RCT villages

load("data/processed/rct_villages/rct_network_adjacency_layers.RData")
jati = read_rds("data/processed/rct_villages/rct_hh_pair_jati_layer.rds")

stats_out = map(net_list_rct, network_stats, dims = c(1:26, 28:71))

network_stats_diff("social", "kerorice", stats_out, net_list_rct)
network_stats_diff("social", "advice", stats_out, net_list_rct)
network_stats_diff("social", "decision", stats_out, net_list_rct)
network_stats_diff("advice", "decision", stats_out, net_list_rct)


df_stats_rct = map(net_list_rct, network_stats_avg, dims = c(1:26, 28:71)) %>% 
  bind_rows() %>% round(3) %>% 
  mutate(Network = c(net_list_rct)) %>% 
  relocate(Network)

# df_stats_rct %>% datasummary_df(output = "tables/tab_aux_descriptive_statistics_rct.tex",
#                             title = "Summary Statistics: Diffusion RCT")


## MF Villages

load("data/processed/mf_villages/mf_network_adjacency_layers.RData")
jati = read_rds("data/processed/mf_villages/mf_hh_pair_jati_layer.rds")

stats_out = map(netlist_mf, network_stats, dims = c(1:12, 14:21, 23:77))

network_stats_diff("social", "kerorice", stats_out, netlist_mf)
network_stats_diff("social", "advice", stats_out, netlist_mf)
network_stats_diff("social", "decision", stats_out, netlist_mf)
network_stats_diff("advice", "decision", stats_out, netlist_mf)

df_stats_mf = map(netlist_mf, network_stats_avg, dims = c(1:12, 14:21, 23:77)) %>%
  bind_rows() %>% round(3) %>% 
  mutate(Network = c(netlist_mf)) %>% 
  relocate(Network)

# df_stats_mf %>% datasummary_df(output = "tables/tab_aux_descriptive_statistics_mf.tex",
#                             title = "Summary Statistics: Microfinance Villages")


bind_rows(df_stats_mf, df_stats_rct) %>% 
  kbl(format = "latex",  booktabs = T) %>% 
  pack_rows("Microfinance villages", 1, nrow(df_stats_mf), hline_before = T, hline_after = T) %>%
  pack_rows("RCT villages", nrow(df_stats_mf) + 1, nrow(df_stats_mf) + nrow(df_stats_rct), hline_before = T, hline_after = T) %>% 
  save_kable(file = "tables/tab_si_01_descriptive_statistics.tex")
