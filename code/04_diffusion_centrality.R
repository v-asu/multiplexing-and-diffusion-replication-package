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

## This file calculates the Seed Set's Diffusion Centrality for the RCT villages
## Also give us input into regressions later

if (!exists("netlist_rct")) {
  netlist_rct = c("social", "kerorice", "advice", "decision")
}
netlist_rct = union(netlist_rct, "jati") ## adding jati to the list

##################### Data ####################################

load("data/processed/rct_villages/rct_network_adjacency_layers.RData") ## loading processed adjmats
jati = readRDS("data/processed/rct_villages/rct_hh_pair_jati_layer.rds")

backbone = readRDS("data/processed/rct_villages/rct_hh_pair_backbone_with_jati.rds") ## continuous backbone layer
backbone_no_jati = readRDS("data/processed/rct_villages/rct_hh_pair_backbone_no_jati.rds")


dat_hh_rct = haven::read_dta("data/raw/rct_villages/rct_hh_covariates.dta") ## household level data
dat_calls_rct = haven::read_dta("data/raw/rct_villages/rct_village_diffusion_outcomes.dta") ## village level data

#################### Lambda 1 ####################################

## Lambda_1 Calc
for (k in c(netlist_rct, "union_link", "intersect_link", "backbone", "union_wtd", "raw_wtd")) {
  x = eval(as.name(k))
  assign(paste0("lambda_1_", k), get_lambda_1(x, c(1:26,28:71)) %>% unlist)
}

lambda_1_rct = do.call(cbind, lapply(c(netlist_rct, "union_link", "intersect_link", "backbone", "union_wtd", "raw_wtd"), \(name) {
  get(paste0("lambda_1_", name))
})) %>%
  as.data.frame() %>%
  setNames(paste0("lambda_1_", c(netlist_rct, "union_link", "intersect_link", "backbone", "union_wtd", "raw_wtd"))) %>%
  cbind(village = c(1:26,28:71)) %>%
  filter(village != 40) ## Village 40 had issues with experiment implementation

write.csv(lambda_1_rct, "data/processed/rct_villages/lambda_1_rct.csv", row.names = F)

#################### Seed Set DC ####################################

## seed vertices
df_seed = dat_hh_rct %>% select(villageid, hhid, seed_dummy, vertex) %>%
  filter(seed_dummy == 1)

# Generate column names
column_names = c("village", paste0("sum", "_", c(netlist_rct, "union_link", "intersect_link", "backbone", "union_wtd","raw_wtd")))

DC_seed_list = lapply(c(1:26, 28:39, 41:71), \(k) {
  values = c(k, sapply(c(netlist_rct, "union_link", "intersect_link", "backbone", "union_wtd", "raw_wtd"), \(net) DC_out_seed(net, k, df_seed, lambda_1_rct)))
  setNames(values, column_names)
})

no_jati_layers = c("union_link_no_jati", "intersect_link_no_jati", "backbone_no_jati")
no_jati_column_names = c("link_no_jati", "intersect_link_no_jati", "backbone_no_jati")

lambda_1_rct_no_jati = do.call(cbind, lapply(no_jati_layers, \(name) {
  get_lambda_1(get(name), c(1:26, 28:71)) %>% unlist()
})) %>%
  as.data.frame() %>%
  setNames(paste0("lambda_1_", no_jati_layers)) %>%
  cbind(village = c(1:26, 28:71)) %>%
  filter(village != 40)

DC_seed_list_no_jati = lapply(c(1:26, 28:39, 41:71), \(k) {
  values = sapply(no_jati_layers, \(net) DC_out_seed(net, k, df_seed, lambda_1_rct_no_jati))
  setNames(c(k, values), c("village", no_jati_column_names))
})

df_DC_seed = bind_rows(DC_seed_list) %>%
  left_join(bind_rows(DC_seed_list_no_jati), by = "village")
df_DC_seed[df_DC_seed$village == 62, -1] = NA ## No seed data in this village

write.csv(df_DC_seed, "data/processed/rct_villages/rct_village_seed_diffusion_centrality.csv", row.names = FALSE)
