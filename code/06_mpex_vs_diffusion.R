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

## Generates table for multiplexing vs diffusion in the paper

##################### Data ####################################

df_DC_seed = read.csv("data/processed/rct_villages/rct_village_seed_diffusion_centrality.csv") ## seed set DC
dat_calls_rct = haven::read_dta("data/raw/rct_villages/rct_village_diffusion_outcomes.dta") ## village level data
df_edges_rct = readRDS("data/processed/rct_villages/rct_hh_pair_edges.rds") ##pair level links data
load("data/processed/rct_villages/rct_network_adjacency_layers.RData") ## All layers
jati = read_rds("data/processed/rct_villages/rct_hh_pair_jati_layer.rds") ## jati layer
dat_hh_rct = haven::read_dta("data/raw/rct_villages/rct_hh_covariates.dta") ## household level data

#################### Processing ###########################

df_seed = dat_hh_rct %>% select(villageid, hhid, seed_dummy, vertex) %>% 
  filter(seed_dummy == 1)

# we want to create m_i = number of unique links (p) /total links (d)
total_links = vector(mode = "list", length = 71)
total_people = vector(mode = "list", length = 71)
M_index = vector(mode = "list", length = 71)
degree_vill = vector(mode = "list", length = 71)

L = 4 ## social, kerorice, advice, decision

for (i in c(1:26, 28:39, 41:71)) {

total_links[[i]] = (rowSums(social[[i]]) + rowSums(kerorice[[i]]) +
                     + rowSums(advice[[i]]) + rowSums(decision[[i]]))/L

total_people[[i]] = map(1:nrow(union_link[[i]]), \(x) length(which((social[[i]] + kerorice[[i]] + decision[[i]] + advice[[i]])[x, ] > 0))) %>% unlist()

M_index[[i]] = total_links[[i]]/total_people[[i]]

degree_vill[[i]] = mean(rowSums(advice[[i]]))
}

m_i = map(M_index, mean, na.rm = T) %>% discard(is.na) %>% unlist()
avg_degree = map(total_links, mean, na.rm = T) %>% discard(is.na) %>% unlist() ## total degree across layers (scaled by number of layers)
degree_vill = unlist(degree_vill) ## info layer degree

## setting up the data for regression
dat_calls_rct = dat_calls_rct %>% 
  filter(villageid %in% c(1:26, 28:39, 41:71)) %>% 
  mutate(calls_per_HH = CallsReceived/num_hh_random,
         dummy_3_seeds = ifelse(num_seeds == 3, 1, 0)) %>% 
  select(CallsReceived, dummy_3_seeds, calls_per_HH, num_hh_random, villageid) %>% 
  mutate(num_hh_random_sq = (num_hh_random)^2,
         num_hh_random_cb = (num_hh_random)^3)

reg_data_rct = dat_calls_rct %>% 
  left_join(df_DC_seed, by = c("villageid" = "village"))

reg_data_rct = reg_data_rct %>% arrange(villageid)

reg_data_rct$M_i = m_i
reg_data_rct$avg_deg = avg_degree

reg_data_rct = reg_data_rct %>% 
  mutate(sum_advice = as.vector(scale(sum_advice)),
         M_50 = ifelse(M_i > quantile(M_i, 0.5), 1, 0))

reg_list = list()
reg_list[["(1)"]] = lm_robust(calls_per_HH ~ M_50*sum_advice + dummy_3_seeds + avg_deg,  reg_data_rct)

modelsummary(reg_list, coef_omit = '(Intercept)|dummy_3_seeds|num_hh_random|avg_deg',
             gof_omit = '[^Num.Obs]',
             coef_rename = c("M_50:sum_advice" = "High Multiplexing X Seed Set Centrality",
                             "sum_advice" = "Seed Set Centrality",
                          "M_50" = "High Multiplexing"
                          ),
             statistic = c("std.error", "[{p.value}]"),
             output = "latex") %>% 
  add_header_above(c(" " = 1, "Calls per Household" = 1)) %>% 
  save_kable("tables/tab_main_02_multiplexing_diffusion.tex")
