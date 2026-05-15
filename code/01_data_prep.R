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

## Loads the raw network data and processes it.

########################## Micro-finance Village Raw Data ############################

dat_graph_mf = readMat("data/raw/mf_villages/mf_network_adjacency_layers.mat") ## Adj Mats (already symmetric)
## 75 villages, seq - 1:12, 14:21, 23:77

## pairwise household distances derived from Wave 2 GPS coordinates.
## The raw GPS coordinates are excluded from the public replication package.
dat_distance_mf = readRDS("data/raw/mf_villages/mf_hh_pair_distance_winsorized.rds")

## w2_vertex to hhid match.
dat_vertex_mf = readRDS("data/raw/mf_villages/mf_hh_vertex_crosswalk_wave2.rds")

## loading jati data
dat_jati_mf = haven::read_dta("data/raw/mf_villages/mf_hh_covariates.dta")
dat_jati_mf = dat_jati_mf %>%
  select(village, newHHID, hhcaste2_0, hhsubc2_1) %>%
  left_join(dat_vertex_mf, by = c("village", "newHHID" = "newhhid")) %>%
  mutate(hhsubc2_1 = ifelse(hhsubc2_1 %in% c(-999, -333), -966, hhsubc2_1),
         hhsubc2_1 = ifelse(is.na(hhsubc2_1), -966, hhsubc2_1),
         hhcaste2_0 = ifelse(hhcaste2_0 %in% c(-999, -888, -333), -966, hhcaste2_0))

########################## RCT Village Raw Data ############################
dat_graph_rct = readMat("data/raw/rct_villages/rct_network_adjacency_layers.mat") ## reading all adjmats
## 70 villages, seq - 1:26, 28:71

dat_keys_rct = readMat("data/raw/rct_villages/rct_hh_vertex_crosswalk.mat") ## keys for adjmat HHs
dat_keys_rct = dat_keys_rct[[1]]

## loading the caste co-variates
raw_covariates = list.files(path = "data/raw/rct_villages/hh_covariates_by_village",
                            pattern = "covariates",
                            full.names = TRUE)
raw_covariates = gtools::mixedsort(raw_covariates)
dat_covariates_rct = map(raw_covariates, fread)
dat_covariates_rct = rbindlist(dat_covariates_rct) ## 13756 obs

## For key to covariate vars see HH census instrument
## putting all keys in one seq
df_keys_rct = map(dat_keys_rct, as.data.frame)
df_keys_rct = map2_df(df_keys_rct, c(1:71), ~mutate(.x, village = .y))
df_keys_rct = df_keys_rct %>% set_names("v_id", "hh_id", "village")

dat_covariates_rct = dat_covariates_rct %>%
  mutate(caste_category = ifelse(caste_category %in% c(-999, -666, -333), -966, caste_category),
         subcaste = ifelse(subcaste %in% c(-888, -666, -333), -966, subcaste))

############################# Micro-finance Processing ################################

# Note - All microfinance village networks are already symmetric

## creating empty containers
visitgo = list()
visitcome = list()
rel = list()
nonrel = list()
medic = list()
keroricego = list()
keroricecome = list()
templecompany = list()
bormoney = list()
lendmoney = list()
decision = list()
advice = list()
jati = list()

## collapsing some of the network layers
social = list() ## visitgo + visitcome + rel + nonrel
kerorice = list() ## keroricecome + keroricego
money = list() ## bormoney + lendmoney
union_link_no_jati = list()
intersect_link_no_jati = list()
union_link = list()
intersect_link = list()
union_wtd = list() ## summing up the collapsed layers
raw_wtd = list() ## summing up the raw layers (directed)

for (i in c(1:12, 14:21, 23:77)) {

  visitgo[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[1]][[1]])
  visitcome[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[2]][[1]])
  nonrel[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[3]][[1]])
  rel[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[4]][[1]])
  medic[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[5]][[1]])
  keroricego[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[6]][[1]])
  keroricecome[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[7]][[1]])
  templecompany[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[8]][[1]])
  bormoney[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[9]][[1]])
  lendmoney[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[10]][[1]])
  decision[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[11]][[1]])
  advice[[i]] = as.matrix(dat_graph_mf[[1]][[i]][[1]][[13]][[1]])

  ## collapsing some of the network layers
  social[[i]] = (visitgo[[i]] + visitcome[[i]] + rel[[i]] + nonrel[[i]] > 0)*1
  kerorice[[i]] = (keroricego[[i]] + keroricecome[[i]] > 0)*1
  money[[i]] = (bormoney[[i]] + lendmoney[[i]] > 0)*1

  ## create a raw weighted network (just sum of all layers, can be directed too)
  raw_wtd[[i]] = (visitgo[[i]] + visitcome[[i]] + rel[[i]] + nonrel[[i]] +
                  medic[[i]] + keroricego[[i]] + keroricecome[[i]] +
                  templecompany[[i]] + bormoney[[i]] + lendmoney[[i]] +
                  decision[[i]] + advice[[i]])

  union_link_no_jati[[i]] = (social[[i]] + medic[[i]] +
                    kerorice[[i]] + templecompany[[i]] + money[[i]] + advice[[i]] + decision[[i]] > 0)*1

  intersect_link_no_jati[[i]] = (social[[i]] * medic[[i]] *
                    kerorice[[i]] * templecompany[[i]] *
                    money[[i]] * advice[[i]] * decision[[i]] > 0)*1

  union_link[[i]] = union_link_no_jati[[i]]
  intersect_link[[i]] = intersect_link_no_jati[[i]]

  union_wtd[[i]] = (social[[i]] + medic[[i]] +
                    kerorice[[i]] + templecompany[[i]] + money[[i]] + advice[[i]] + decision[[i]])

}

## saving all the processed adjmat list for future use
save(visitgo, visitcome, nonrel, rel, medic, keroricego, keroricecome,
     templecompany, bormoney, lendmoney, advice, social, decision,
     kerorice, money, union_link_no_jati, intersect_link_no_jati,
     union_link, intersect_link, union_wtd, raw_wtd,
     file = "data/processed/mf_villages/mf_network_adjacency_layers.RData")

netlist = c("visitgo", "visitcome", "nonrel", "rel", "medic",
            "keroricego", "keroricecome", "templecompany",
            "bormoney", "lendmoney", "decision", "advice",
            "social", "kerorice", "money")

## Reshaping each AdjMat to edgelist dataframe
for (i in netlist) {
  assign(i, mat_to_df(i, dims = c(1:12, 14:21, 23:77)))
}

## binding all networks together
all_links = bind_cols(bormoney, advice, decision, keroricecome,
                      keroricego, lendmoney, medic, nonrel,
                      rel, templecompany, visitcome, visitgo,
                      social, kerorice, money)

all_links = all_links %>%
  select(contains(netlist), village...4, id1...1, id2...2) %>%
  rename(id1 = id1...1, id2 = id2...2, village = village...4) %>%
  relocate(id1, id2, village) %>%
  mutate(across(everything(), as.integer)) ##1932224 unique edges

rm(list = netlist)
gc()

## Switching to data table for SPEED
setDT(all_links)
setDT(dat_distance_mf)
setDT(dat_vertex_mf)
setDT(dat_jati_mf)

## Merging in pairwise distance used by the analysis.
all_links = merge(all_links, dat_distance_mf, all.x = TRUE,
                  by = c("village", "id1", "id2"))

## Merging Caste and Jati Info
all_links = merge(all_links, dat_jati_mf, all.x = TRUE, by.x = c("village","id1"),
                  by.y = c("village","w2_vertex"))
setnames(all_links, old = c("hhcaste2_0", "hhsubc2_1", "newHHID"), new = c("caste1", "jati1", "newhhid1"))
all_links = merge(all_links, dat_jati_mf, all.x = TRUE, by.x = c("village","id2"),
                  by.y = c("village","w2_vertex"))
setnames(all_links, old = c("hhcaste2_0", "hhsubc2_1", "newHHID"), new = c("caste2", "jati2", "newhhid2"))

all_links[, caste := fifelse(caste1 == caste2, 1, 0)]
all_links[, jati := fifelse(jati1 == jati2, 1, 0)]

## give anyone with -966 caste/jati a 0 link
all_links[, caste := fifelse(caste1 == -966 | caste2 == -966, 0, caste)]
all_links[, jati := fifelse(jati1 == -966 | jati2 == -966, 0, jati)]

all_links = all_links[, c("id1", "id2", "village", netlist, "caste", "jati", "distance_w"),
                      with = FALSE]

setDF(all_links)

saveRDS(all_links, "data/processed/mf_villages/mf_hh_pair_edges.rds")

## save jati as a network layer
for (i in c(1:12, 14:21, 23:77)){

  ## Generating Jati Adjmat
  jati_edge = all_links %>% filter(village == i) %>%
    select(id1, id2, jati)

  jati[[i]] = graph_from_data_frame(jati_edge, directed = FALSE) %>%
    as_adjacency_matrix(sparse = FALSE, attr = "jati")

  diag(jati[[i]]) = 0

}

saveRDS(jati, "data/processed/mf_villages/mf_hh_pair_jati_layer.rds")

## load the saved adjmats again
load("data/processed/mf_villages/mf_network_adjacency_layers.RData")

for (i in c(1:12, 14:21, 23:77)) {

  ## add jati layer to the default aggregate layers
  union_link[[i]] = (union_link[[i]] + jati[[i]] > 0)*1

  intersect_link[[i]] = (intersect_link[[i]] * jati[[i]] > 0)*1

  union_wtd[[i]] = union_wtd[[i]] + jati[[i]]

}

save(visitgo, visitcome, nonrel, rel, medic, keroricego, keroricecome,
   templecompany, bormoney, lendmoney, advice, social, decision,
   kerorice, money, union_link_no_jati, intersect_link_no_jati,
   union_link, intersect_link, union_wtd, raw_wtd,
   file = "data/processed/mf_villages/mf_network_adjacency_layers.RData")

rm(visitgo, visitcome, nonrel, rel, medic, keroricego, keroricecome,
   templecompany, bormoney, lendmoney, advice, social, decision,
   kerorice, money, union_link_no_jati, intersect_link_no_jati,
   union_link, intersect_link, union_wtd, raw_wtd)

################################ RCT Village Processing ################################

visitgo = list()
visitcome = list()
keroricego = list()
keroricecome = list()
decision = list()
advice = list()

union_link = list()
intersect_link = list()
union_wtd = list()
raw_wtd = list()
union_link_no_jati = list()
intersect_link_no_jati = list()

jati = list()
social = list()
kerorice = list()
## no money here since we dont record bor/lend relations

for (i in c(1:26, 28:71)) {

  keroricego[[i]] = dat_graph_rct[[1]][[i]][[1]][[1]][[1]]
  keroricego[[i]] = (keroricego[[i]] + t(keroricego[[i]]))>0
  diag(keroricego[[i]]) = 0

  keroricecome[[i]] = dat_graph_rct[[1]][[i]][[1]][[2]][[1]]
  keroricecome[[i]] = (keroricecome[[i]] + t(keroricecome[[i]]))>0
  diag(keroricecome[[i]]) = 0

  visitcome[[i]] = dat_graph_rct[[1]][[i]][[1]][[3]][[1]]
  visitcome[[i]] = (visitcome[[i]] + t(visitcome[[i]]))>0
  diag(visitcome[[i]]) = 0

  visitgo[[i]] = dat_graph_rct[[1]][[i]][[1]][[4]][[1]]
  visitgo[[i]] = (visitgo[[i]] + t(visitgo[[i]]))>0
  diag(visitgo[[i]]) = 0

  decision[[i]] = dat_graph_rct[[1]][[i]][[1]][[5]][[1]]
  decision[[i]] = (decision[[i]] + t(decision[[i]]))>0
  diag(decision[[i]]) = 0

  advice[[i]] = dat_graph_rct[[1]][[i]][[1]][[6]][[1]]
  advice[[i]] = (advice[[i]] + t(advice[[i]]))>0
  diag(advice[[i]]) = 0

  ## Social
  social[[i]] = (visitcome[[i]] + visitgo[[i]] > 0)*1

  ## Favors
  kerorice[[i]] =  (keroricego[[i]] + keroricecome[[i]] > 0)*1

  ## Raw Weighted Network
  raw_wtd[[i]] = (dat_graph_rct[[1]][[i]][[1]][[1]][[1]] +
                  dat_graph_rct[[1]][[i]][[1]][[2]][[1]] +
                  dat_graph_rct[[1]][[i]][[1]][[3]][[1]] +
                  dat_graph_rct[[1]][[i]][[1]][[4]][[1]] +
                  dat_graph_rct[[1]][[i]][[1]][[5]][[1]] +
                  dat_graph_rct[[1]][[i]][[1]][[6]][[1]])

  diag(raw_wtd[[i]]) = 0

  ## Union Link
  union_link_no_jati[[i]] = (kerorice[[i]] + social[[i]] + decision[[i]] + advice[[i]] > 0)*1

  ## intersection
  intersect_link_no_jati[[i]] = (kerorice[[i]] * social[[i]] * decision[[i]] * advice[[i]] > 0)*1

  union_link[[i]] = union_link_no_jati[[i]]
  intersect_link[[i]] = intersect_link_no_jati[[i]]

  ## union weighted
  union_wtd[[i]] = (kerorice[[i]] + social[[i]] + decision[[i]] + advice[[i]])
}

## saving all the processed adjmat list for future use
save(keroricecome, keroricego, visitcome, visitgo, advice, decision,
     union_link_no_jati, intersect_link_no_jati,
     union_link, intersect_link, social, kerorice, union_wtd, raw_wtd,
     file = "data/processed/rct_villages/rct_network_adjacency_layers.RData")

netlist = c("keroricego", "keroricecome","visitcome", "visitgo",
            "decision", "advice", "social",
            "kerorice")

## reshaping each adj matrix to pair level dataframe
for (i in netlist) {
  assign(i, mat_to_df(i, dims = c(1:26, 28:71)))
}

## Appending everything together
all_links = bind_cols(advice, decision, keroricecome,
                 keroricego, visitcome, visitgo,
                 social, kerorice)

all_links = all_links %>%
  select(contains(netlist), village...4, id1...1, id2...2) %>%
  rename(id1 = id1...1, id2 = id2...2, village = village...4) %>%
  relocate(id1, id2) %>%
  mutate(across(everything(), as.integer)) ## 1450245 unique edges

gc()

setDT(all_links)
setDT(df_keys_rct)
setDT(dat_covariates_rct)

## binding the keys to hhid matching here
all_links = merge(all_links, df_keys_rct, all.x = TRUE, by.x = c("village", "id1"),
                  by.y = c("village", "v_id"))
setnames(all_links, old = "hh_id", new = "hh_id1")
all_links = merge(all_links, df_keys_rct, all.x = TRUE, by.x = c("village", "id2"),
                  by.y = c("village", "v_id"))
setnames(all_links, old = "hh_id", new = "hh_id2")

## merging the caste data now
dat_covariates_rct = dat_covariates_rct[, .(caste_category, subcaste, villageid, hh_id)]

all_links = merge(all_links, dat_covariates_rct, all.x = TRUE, by.x = c("village", "hh_id1"),
                  by.y = c("villageid", "hh_id"))
setnames(all_links, old = c("caste_category", "subcaste"), new = c("caste1", "jati1"))
all_links = merge(all_links, dat_covariates_rct, all.x = TRUE, by.x = c("village", "hh_id2"),
                  by.y = c("villageid", "hh_id"))
setnames(all_links, old = c("caste_category", "subcaste"), new = c("caste2", "jati2"))

all_links[, caste := fifelse(caste1 == caste2, 1, 0)]
all_links[, jati := fifelse(jati1 == jati2, 1, 0)]

## give anyone with -966 caste/jati a 0 link
all_links[, caste := fifelse(caste1 == -966 | caste2 == -966, 0, caste)]
all_links[, jati := fifelse(jati1 == -966 | jati2 == -966, 0, jati)]

all_links = all_links[, c("id1", "id2", "village", netlist, "caste", "jati"),
                      with = FALSE]

setDF(all_links)

## constructing jati level adjmats
for (i in c(1:26, 28:71)){

  ## Generating Jati Adjmat
  jati_edge = all_links %>% filter(village == i) %>%
    select(id1, id2, jati)

  jati[[i]] = graph_from_data_frame(jati_edge, directed = FALSE) %>%
    as_adjacency_matrix(sparse = FALSE, attr = "jati")

  diag(jati[[i]]) = 0

}

##saving jati layer
saveRDS(jati, "data/processed/rct_villages/rct_hh_pair_jati_layer.rds")

## saving edgelist
saveRDS(all_links, "data/processed/rct_villages/rct_hh_pair_edges.rds")

## load the saved adjmats again
load("data/processed/rct_villages/rct_network_adjacency_layers.RData")

for (i in c(1:26, 28:71)) {

  ## add jati layer to the default aggregate layers
  union_link[[i]] = (union_link[[i]] + jati[[i]] > 0)*1

  intersect_link[[i]] = (intersect_link[[i]] * jati[[i]] > 0)*1

  union_wtd[[i]] = union_wtd[[i]] + jati[[i]]
}

save(keroricecome, keroricego, visitcome, visitgo, advice, decision,
   union_link_no_jati, intersect_link_no_jati,
   union_link, intersect_link, social, kerorice, union_wtd, raw_wtd,
   file = "data/processed/rct_villages/rct_network_adjacency_layers.RData")

rm(keroricecome, keroricego, visitcome, visitgo, advice, decision,
   union_link_no_jati, intersect_link_no_jati,
   union_link, intersect_link, social, kerorice, union_wtd, raw_wtd)
