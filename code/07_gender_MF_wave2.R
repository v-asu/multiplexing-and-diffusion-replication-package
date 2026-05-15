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

## Using wave 2 individual data to look at how multiplexing varies with gender.


# Data ---------------------------------------------------------------------
netlist = c(
  "giveadvice", "helpdecision", "keroricecome", "keroricego",
  "lendmoney", "borrowmoney", "medic", "rel", "nonrel",
  "templecompany", "visitcome", "visitgo"
)

# Load processed adjacency matrices
adjmat_wave2_individual = readRDS("data/raw/mf_villages/gender/mf_individual_wave2.rds")

# Covariates (for gender)
df_covs = read_dta("data/raw/mf_villages/mf_ind_covariates_wave2.dta")

# Load networks ------------------------------------------------------------
for (net in netlist) {
    assign(net, adjmat_wave2_individual[[net]])
}

# Processing ----------------------------------------------------------------

## This gives us a list of matrices for each network
for (net in netlist) {
    assign(net, adjmat_wave2_individual[[net]])
}

## create collapsed networks (to be consistent with rest of the paper)

social = list()
kerorice = list()
money = list()
advice = list()
decision = list()
union_link = list()

for (i in 1:75) {

social[[i]] = (visitgo[[i]] + visitcome[[i]] + rel[[i]] + nonrel[[i]] > 0)*1
kerorice[[i]] = (keroricego[[i]] + keroricecome[[i]] > 0)*1
money[[i]] = (borrowmoney[[i]] + lendmoney[[i]] > 0)*1
advice[[i]] = giveadvice[[i]]
decision[[i]] = helpdecision[[i]]

union_link[[i]] = (social[[i]] + kerorice[[i]] + money[[i]] + advice[[i]] + decision[[i]] + templecompany[[i]] + medic[[i]] > 0)*1

}

# M index ------------------------------------------------------------------


total_links = vector(mode = "list", length = 75)
total_people = vector(mode = "list", length = 75)
M_index = vector(mode = "list", length = 75)

main_layers = c("social", "kerorice", "money", "advice", "decision", "medic", "templecompany")

L = length(main_layers)

for (i in 1:75) {

  total_links[[i]] = Reduce(`+`, map(main_layers, ~rowSums(get(.x)[[i]]))) / L

  total_people[[i]] = rowSums(union_link[[i]])

  M_index[[i]] = total_links[[i]]/total_people[[i]]

}

## Now we want to map this back to the gender data from individual covariates
df_M = map(M_index, \(x) data.frame(M_i = x, adjmatrix_key = 1:length(x), UniqueID = as.numeric(names(x)))) %>%
    map2_df(c(1:12, 14:21, 23:77), ~mutate(.x, village = .y))

# degree in the union network
df_degree = map(total_people, \(x) data.frame(degree = x, adjmatrix_key = 1:length(x), UniqueID = as.numeric(names(x)))) %>%
    map2_df(c(1:12, 14:21, 23:77), ~mutate(.x, village = .y)) 

df_degree = df_degree %>% select(village, UniqueID, degree)

df_M = df_M %>%
    left_join(df_degree, by = c("village", "UniqueID"))

## Village is encoded in newHHID
df_covs = df_covs %>%
    mutate(village = floor(newHHID / 1000))


## now we can merge in our M index data (note resp_gender = 1 male)
df_M = df_M %>%
    left_join(df_covs, by = c("village", "UniqueID"))

df_M = df_M %>%
    filter(!is.na(respgender0_3))

df_M = df_M %>%
mutate(gender = ifelse(respgender0_3 == 1, "Male", "Female"))

# Plots -----------------------------------------------------------------

df_M_collapsed = df_M %>%
    group_by(village, gender) %>%
    summarise(M_i = mean(M_i, na.rm = TRUE)) %>%
    ungroup()

## collapsed plot
p1 = ggplot(df_M_collapsed) +
    geom_density(aes(x = M_i, fill = gender), alpha = 0.8) +
    labs(x = "Multiplexing Index", y = "Density", fill = " ", title = " ") +
    theme_bw()+
    xlim(0, 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = c(0.9,0.7),
     axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

ggsave("figures/fig_si_04a_gender_multiplexing_density_aggregate.pdf", p1, width = 7, height = 5, units = "in", dpi = 300)

ks_out = ks.test(df_M_collapsed %>% filter(gender == "Male") %>% pull(M_i),
        df_M_collapsed %>% filter(gender == "Female") %>% pull(M_i),
        alternative = "greater")

## plot ecdf
p2 = ggplot(df_M_collapsed) +
    stat_ecdf(aes(x = M_i, color = gender), show.legend = F) +
    labs(x = "Multiplexing Index", y = "Cumulative Density",
            color = "Gender", title = " ") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("figures/fig_si_04b_gender_multiplexing_cdf_aggregate.pdf", p2, width = 7, height = 5, units = "in", dpi = 300)

## Raw plot

p3 = ggplot(df_M) +
    geom_density(aes(x = M_i, fill = gender), alpha = 0.8) +
    labs(x = "Multiplexing Index", y = "Density", fill = " ", title = " ") +
    theme_bw()+
    xlim(0, 1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = c(0.9,0.7),
     axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

ggsave("figures/fig_main_03d_multiplexing_by_gender.pdf", p3, width = 7, height = 5, units = "in", dpi = 300)
ggsave("figures/fig_si_04c_gender_multiplexing_density_individual.pdf", p3, width = 7, height = 5, units = "in", dpi = 300)

ks_out = ks.test(df_M %>% filter(gender == "Male") %>% pull(M_i),
        df_M %>% filter(gender == "Female") %>% pull(M_i),
        alternative = "greater")

## plot ecdf
p4 = ggplot(df_M) +
    stat_ecdf(aes(x = M_i, color = gender), show.legend = F) +
    labs(x = "Multiplexing Index", y = "Cumulative Density",
            color = "Gender", title = " ") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("figures/fig_si_04d_gender_multiplexing_cdf_individual.pdf", p4, width = 7, height = 5, units = "in", dpi = 300)
