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

# Data ---------------------------------------------------------------------

adjmat_wave1_individual = readRDS("data/raw/mf_villages/gender/mf_individual_wave1.rds")

## covariate data (individual level) to get the gender variable (keys were already merged in here)
df_covs = haven::read_dta("data/raw/mf_villages/mf_ind_covariates_wave1.dta")
df_covs = df_covs %>% select(village, adjmatrix_key, pid, resp_gend)

main_layers = c("social", "kerorice", "money", "advice", "decision", "medic", "templecompany")

## This gives us a list of matrices for each network
for (net in main_layers) {
    assign(net, adjmat_wave1_individual[[net]])
}

# Analysis ----------------------------------------------------------------

union_link = list()

for (i in 1:75) {
union_link[[i]] = (social[[i]] + kerorice[[i]] + money[[i]] + advice[[i]] + decision[[i]] + templecompany[[i]] + medic[[i]] > 0)*1
}

## Now we calculate the multiplexing stats for each individual

total_links = vector(mode = "list", length = 75)
total_people = vector(mode = "list", length = 75)
M_index = vector(mode = "list", length = 75)


L = length(main_layers)

for (i in 1:75) {

  total_links[[i]] = Reduce(`+`, map(main_layers, ~rowSums(get(.x)[[i]]))) / L

    total_people[[i]] = rowSums(union_link[[i]])

  M_index[[i]] = total_links[[i]]/total_people[[i]]

}

## Now we want to map this back to the gender data from individual covariates
df_M = map(M_index, \(x) data.frame(M_i = x, adjmatrix_key = 1:length(x))) %>%
    map2_df(c(1:12, 14:21, 23:77), ~mutate(.x, village = .y))

## convert rownames to a column called adjmatrix_key
rownames(df_M) = NULL

df_M = df_M %>%
    left_join(df_covs, by = c("village", "adjmatrix_key"))

## keeping only those observations for which we have covariate data, left with 16984 obs
df_M = df_M %>% filter(!is.na(resp_gend))

df_M = df_M %>%
    mutate(gender = ifelse(resp_gend == 1, "Male", "Female"))

# Plotting ----------------------------------------------------------------


df_M_collapsed = df_M %>%
    group_by(village, gender) %>%
    summarise(M_i = mean(M_i, na.rm = TRUE)) %>%
    ungroup()

gender_mf_wave1_plot = ggplot(df_M_collapsed) +
    geom_density(aes(x = M_i, fill = gender), alpha = 0.8) +
    labs(x = "Multiplexing Index", y = "Density", fill = "Gender") +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("figures/fig_si_05_gender_multiplexing_mf_wave1.pdf", gender_mf_wave1_plot, width = 7, height = 4)
