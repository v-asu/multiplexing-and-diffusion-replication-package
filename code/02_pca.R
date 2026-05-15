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

## This file runs PCA on both Microfinance Village and RCT Village Dataset

## network list to be used (this continues to other scripts post this as well)
netlist_mf = c("social", "kerorice", "money", "advice", "decision", "medic", "templecompany", "jati", "distance_w")
netlist_rct = c("social", "kerorice", "advice", "decision", "jati")

################### Microfinance Data ##############################
df_edges_mf = readRDS("data/processed/mf_villages/mf_hh_pair_edges.rds") ## 1932224 obs

################### RCT Data ####################################
df_edges_rct = readRDS("data/processed/rct_villages/rct_hh_pair_edges.rds") ## 1450245 obs

################## PCA Microfinance Villages #######################

## All Layers

X_mf = df_edges_mf %>%
  select(all_of(netlist_mf)) %>%
  drop_na(jati, distance_w) %>%
  rename("temple" = "templecompany",
         "distance" = "distance_w")

all_pca = prcomp(X_mf, scale. = TRUE, center = TRUE)
all_pca$rotation = all_pca$rotation * -1

## Next we output the plots for each of the above PCA
options(ggrepel.max.overlaps = Inf)

## PC plot Microfinance
loading_plot_mf = factoextra::fviz_pca_var(all_pca, col.var = "red",
                                           axes = c(1, 2),
                                           labelsize = 3,
                                           arrowsize = 0.4,
                                           repel = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        labs(title = " ")

ggsave("figures/fig_main_01a_pca_mf_all.pdf", loading_plot_mf, units = "in",
       width = 5, height = 5, dpi = 300)

## Scree plot Microfinance (remove grid lines)
scree_plot_mf = factoextra::fviz_eig(all_pca) +
  labs(title = " ") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(0,90))

# ggsave("figures/fig_aux_pca_mf_all_scree.pdf", scree_plot_mf, units = "in",
#        width = 5, height = 5, dpi = 300)

## Loadings Table Microfinance
loading_mf_table = rownames_to_column(all_pca$rotation %>% as.data.frame()) %>%
 rename("Network" = "rowname")

datasummary_df(loading_mf_table, output = "tables/tab_si_05_pca_loadings_mf.tex")

## without Jati and Distance -------

X_mf = df_edges_mf %>%
  select(all_of(netlist_mf)) %>%
  select(-jati, -distance_w) %>%
  rename("temple" = "templecompany")

all_pca = prcomp(X_mf, scale. = TRUE, center = TRUE)
all_pca$rotation = all_pca$rotation * -1

## Next we output the plots for each of the above PCA
options(ggrepel.max.overlaps = Inf)

## PC plot Microfinance
loading_plot_mf = factoextra::fviz_pca_var(all_pca, col.var = "red",
                                           axes = c(1, 2),
                                           labelsize = 3,
                                           arrowsize = 0.4,
                                           repel = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        labs(title = " ")

ggsave("figures/fig_si_03b_pca_mf_temple_loadings.pdf", loading_plot_mf, units = "in",
       width = 5, height = 5, dpi = 300)

## Scree plot Microfinance (remove grid lines)
scree_plot_mf = factoextra::fviz_eig(all_pca) +
  labs(title = " ") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(0,90))

ggsave("figures/fig_si_03a_pca_mf_temple_scree.pdf", scree_plot_mf, units = "in",
       width = 5, height = 5, dpi = 300)

## Loadings Table Microfinance
loading_mf_table = rownames_to_column(all_pca$rotation %>% as.data.frame()) %>%
 rename("Network" = "rowname")

# datasummary_df(
#   loading_mf_table,
#   output = "tables/tab_aux_pca_loadings_mf_without_jati_distance.tex"
# )

## Without Temple ---------

X_mf = df_edges_mf %>%
  select(all_of(netlist_mf)) %>%
  select(-templecompany, -jati, -distance_w)

all_pca = prcomp(X_mf, scale. = TRUE, center = TRUE)
all_pca$rotation = all_pca$rotation * -1

## Next we output the plots for each of the above PCA
options(ggrepel.max.overlaps = Inf)

## PC plot Microfinance
loading_plot_mf = factoextra::fviz_pca_var(all_pca, col.var = "red",
                                           axes = c(1, 2),
                                           labelsize = 3,
                                           arrowsize = 0.4,
                                           repel = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        labs(title = " ")

ggsave("figures/fig_main_01b_pca_mf_subset.pdf", loading_plot_mf, units = "in",
       width = 5, height = 5, dpi = 300)

## Scree plot Microfinance (remove grid lines)
scree_plot_mf = factoextra::fviz_eig(all_pca) +
  labs(title = " ") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(0,90))

# ggsave("figures/fig_aux_pca_mf_subset_scree.pdf", scree_plot_mf, units = "in",
#        width = 5, height = 5, dpi = 300)

## Loadings Table Microfinance
loading_mf_table = rownames_to_column(all_pca$rotation %>% as.data.frame()) %>%
 rename("Network" = "rowname")

# datasummary_df(
#   loading_mf_table,
#   output = "tables/tab_aux_pca_loadings_mf_without_temple_jati_distance.tex"
# )


####################### PCA RCT Villages ########################

## with jati

X_rct = df_edges_rct %>%
  select(all_of(netlist_rct)) %>%
      drop_na(jati)

rct_pca = prcomp(X_rct, scale. = TRUE, center = TRUE)
rct_pca$rotation = rct_pca$rotation * -1

loading_plot_rct = factoextra::fviz_pca_var(rct_pca, col.var = "red",
                                              axes = c(1, 2),
                                              labelsize = 3,
                                              arrowsize = 0.4,
                                              repel = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        labs(title = " ")

ggsave("figures/fig_main_01c_pca_rct_all.pdf", loading_plot_rct, units = "in",
       width = 5, height = 5, dpi = 300)

## Scree plot

scree_plot_rct = factoextra::fviz_eig(rct_pca) +
  labs(title = " ") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(0,90))

# ggsave("figures/fig_aux_pca_rct_all_scree.pdf", scree_plot_rct, units = "in",
#        width = 5, height = 5, dpi = 300)

loading_rct_table = rownames_to_column(rct_pca$rotation %>% as.data.frame()) %>%
 rename("Network" = "rowname")

datasummary_df(loading_rct_table, output = "tables/tab_si_06_pca_loadings_rct.tex")

## without jati -----------

X_rct = df_edges_rct %>%
  select(all_of(netlist_rct)) %>%
  select(-jati)

rct_pca = prcomp(X_rct, scale. = TRUE, center = TRUE)
rct_pca$rotation = rct_pca$rotation * -1

loading_plot_rct = factoextra::fviz_pca_var(rct_pca, col.var = "red",
                                              axes = c(1, 2),
                                              labelsize = 3,
                                              arrowsize = 0.4,
                                              repel = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        labs(title = " ")

ggsave("figures/fig_main_01d_pca_rct_subset.pdf", loading_plot_rct, units = "in",
       width = 5, height = 5, dpi = 300)

## Scree plot

scree_plot_rct = factoextra::fviz_eig(rct_pca) +
  labs(title = " ") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(0,90))

# ggsave("figures/fig_aux_pca_rct_subset_scree.pdf", scree_plot_rct, units = "in",
#        width = 5, height = 5, dpi = 300)

loading_rct_table = rownames_to_column(rct_pca$rotation %>% as.data.frame()) %>%
 rename("Network" = "rowname")

# datasummary_df(
#   loading_rct_table,
#   output = "tables/tab_aux_pca_loadings_rct_without_jati.tex"
# )
