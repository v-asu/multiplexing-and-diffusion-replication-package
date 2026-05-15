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

## heatmaps for correlation across layers
# Loading Data ----------------------------------------------------------------------

# Microfinance Data --------------------------------------------
df_edges_mf = readRDS("data/processed/mf_villages/mf_hh_pair_edges.rds") ## 1932224 obs

# RCT Data --------------------------------------------------
df_edges_rct =  readRDS("data/processed/rct_villages/rct_hh_pair_edges.rds") ## 1450245 obs

## Microfinance Villages Heat Maps

X_edges_mf_2 = df_edges_mf %>% 
  select(social, kerorice, advice, decision, money, medic, templecompany, jati, distance_w) %>%
  rename("distance" = "distance_w") %>%
  rename("temple" = "templecompany") %>% 
  drop_na(distance, jati) %>% 
  as.matrix()

cor_mf_2 = cor(X_edges_mf_2) %>% 
  reshape2::melt() %>% 
  filter(Var1 != Var2)

get_text_color = function(value) {
  if (value > 0.5) {
    return("white")
  } else {
    return("black")
  }
}

cor_mf_2$text_color = sapply(cor_mf_2$value, get_text_color)

heatmap_mf = ggplot(data = cor_mf_2, 
                    aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", direction = -1)+
  geom_text(aes(Var2, Var1, label = value %>% round(3), color = text_color), size = 4) +
  scale_color_identity() + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
         panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(title = "", x = "", y = "")

ggsave("figures/fig_si_01a_layer_correlations_mf.pdf", heatmap_mf, units = "in",
       width = 8, height = 5, dpi = 300)

## RCT Villages Heat map

X_edges_rct_2 = df_edges_rct %>% 
  select(social, kerorice, advice, decision, jati) %>% 
  drop_na(jati) %>% 
  as.matrix()

cor_rct_2 = cor(X_edges_rct_2) %>% 
  reshape2::melt() %>% 
  filter(Var1 != Var2)

cor_rct_2$text_color = sapply(cor_rct_2$value, get_text_color)

heatmap_rct = ggplot(data = cor_rct_2, 
                     aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", direction = -1)+
  geom_text(aes(Var2, Var1, label = value %>% round(3), color = text_color), size = 4) +
  scale_color_identity() + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(title = "", x = "", y = "", legend = "")

ggsave("figures/fig_si_01b_layer_correlations_rct.pdf", heatmap_rct, units = "in",
       width = 8, height = 5, dpi = 300)

