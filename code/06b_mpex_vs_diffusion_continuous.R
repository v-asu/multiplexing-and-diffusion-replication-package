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

## Regression with multiplexing x DC terciles + controls x DC terciles
## Plus marginal effects and three-panel plot

##################### Data ####################################

df_DC_seed = read.csv("data/processed/rct_villages/rct_village_seed_diffusion_centrality.csv")
dat_calls_rct = haven::read_dta("data/raw/rct_villages/rct_village_diffusion_outcomes.dta")
load("data/processed/rct_villages/rct_network_adjacency_layers.RData")

#################### Processing ###########################

total_links = vector(mode = "list", length = 71)
total_people = vector(mode = "list", length = 71)
M_index = vector(mode = "list", length = 71)

L = 4

village_list = c(1:26, 28:39, 41:71)

for (i in village_list) {
    if (is.null(social[[i]]) || is.null(kerorice[[i]]) ||
        is.null(advice[[i]]) || is.null(decision[[i]]) ||
        is.null(union_link[[i]])) {
        total_links[[i]] = NA
        total_people[[i]] = NA
        M_index[[i]] = NA
        next
    }

    total_links[[i]] = (rowSums(social[[i]]) + rowSums(kerorice[[i]]) +
                     + rowSums(advice[[i]]) + rowSums(decision[[i]]))/L

    combined = social[[i]] + kerorice[[i]] + decision[[i]] + advice[[i]]
    total_people[[i]] = sapply(1:nrow(combined), function(x) length(which(combined[x, ] > 0)))

    M_index[[i]] = total_links[[i]] / total_people[[i]]
    M_index[[i]][total_people[[i]] == 0] = NA
}

m_i_list = map(M_index[village_list], ~mean(.x, na.rm = TRUE))
m_i = unlist(m_i_list)
names(m_i) = village_list

avg_degree_list = map(total_links[village_list], ~mean(.x, na.rm = TRUE))
avg_degree = unlist(avg_degree_list)
names(avg_degree) = village_list

dat_calls_rct = dat_calls_rct %>%
  filter(villageid %in% c(1:26, 28:39, 41:71)) %>%
  mutate(calls_per_HH = CallsReceived/num_hh_random,
         dummy_3_seeds = ifelse(num_seeds == 3, 1, 0)) %>%
  select(CallsReceived, dummy_3_seeds, calls_per_HH, num_hh_random, villageid)

reg_data_rct = dat_calls_rct %>%
  left_join(df_DC_seed, by = c("villageid" = "village"))

reg_data_rct = reg_data_rct %>% arrange(villageid)

reg_data_rct$M_i = m_i[as.character(reg_data_rct$villageid)]
reg_data_rct$avg_deg = avg_degree[as.character(reg_data_rct$villageid)]

reg_data_rct = reg_data_rct %>%
  mutate(sum_advice_std = as.vector(scale(sum_advice)),
         M_i_std = as.vector(scale(M_i)),
         avg_deg_std = as.vector(scale(avg_deg)))

# Tertile dummies on seed centrality (advice layer)
reg_data_rct$DC_tercile = cut(reg_data_rct$sum_advice,
                              breaks = quantile(reg_data_rct$sum_advice, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                              labels = c("Low", "Medium", "High"),
                              include.lowest = TRUE)

# Set Low as reference category
reg_data_rct$DC_tercile = relevel(reg_data_rct$DC_tercile, ref = "Low")

# Drop rows with NA in relevant variables
reg_data_rct = reg_data_rct %>%
  filter(!is.na(M_i),
         !is.na(calls_per_HH),
         !is.na(sum_advice),
         !is.na(avg_deg),
         !is.na(dummy_3_seeds))

#################### Table: Multiplexing x DC Terciles ###########################

# Helper function to format marginal effect with p-value
format_me_pval <- function(est, pval) {
  sprintf("%.3f (p=%.3f)", est, pval)
}

# Regression: Multiplexing x DC terciles WITH control interactions (but only show M in table)
reg_table = lm_robust(calls_per_HH ~ M_i_std * DC_tercile + dummy_3_seeds * DC_tercile + avg_deg_std * DC_tercile,
                      data = reg_data_rct)

# Get marginal effects for the table
me_table = slopes(reg_table, variables = "M_i_std", by = "DC_tercile", vcov = TRUE)

# Add marginal effects as additional rows
me_add_tercile <- data.frame(
  term = c("ME: Low DC", "ME: Medium DC", "ME: High DC"),
  "(1)" = c(format_me_pval(me_table$estimate[1], me_table$p.value[1]),
            format_me_pval(me_table$estimate[2], me_table$p.value[2]),
            format_me_pval(me_table$estimate[3], me_table$p.value[3]))
)

# Footnote
note_tercile <- "Seed set diffusion centrality (DC) is split into three terciles (Low, Medium, High). Multiplexing is a continuous measure. Controls included and interacted with DC terciles. Marginal effects of Multiplexing for each DC tercile are reported at the bottom."

# Configure modelsummary for kableExtra
options(modelsummary_factory_latex = 'kableExtra')

reg_list = list()
reg_list[["(1)"]] = reg_table

modelsummary(reg_list,
             coef_omit = '(Intercept)|dummy_3_seeds.*|avg_deg.*',
             gof_omit = '[^Num.Obs]',
             add_rows = me_add_tercile,
             coef_rename = c("M_i_std:DC_tercileMedium" = "Multiplexing x Medium DC",
                             "M_i_std:DC_tercileHigh" = "Multiplexing x High DC",
                             "DC_tercileMedium" = "Medium DC (Tertile)",
                             "DC_tercileHigh" = "High DC (Tertile)",
                             "M_i_std" = "Multiplexing"),
             statistic = c("std.error", "[{p.value}]"),
             note = note_tercile,
             output = "latex") %>%
  save_kable("tables/tab_si_08_multiplexing_by_dc_tercile.tex")

cat("\nTable saved to tables/tab_si_08_multiplexing_by_dc_tercile.tex\n")

#################### Plot: Three-Panel ###########################

# Create DC terciles for plotting (convert to labels for facets)
reg_data_rct$DC_tercile_plot = factor(reg_data_rct$DC_tercile,
                                      levels = c("Low", "Medium", "High"),
                                      labels = c("Low Seed Set Centrality", "Medium Seed Set Centrality", "High Seed Set Centrality"))

p3 = ggplot(reg_data_rct, aes(x = M_i, y = calls_per_HH)) +
  geom_point(size = 3, alpha = 0.7, color = "#333333") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray90", linewidth = 0.8) +
  facet_wrap(~DC_tercile_plot, nrow = 1, dir = "v") +
  labs(
    x = "Multiplexing Index (M)",
    y = "Calls per Household"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", hjust = 0),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

ggsave("figures/fig_si_07_multiplexing_by_dc_tercile.pdf", p3, width = 10, height = 6, dpi = 300)
