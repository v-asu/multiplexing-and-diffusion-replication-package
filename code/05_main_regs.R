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

## Main regressions for RCT Data

##################### Data ####################################
dat_calls_rct = haven::read_dta("data/raw/rct_villages/rct_village_diffusion_outcomes.dta") ## village level data
df_DC_seed = read.csv("data/processed/rct_villages/rct_village_seed_diffusion_centrality.csv") ## Sum of Seed set centrality

df_DC_seed = df_DC_seed %>% select(-sum_union_wtd)
no_jati_column_names = c("link_no_jati", "intersect_link_no_jati", "backbone_no_jati")
missing_no_jati_columns = setdiff(no_jati_column_names, names(df_DC_seed))
if (length(missing_no_jati_columns) > 0) {
  stop(
    "rct_village_seed_diffusion_centrality.csv is missing no-Jati columns: ",
    paste(missing_no_jati_columns, collapse = ", "),
    ". Run code/04_diffusion_centrality.R first."
  )
}

##################### Processing ###########################

## merging data
dat_calls_rct = dat_calls_rct %>%
  filter(villageid %in% c(1:26, 28:39, 41:71)) %>% ## village 40 dropped (implementation issues)
  mutate(calls_per_HH = CallsReceived/num_hh_random,
         dummy_3_seeds = ifelse(num_seeds == 3, 1, 0)) %>%
  select(CallsReceived, dummy_3_seeds, calls_per_HH, num_hh_random, villageid)

reg_data_rct = dat_calls_rct %>%
  left_join(df_DC_seed, by = c("villageid" = "village")) %>%
  mutate(num_hh_random_sq = (num_hh_random)^2,
         num_hh_random_cb = (num_hh_random)^3)

## Standardizing exog variables
reg_data_rct = reg_data_rct %>% mutate(across(.cols = -c(CallsReceived, dummy_3_seeds, calls_per_HH, villageid),
                                               \(x) as.numeric(scale(x))))


############### Individual Regs ############################

var_names = names(reg_data_rct)

mod_in_DC_seed = data.frame(dep_var = rep("CallsReceived"),
                            exog_vars = var_names[grep("sum", var_names)],
                            controls = I(rep(list(grep("dummy_3_seeds|num_hh", var_names, value = T)),
                             length(var_names[grep("sum", var_names)]))))

mod_out_DC_seed = pmap(mod_in_DC_seed, model_out, dat = reg_data_rct)
names(mod_out_DC_seed) =  1:length(mod_out_DC_seed)

## Dep var mean
mean_val = reg_data_rct %>% filter(villageid != 62) %>% pull(CallsReceived) %>% mean() %>% round(3)

mean_val_df = tribble(~term, ~"(1)", ~"(2)", ~"(3)", ~"(4)", ~"(5)", ~"(6)", ~"(7)", ~"(8)", ~"(9)",
                      "Dep Var mean", mean_val, mean_val,mean_val, mean_val, mean_val, mean_val, mean_val, mean_val, mean_val)

mean_val_df = data.frame(
  term = "Dep Var mean",
  matrix(mean_val, nrow = 1, ncol = length(mod_out_DC_seed))
)

modelsummary(mod_out_DC_seed, coef_omit = 'num_hh_random|(Intercept)|dummy_3_seeds',
             gof_omit = '[^Num.Obs.|R2]',
             coef_rename = c("sum_union_link" = "Union",
                          "sum_intersect_link" = "Intersection",
                          "sum_advice" = "Advice",
                          "sum_kerorice" = "Kero/Rice",
                          "sum_social" = "Social",
                          "sum_info" = "Information",
                          "sum_backbone" = "Backbone",
                          "sum_jati" = "Jati",
                          "sum_decision" = "Decision",
                          "sum_raw_wtd" = "Total Links"),
             statistic = c("std.error", "[{p.value}]"),
             title = "Seed Set Diffusion Centrality",
             output = "latex",
             add_rows = mean_val_df) %>%
  add_header_above(c(" " = 1, "No. Calls Received" = length(mod_out_DC_seed))) %>%
  footnote(general = c("Robust Std.Err are given in paranthesis, while p-values are given in square brackets.",
             "Controls added for number of Households and its powers, and a dummy for number of seeds in the village",
             "Exog variables are the sum of Diffusion Centrality for seeds in each village for the layer",
             "Exog variables have been standardized")) %>%
  save_kable("tables/tab_si_02_seed_dc_by_layer.tex")

############### Individual Regs: no-jati composite layers ############################

df_DC_seed_no_jati_composites = df_DC_seed %>%
  select(-sum_union_link, -sum_intersect_link, -sum_backbone)

reg_data_rct_no_jati_composites = dat_calls_rct %>%
  left_join(df_DC_seed_no_jati_composites, by = c("villageid" = "village")) %>%
  mutate(num_hh_random_sq = (num_hh_random)^2,
         num_hh_random_cb = (num_hh_random)^3)

reg_data_rct_no_jati_composites = reg_data_rct_no_jati_composites %>%
  mutate(across(.cols = -c(CallsReceived, dummy_3_seeds, calls_per_HH, villageid),
                \(x) as.numeric(scale(x))))

var_names_no_jati_composites = names(reg_data_rct_no_jati_composites)
sum_vars_no_jati_composites = var_names_no_jati_composites[grep("sum", var_names_no_jati_composites)]
exog_vars_no_jati_composites = c(
  setdiff(sum_vars_no_jati_composites, "sum_raw_wtd"),
  no_jati_column_names,
  "sum_raw_wtd"
)

mod_in_DC_seed_no_jati_composites = data.frame(
  dep_var = rep("CallsReceived"),
  exog_vars = exog_vars_no_jati_composites,
  controls = I(rep(
    list(grep("dummy_3_seeds|num_hh", var_names_no_jati_composites, value = TRUE)),
    length(exog_vars_no_jati_composites)
  ))
)

mod_out_DC_seed_no_jati_composites = pmap(
  mod_in_DC_seed_no_jati_composites,
  model_out,
  dat = reg_data_rct_no_jati_composites
)
names(mod_out_DC_seed_no_jati_composites) = seq_along(mod_out_DC_seed_no_jati_composites)

mean_val_no_jati_composites = reg_data_rct_no_jati_composites %>%
  filter(villageid != 62) %>%
  pull(CallsReceived) %>%
  mean() %>%
  round(3)

mean_val_df_no_jati_composites = data.frame(
  term = "Dep Var mean",
  matrix(mean_val_no_jati_composites, nrow = 1, ncol = length(mod_out_DC_seed_no_jati_composites))
)

modelsummary(mod_out_DC_seed_no_jati_composites, coef_omit = 'num_hh_random|(Intercept)|dummy_3_seeds',
             gof_omit = '[^Num.Obs.|R2]',
             coef_rename = c("sum_union_link" = "Union",
                          "sum_intersect_link" = "Intersection",
                          "link_no_jati" = "Union",
                          "intersect_link_no_jati" = "Intersection",
                          "sum_advice" = "Advice",
                          "sum_kerorice" = "Kero/Rice",
                          "sum_social" = "Social",
                          "sum_info" = "Information",
                          "sum_backbone" = "Backbone",
                          "backbone_no_jati" = "Backbone",
                          "sum_jati" = "Jati",
                          "sum_decision" = "Decision",
                          "sum_raw_wtd" = "Total Links"),
             statistic = c("std.error", "[{p.value}]"),
             title = "Seed Set Diffusion Centrality (No-jati aggregate layers)",
             output = "latex",
             add_rows = mean_val_df_no_jati_composites) %>%
  add_header_above(c(" " = 1, "No. Calls Received" = length(mod_out_DC_seed_no_jati_composites))) %>%
  footnote(general = c("Robust Std.Err are given in paranthesis, while p-values are given in square brackets.",
             "Controls added for number of Households and its powers, and a dummy for number of seeds in the village",
             "Union, intersection, and backbone layers are constructed without jati; the other layers are unchanged.",
             "Exog variables are the sum of Diffusion Centrality for seeds in each village for the layer",
             "Exog variables have been standardized")) %>%
  save_kable("tables/tab_si_03_seed_dc_no_jati_aggregates.tex")

################## Puffer Lasso ############################

puffer_DC = list()
puffer_DC[["1"]] =  post_lasso_ols_glmnet(
  exog_selection = "sum",
  lambda_quantile = 0.95,
  data = reg_data_rct,
  var_names = var_names
)

modelsummary(puffer_DC, coef_omit = '(Intercept)|dummy_3_seeds',
             gof_omit = '[^Num.Obs.|R2]',
             coef_map = c("dummy_3_seeds" = "No. Seeds Dummy",
                          "sum_advice" = "Advice",
                          "sum_backbone" = "Backbone",
                          "sum_kerorice" = "Kero/Rice"),
             statistic = c("std.error", "[{p.value}]"),
             title = "Post Puffer LASSO OLS: Seed Set Diffusion Centrality",
             output = "latex",
             add_rows = mean_val_df[, 1:2]) %>%
  add_header_above(c(" " = 1, "No. Calls Received" = 1)) %>%
  save_kable("tables/tab_si_04_post_lasso_seed_dc.tex")

################## OLS Plot ################################

dc_plot_df = map_dfr(c(0.9, 0.95), function(x) {
  modelplot(mod_out_DC_seed,
            coef_omit = 'num_hh_random|(Intercept)|dummy_3_seeds',
            coef_map = c("sum_union_link" = "Union",
            "sum_union_wtd" = "Total Links",
            "sum_raw_wtd" = "Total Links",
                         "sum_intersect_link" = "Intersection",
                         "sum_advice" = "Advice",
                         "sum_kerorice" = "Kerorice",
                         "sum_social" = "Social",
                         "sum_backbone" = "Backbone",
                         "sum_jati" = "Jati",
                         "sum_decision" = "Decision",
                         "sum_info" = "Information"
                         ),
            conf_level = x,
            draw = FALSE
            ) %>%
    mutate(.width = x)
})


ols_plot = ggplot(dc_plot_df, aes(
  y = term, x = estimate,
  xmin = conf.low, xmax = conf.high,
  color = model,
  point_alpha = 0.2
)) +
  geom_pointinterval(
    position = "identity",
    interval_size_range = c(0.5,1.2),
    point_size = 2,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  scale_y_discrete(limits=rev) +
  geom_text(aes(label = c(paste0("p val = ", dc_plot_df %>% filter(.width == 0.95) %>% pull(p.value) %>% round(3)),
                         rep(" ", length(mod_out_DC_seed)))), show.legend = FALSE,
            nudge_y = 0.3, size = 3) +
  labs(x = "Estimate", y = "Regression", title = "Seed Set Diffusion Centrality") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_npg()

ggsave("figures/fig_main_02_seed_dc_by_layer.pdf", ols_plot, units = "in",
       width = 9, height = 5, dpi = 300)

################## Lasso Plot ################################

# Create a named vector for layer renaming
layer_names = c(
  "sum_advice" = "Advice",
  "sum_decision" = "Decision",
  "sum_social" = "Social",
  "sum_kerorice" = "Kero/Rice",
  "sum_backbone" = "Backbone",
  "sum_jati" = "Jati",
  "sum_info" = "Information",
  "sum_union_wtd" = "Total Links",
  "sum_intersect_link" = "Intersection",
  "sum_union_link" = "Union",
  "sum_raw_wtd" = "Total Links"
)

puffer_data = puffer_N_transform(data = reg_data_rct, index_dep = 1, index_exog = grep("sum", var_names))
lasso_glmnet = glmnet(x = puffer_data[, -1], y = puffer_data[, 1], nlambda = 1000)

## We want to get the sequence of rsquared as we keep adding exogenous variables
rhs_list = c(coef(lasso_glmnet, lasso_glmnet$lambda[1]) %>% rownames())
rhs_list = unique(map(lasso_glmnet$lambda, ~rhs_list[which(coef(lasso_glmnet, .x) != 0)][-1]))[-1]

## regress all of the rhs list one by one and get the rsquared
rsq_list_post_ols = map_dbl(rhs_list, ~lm_robust(reformulate(response = "CallsReceived",
                                   termlabels = .x),
                       data = reg_data_rct)$r.squared)

## get the support list
support_list = map(lasso_glmnet$lambda, lasso_plot_vals)

lasso_layers = unique(unlist(rhs_list))

df_lasso = data.frame(penalty = lasso_glmnet$lambda[-1])

# Add columns for each layer
for (layer in lasso_layers) {
  df_lasso[[layer]] = 0
}

# Populate df_lasso
support_list = map(lasso_glmnet$lambda[-1], lasso_plot_vals)
for (i in seq_along(df_lasso$penalty)) {
  for (layer in lasso_layers) {
    df_lasso[i, layer] = as.integer(any(str_detect(support_list[[i]], layer)))
  }
}

# Remove duplicate rows and the row where all are dropped
df_lasso = df_lasso %>%
  distinct(across(all_of(lasso_layers)), .keep_all = TRUE) %>%
  filter(rowSums(select(., all_of(lasso_layers))) > 0)

df_lasso$rsq_post_ols = rsq_list_post_ols

# Reshape data
df_lasso_long = df_lasso %>%
  pivot_longer(cols = all_of(lasso_layers), names_to = "layer", values_to = "picked") %>%
  mutate(picked = ifelse(picked == 1, "Yes", "No"))

# Rename layers
df_lasso_long = df_lasso_long %>%
  mutate(layer = recode(layer, !!!layer_names))

# Order layers based on their selection at penalty levels
layer_order = df_lasso_long %>%
  filter(picked == "Yes") %>%
  group_by(layer) %>%
  summarise(max_penalty = max(penalty)) %>%
  arrange(desc(max_penalty)) %>%
  pull(layer)

df_lasso_long$layer = factor(df_lasso_long$layer, levels = rev(layer_order))

## Raw plot without annotations
lasso_plot = ggplot(df_lasso_long %>% filter(picked == "Yes"), aes(x = penalty, y = layer)) +
  geom_point(shape = 4) +
  scale_x_reverse() +
  theme_bw() +
  labs(x = "Penalty level (lambda)", y = "Layer",
       title = "'x' implies that the layer was selected")

ggsave("figures/fig_si_02_lasso_layer_selection.pdf", lasso_plot, width = 10, height = 5, units = "in")

####### F test table

f_list = map2_dfr(rhs_list[[1]], rhs_list, ~F_calc(.x, .y, reg_data_rct))
f_marginal = map2_dfr(rhs_list[1:(length(rhs_list)-1)], rhs_list[2:length(rhs_list)], ~F_calc(.x, .y, reg_data_rct))

f_test_table = data.frame(
  layer = lasso_layers,
  k = seq_along(lasso_layers),
  raw_r_squared = rsq_list_post_ols[seq_along(lasso_layers)],
  F = f_list$F[seq_along(lasso_layers)],
  p_val = f_list$p_val[seq_along(lasso_layers)],
  F_marginal = c(NA, f_marginal$F),
  p_val_marginal = c(NA, f_marginal$p_val)
)

# Rename layers in the F-test table
f_test_table$layer = recode(f_test_table$layer, !!!layer_names)

f_test_table[-1] = round(f_test_table[-1], 3)

## write it as a latex table
kable(f_test_table %>%
  set_names(c("layer", "df", "R.sq.", "F-stat", "p-val", "F-stat marginal", "p-val marginal")), format = "latex", booktabs = T, digits = 3, caption = "F-test for the layers", na = "") %>%
  save_kable(file = "tables/tab_main_01_f_test_layers.tex")

###### New F-test table version (without composite layers) #######

reg_data_rct = reg_data_rct %>%
  select(-sum_union_link, -sum_intersect_link, -sum_backbone, -sum_raw_wtd)

var_names = names(reg_data_rct)

puffer_data = puffer_N_transform(data = reg_data_rct, index_dep = 1, index_exog = grep("sum", var_names))
lasso_glmnet = glmnet(x = puffer_data[, -1], y = puffer_data[, 1], nlambda = 1000)

## We want to get the sequence of rsquared as we keep adding exogenous variables
rhs_list = c(coef(lasso_glmnet, lasso_glmnet$lambda[1]) %>% rownames())
rhs_list = unique(map(lasso_glmnet$lambda, ~rhs_list[which(coef(lasso_glmnet, .x) != 0)][-1]))[-1]

f_list = map2_dfr(rhs_list[[1]], rhs_list, ~F_calc(.x, .y, reg_data_rct))
f_marginal = map2_dfr(rhs_list[1:(length(rhs_list)-1)], rhs_list[2:length(rhs_list)], ~F_calc(.x, .y, reg_data_rct))

rsq_list_post_ols = map_dbl(rhs_list, ~lm_robust(reformulate(response = "CallsReceived",
                                   termlabels = .x),
                       data = reg_data_rct)$r.squared)

lasso_layers_short = unique(unlist(rhs_list))

f_test_table_short = data.frame(
  layer = lasso_layers_short,
  k = seq_along(lasso_layers_short),
  raw_r_squared = rsq_list_post_ols,
  F = f_list$F,
  p_val = f_list$p_val,
  F_marginal = c(NA, f_marginal$F),
  p_val_marginal = c(NA, f_marginal$p_val)
)

# Rename layers in the F-test table
f_test_table_short$layer = recode(f_test_table_short$layer, !!!layer_names)

f_test_table_short[-1] = round(f_test_table_short[-1], 3)

## write it as a latex table
kable(f_test_table_short %>%
  set_names(c("layer", "df", "R.sq.", "F-stat", "p-val", "F-stat marginal", "p-val marginal")), format = "latex", booktabs = T, digits = 3, caption = "F-test for the layers", na = "") %>%
  save_kable(file = "tables/tab_si_07_f_test_basic_layers.tex")

############### END ########################
