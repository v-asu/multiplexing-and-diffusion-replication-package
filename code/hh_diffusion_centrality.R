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

# Analysis: What household characteristics predict diffusion centrality?
# Using RCT village data - Multiple network layers

# Load processed adjacency matrices
load("data/processed/rct_villages/rct_network_adjacency_layers.RData")

# Load household data
hh = read_dta("data/raw/rct_villages/rct_hh_covariates.dta")

# Load lambda_1 data
lambda_data = read.csv("data/processed/rct_villages/lambda_1_rct.csv")

# Load covariates (caste and occupation)
raw_covariates = list.files(path = "data/raw/rct_villages/hh_covariates_by_village",
                            pattern = "covariates",
                            full.names = TRUE)
raw_covariates = gtools::mixedsort(raw_covariates)
dat_covariates = map(raw_covariates, data.table::fread)
dat_covariates = data.table::rbindlist(dat_covariates)

# Clean covariates - select relevant columns and rename
dat_covariates = dat_covariates %>%
    select(villageid, hh_id, caste_category, occupation) %>%
    rename(hhid = hh_id, village = villageid) %>%
    mutate(
        # Clean caste: 1=SC, 2=ST, 3=OBC, 4=General, 5=Religious Minority
        caste = case_when(
            caste_category == 1 ~ "SC",
            caste_category == 2 ~ "ST",
            caste_category == 3 ~ "OBC",
            caste_category == 4 ~ "General",
            caste_category == 5 ~ "Religious Minority",
            TRUE ~ NA_character_
        ),
        # Recode occupation to 5 groups
        occupation_recoded = case_when(
            occupation %in% c(2, 3, 4, 5, 6) ~ "Agriculture",
            occupation %in% c(7, 10, 12, 15) ~ "Casual labor",
            occupation %in% c(8, 9, 17) ~ "Salaried/Services",
            occupation %in% c(1, 11, 13, 14) ~ "Self-employed/Business",
            occupation == 16 ~ "Housewife",
            TRUE ~ NA_character_
        )
    )

# Village IDs
village_ids = c(1:26, 28:39, 41:71)

# Network layers to analyze
netlist_rct = c("social", "kerorice", "advice", "decision")

# Function to calculate DC for all nodes
compute_DC = function(A, q) {
    n = nrow(A)
    g = graph_from_adjacency_matrix(A, mode = "max")
    diam = diameter(g, directed = FALSE)

    DC = numeric(n)
    A_power = diag(n)

    for (i in 0:diam) {
        DC = DC + as.numeric(q^i * rowSums(A_power))
        A_power = A_power %*% A
    }

    return(DC)
}

# Calculate DC for each layer and village
# Create dataframe for first network
net = netlist_rct[1]

net_dc = data.frame(village = integer(), vertex = integer(), dc_value = numeric())

for (idx in 1:length(village_ids)) {
    v = village_ids[idx]
    A = get(net)[[idx]]

    if (!is.null(A) && nrow(A) > 0) {
        lambda_1 = lambda_data[lambda_data$village == v, paste0("lambda_1_", net)]

        if (!is.na(lambda_1) && lambda_1 > 0) {
            q = 1/lambda_1
            dc = compute_DC(A, q)
            dc = dc / sum(dc)

            net_dc = bind_rows(net_dc, data.frame(
                village = v,
                vertex = 1:length(dc),
                dc_value = dc
            ))
        }
    }
}

df_dc = net_dc %>% rename(!!paste0("dc_", net) := dc_value)

# Process remaining networks
for (net in netlist_rct[-1]) {
    net_dc = data.frame(village = integer(), vertex = integer(), dc_value = numeric())

    for (idx in 1:length(village_ids)) {
        v = village_ids[idx]
        A = get(net)[[idx]]

        if (!is.null(A) && nrow(A) > 0) {
            lambda_1 = lambda_data[lambda_data$village == v, paste0("lambda_1_", net)]

            if (!is.na(lambda_1) && lambda_1 > 0) {
                q = 1/lambda_1
                dc = compute_DC(A, q)
                dc = dc / sum(dc)

                net_dc = bind_rows(net_dc, data.frame(
                    village = v,
                    vertex = 1:length(dc),
                    dc_value = dc
                ))
            }
        }
    }

    net_dc = net_dc %>% rename(!!paste0("dc_", net) := dc_value)
    df_dc = df_dc %>% left_join(net_dc, by = c("village", "vertex"))
}

# Clean household data
hh_clean = hh %>%
    mutate(
        village = villageid,
        own_hh = ifelse(hh_owned == 1, 1, ifelse(hh_owned < 0, NA, 0)),
        num_rooms = ifelse(num_rooms < 0, NA, num_rooms),
        leader = leader_dummy,
        electricity = have_electricity
    ) %>%
    select(hhid, village, vertex, own_hh, num_rooms,
           leader, electricity, seed_dummy) %>%
    left_join(dat_covariates %>% select(village, hhid, caste, occupation_recoded),
              by = c("village", "hhid"))

# Merge with DC
df = df_dc %>%
    left_join(hh_clean, by = c("village", "vertex"))

# Scale DC at the village level (z-score within each village)
for (net in netlist_rct) {
    dv = paste0("dc_", net)
    df = df %>% group_by(village) %>% mutate(!!dv := as.numeric(scale(.data[[dv]]))) %>% ungroup()
}

# Run regressions for each layer
# Store models using fixest::feols (works better with modelsummary)
models_by_layer = list()

for (net in netlist_rct) {
    dv = paste0("dc_", net)

    df_reg = df %>%
        filter(!is.na(.data[[dv]]), !is.na(caste), !is.na(occupation_recoded), !is.na(electricity),
               !is.na(num_rooms), !is.na(leader))

    # Use fixest::feols with village fixed effects and clustered SE
    m = feols(as.formula(paste0(dv, " ~ factor(caste) + factor(occupation_recoded) + electricity + num_rooms + own_hh + leader | village")),
              data = df_reg, vcov = list(~village))

    models_by_layer[[net]] = m
}

# Rename models for display
model_names = c(
    "social" = "Social",
    "kerorice" = "Kero/Rice",
    "advice" = "Advice",
    "decision" = "Decision",
    "union_link" = "Union"
)

# Coefficient rename map
coef_rename_map = c(
    "factor(caste)SC" = "Caste: SC",
    "factor(caste)ST" = "Caste: ST",
    "factor(caste)OBC" = "Caste: OBC",
    "factor(caste)Religious Minority" = "Caste: Religious Minority",
    "factor(occupation_recoded)Casual labor" = "Occupation: Casual labor",
    "factor(occupation_recoded)Salaried/Services" = "Occupation: Salaried/Services",
    "factor(occupation_recoded)Self-employed/Business" = "Occupation: Self-employed/Business",
    "factor(occupation_recoded)Housewife" = "Occupation: Housewife",
    "electricity" = "Has Electricity",
    "num_rooms" = "Number of Rooms",
    "own_hh" = "Owns Household",
    "leader" = "Village Leader"
)

# Get mean DV for each model
mean_dv = sapply(netlist_rct, function(net) {
    dv = paste0("dc_", net)
    df_reg = df %>% filter(!is.na(.data[[dv]]))
    mean(df_reg[[dv]], na.rm = TRUE)
})

# Create mean row
mean_row = data.frame(
    term = "Dep Var mean",
    setNames(as.list(mean_dv), netlist_rct)
)

# Save table using modelsummary
# Set option to use plain latex formatting
options("modelsummary_format_numeric_latex" = "plain")

# Generate and save table
modelsummary(models_by_layer,
            coef_omit = "factor\\(village\\)",
            coef_map = coef_rename_map,
            gof_omit = "IC|Log|Lik|Adj|Within|RMSE",
            statistic = c("std.error", "[{p.value}]"),
            title = "Determinants of Household Diffusion Centrality",
            notes = "Reference categories: Caste = General; Occupation = Agriculture. Standard errors clustered at the village level.",
            output = "latex") %>% 
    add_header_above(c(" " = 1, "Diffusion Centrality (std)" = length(netlist_rct))) %>%
    save_kable("tables/tab_si_10_hh_diffusion_centrality.tex")

cat("\nTable saved to tables/tab_si_10_hh_diffusion_centrality.tex\n")
