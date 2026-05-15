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

# Gender-Stratified Multiplexing Index
# M_i_f: multiplexing with female contacts only
# M_i_m: multiplexing with male contacts only

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
 
# Create collapsed networks -------------------------------------------------
social = list(); kerorice = list(); money = list(); advice = list(); decision = list()

for (i in 1:75) {
    social[[i]] = (visitgo[[i]] + visitcome[[i]] + rel[[i]] + nonrel[[i]] > 0)*1
    kerorice[[i]] = (keroricego[[i]] + keroricecome[[i]] > 0)*1
    money[[i]] = (borrowmoney[[i]] + lendmoney[[i]] > 0)*1
    advice[[i]] = giveadvice[[i]]
    decision[[i]] = helpdecision[[i]]
}

main_layers = c("social", "kerorice", "money", "advice", "decision", "medic", "templecompany")
L = length(main_layers)

# Build gender vector for each village -------------------------------------
village_ids = c(1:12, 14:21, 23:77)

dat_keys = map2_df(
    social,
    village_ids,
    ~data.frame(
        village = .y,
        v_id = seq_len(nrow(.x)),
        UniqueID = as.numeric(rownames(.x))
    )
)

df_covs_clean = df_covs %>%
    mutate(village = floor(newHHID / 1000))

dat_keys = dat_keys %>%
    left_join(df_covs_clean %>% select(village, UniqueID, respgender0_3),
              by = c("village", "UniqueID"))

# Create gender vector (1 = female, 0 = male) for each village
gender_by_village = list()

for (idx in seq_along(village_ids)) {
    v = village_ids[idx]
    sub = dat_keys %>% filter(village == v) %>% arrange(v_id)
    gender = ifelse(sub$respgender0_3 == 2, 1, 0)  # 2 = female -> 1, 1 = male -> 0
    gender[is.na(sub$respgender0_3)] = NA
    gender_by_village[[v]] = gender
}

# Calculate M_i_f and M_i_m -----------------------------------------------

village_ids <- c(1:12, 14:21, 23:77)
stopifnot(length(village_ids) == length(social))

M_f <- vector("list", length = length(village_ids))
M_m <- vector("list", length = length(village_ids))

for (idx in seq_along(village_ids)) {
  v <- village_ids[idx]
  n <- nrow(social[[idx]])

  gender <- gender_by_village[[v]]

  if (is.null(gender) || length(gender) != n || all(is.na(gender))) {
    M_f[[idx]] <- rep(NA, n)
    M_m[[idx]] <- rep(NA, n)
    next
  }

  # alter-gender mask repeated across rows (works for A_ij ties to j)
  female_mat <- matrix(gender, nrow = n, ncol = n, byrow = TRUE)
  male_mat   <- 1 - female_mat

  female_mat[is.na(female_mat)] <- 0
  male_mat[is.na(male_mat)] <- 0

  union_female <- matrix(0, n, n)
  union_male   <- matrix(0, n, n)

  sum_deg_female <- numeric(n)
  sum_deg_male   <- numeric(n)

  for (layer in main_layers) {
    A <- get(layer)[[idx]]

    deg_female <- rowSums(A * female_mat)
    deg_male   <- rowSums(A * male_mat)

    sum_deg_female <- sum_deg_female + deg_female
    sum_deg_male   <- sum_deg_male + deg_male

    union_female <- union_female + (A * female_mat > 0) * 1
    union_male   <- union_male   + (A * male_mat   > 0) * 1
  }

  deg_union_female <- rowSums(union_female > 0)
  deg_union_male   <- rowSums(union_male > 0)

  M_f[[idx]] <- ifelse(deg_union_female > 0, (sum_deg_female / L) / deg_union_female, NA)
  M_m[[idx]] <- ifelse(deg_union_male   > 0, (sum_deg_male   / L) / deg_union_male,   NA)
}

# Calculate edge counts, densities, and node-level stats per village ----------
edge_summary <- data.frame()

for (idx in seq_along(village_ids)) {
  v <- village_ids[idx]
  n <- nrow(social[[idx]])
  gender <- gender_by_village[[v]]

  if (is.null(gender) || length(gender) != n || all(is.na(gender))) next

  # Build union matrix
  union_mat <- matrix(0, n, n)
  for (layer in main_layers) {
    A <- get(layer)[[idx]]
    union_mat <- union_mat + (A > 0) * 1
  }
  union_mat <- (union_mat > 0) * 1

  # Gender matrices
  female_mat <- matrix(gender, nrow = n, ncol = n, byrow = TRUE)
  male_mat <- 1 - female_mat
  female_mat[is.na(female_mat)] <- 0
  male_mat[is.na(male_mat)] <- 0

  # Count actual edges (undirected - upper triangle only)
  upper_idx <- upper.tri(union_mat, diag = FALSE)

  n_ff <- sum(union_mat * female_mat * t(female_mat) * upper_idx, na.rm = TRUE)
  n_mm <- sum(union_mat * male_mat * t(male_mat) * upper_idx, na.rm = TRUE)
  n_cross <- sum(union_mat * (female_mat * t(male_mat) + male_mat * t(female_mat)) * upper_idx, na.rm = TRUE)

  # Count potential edges by type
  n_female <- sum(!is.na(gender) & gender == 1)
  n_male <- sum(!is.na(gender) & gender == 0)

  pot_ff <- n_female * (n_female - 1) / 2
  pot_mm <- n_male * (n_male - 1) / 2
  pot_cross <- n_female * n_male

  # Densities
  dens_ff <- ifelse(pot_ff > 0, n_ff / pot_ff, NA)
  dens_mm <- ifelse(pot_mm > 0, n_mm / pot_mm, NA)
  dens_cross <- ifelse(pot_cross > 0, n_cross / pot_cross, NA)

  # Node-level: fraction with no links, no cross-gender links
  deg_total <- rowSums(union_mat)
  deg_female_alter <- rowSums(union_mat * female_mat)
  deg_male_alter <- rowSums(union_mat * male_mat)

  is_female <- !is.na(gender) & gender == 1
  is_male <- !is.na(gender) & gender == 0

  frac_female_no_links <- ifelse(sum(is_female) > 0, mean(is_female & deg_total == 0), NA)
  frac_female_no_male_links <- ifelse(sum(is_female) > 0, mean(is_female & deg_male_alter == 0), NA)
  frac_male_no_links <- ifelse(sum(is_male) > 0, mean(is_male & deg_total == 0), NA)
  frac_male_no_female_links <- ifelse(sum(is_male) > 0, mean(is_male & deg_female_alter == 0), NA)

  edge_summary <- rbind(edge_summary, data.frame(
    village = v,
    n_ff = n_ff, n_mm = n_mm, n_cross = n_cross,
    pot_ff = pot_ff, pot_mm = pot_mm, pot_cross = pot_cross,
    dens_ff = dens_ff, dens_mm = dens_mm, dens_cross = dens_cross,
    frac_female_no_links = frac_female_no_links,
    frac_female_no_male_links = frac_female_no_male_links,
    frac_male_no_links = frac_male_no_links,
    frac_male_no_female_links = frac_male_no_female_links
  ))
}

# Build DF aligned to matrix indices, then attach true village id
df_M_gender <- map2(M_f, M_m, ~data.frame(M_f = .x, M_m = .y, adjmatrix_key = seq_along(.x))) %>%
  bind_rows(.id = "idx") %>%
  mutate(idx = as.integer(idx),
         village = village_ids[idx]) %>%
  select(village, idx, adjmatrix_key, M_f, M_m)

# Merge with covariates
df_M_gender = df_M_gender %>%
    left_join(dat_keys %>% select(village, v_id, respgender0_3),
              by = c("village", "adjmatrix_key" = "v_id"))

df_M_gender = df_M_gender %>%
    filter(!is.na(respgender0_3)) %>%
    mutate(gender = ifelse(respgender0_3 == 1, "Male", "Female"))

# Density plot --------------------------------------------------------------
# Reshape to long format for ggplot2
df_density <- df_M_gender %>%
    select(gender, M_f, M_m) %>%
    pivot_longer(cols = c(M_f, M_m), names_to = "contact_gender", values_to = "multiplexing") %>%
    mutate(
        ego_abbrev = ifelse(tolower(gender) == "female", "F", "M"),
        contact_abbrev = ifelse(contact_gender == "M_f", "F", "M"),
        label = paste0(ego_abbrev, contact_abbrev)
    )

# Create density plot
p <- ggplot(df_density, aes(x = multiplexing, color = label, linetype = label)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    xlim(0, 1) +
    scale_fill_manual(
        values = c(
            "FF" = "#E69F00",
            "FM" = "#56B4E9",
            "MF" = "#009E73",
            "MM" = "#CC79A7"
        )
    ) +
    scale_color_manual(
        values = c(
            "FF" = "#E69F00",
            "FM" = "#56B4E9",
            "MF" = "#009E73",
            "MM" = "#CC79A7"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "FF" = "solid",
            "FM" = "dashed",
            "MF" = "solid",
            "MM" = "dashed"
        )
    ) +
    labs(
        title = " ",
        x = "Multiplexing Index",
        y = "Density",
        fill = "link type",
        color = "link type",
        linetype = "link type"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

# Save plot
ggsave("figures/fig_si_06_same_cross_gender_multiplexing.pdf", p, width = 8, height = 8, dpi=300)
