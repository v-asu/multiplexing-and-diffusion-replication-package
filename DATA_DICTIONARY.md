# Data Dictionary

This file describes the raw inputs used by the replication code. Paths are
relative to the `replication_package/` root.

## Microfinance Inputs

| Path | Unit | Description |
|---|---|---|
| `data/raw/mf_villages/mf_network_adjacency_layers.mat` | Village network matrices | Household network layers for the microfinance villages. |
| `data/raw/mf_villages/mf_hh_covariates.dta` | Household | Household covariates used to construct same-caste and same-jati pair indicators. |
| `data/raw/mf_villages/mf_hh_vertex_crosswalk_wave2.rds` | Household | Vertex-to-household key with `village`, `w2_vertex`, and `newhhid`. |
| `data/raw/mf_villages/mf_hh_pair_distance_winsorized.rds` | Household pair | Pairwise distance input with `village`, `id1`, `id2`, and `distance_w`; used as the distance measure in the analysis. |
| `data/raw/mf_villages/gender/mf_individual_wave1.rds` | Individual network data | Wave 1 individual-level network object. |
| `data/raw/mf_villages/gender/mf_individual_wave2.rds` | Individual network data | Wave 2 individual-level network object used for gender analyses. |
| `data/raw/mf_villages/mf_ind_covariates_wave1.dta` | Individual | Wave 1 individual IDs and gender used for gender analyses. |
| `data/raw/mf_villages/mf_ind_covariates_wave2.dta` | Individual | Wave 2 individual IDs, household link, village ID, and gender used for gender analyses. |

## RCT Inputs

| Path | Unit | Description |
|---|---|---|
| `data/raw/rct_villages/rct_network_adjacency_layers.mat` | Village network matrices | RCT household network matrices. |
| `data/raw/rct_villages/rct_hh_vertex_crosswalk.mat` | Vertex-household key | Maps RCT network vertices to household IDs. |
| `data/raw/rct_villages/hh_covariates_by_village/covariates_*.csv` | Household | Village-level household covariates with caste, subcaste, occupation, roof type, and retained caste/occupation labels. |
| `data/raw/rct_villages/rct_hh_covariates.dta` | Household | Household-level RCT covariates with village ID, network vertex, asset measures, leader status, and seed indicator. |
| `data/raw/rct_villages/rct_village_diffusion_outcomes.dta` | Village | Village-level RCT diffusion outcomes, including call totals, randomized-household counts, and seed counts. |

## Simulation Inputs

| Path | Unit | Description |
|---|---|---|
| `data/raw/sims/sim_simple_contagion_village_*.npy` | Simulation grid by village | Simple-contagion simulation arrays. Columns are read as `q`, `delta`, `mpex1`, `inf1`, `mpex2`, `inf2`. |
| `data/raw/sims/sim_complex_contagion_village_*.npy` | Simulation grid by village | Complex-contagion simulation arrays with the same column order as the simple-contagion files. |

## Raw Adjacency Matrix Layer Order

The raw MATLAB adjacency-matrix files store village network layers as ordered arrays. The replication code reads the following layer order. The MF layer descriptions correspond to the paper's survey categories; the RCT village layers are the corresponding subset of these categories.

For `data/raw/mf_villages/mf_network_adjacency_layers.mat`:

| Position | Raw layer | Description |
|---:|---|---|
| 1 | `visitgo` | Respondent goes to another household's home. |
| 2 | `visitcome` | Another household comes to the respondent's home. |
| 3 | `nonrel` | Close non-relative social tie. |
| 4 | `rel` | Close relative living outside the household. |
| 5 | `medic` | Help in a medical emergency. |
| 6 | `keroricego` | Respondent borrows kerosene/rice. |
| 7 | `keroricecome` | Respondent lends kerosene/rice. |
| 8 | `templecompany` | Companion for temple, church, or mosque. |
| 9 | `bormoney` | Respondent borrows short-term money. |
| 10 | `lendmoney` | Respondent lends short-term money. |
| 11 | `decision` | Whom respondent turns to for help with important decisions. |
| 12 | not used | Not used by the replication code. |
| 13 | `advice` | Respondent gives information/advice. |

For `data/raw/rct_villages/rct_network_adjacency_layers.mat`:

1. `keroricego`
2. `keroricecome`
3. `visitcome`
4. `visitgo`
5. `decision`
6. `advice`
