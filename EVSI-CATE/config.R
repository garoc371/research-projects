# config.R
library(tidyverse)
library(data.table)
library(dbarts)
library(rstanarm)
library(lhs)
library(future.apply)
library(randtoolbox)
library(MCMCpack)
library(future)


# Simulation parameters
N_bhm <- 600        # Total sample size for trial
n_target <- 20000
age_grid <- 40:80
extended_age_grid <- 40:120

# Utility values
utility_healthy  <- 1
utility_diseased <- 0.7
utility_dead <- 0

# Cost values
cost_healthy <- 100
cost_diseased <- 500
cost_dead <- 0
treatment_cost_per_cycle <- 5000

# Time horizon and number of cycles
time_horizon <- 40

# Simulation counts
n_simulations <- 4000
n_patients    <- n_target
n_cycles      <- time_horizon

# Baseline probability and subgroup settings
baseline_probs <- 0.4

# Directories for results and logs
results_dir <- "sim/results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# Directory for storing generated data
data_dir <- "sim/data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
