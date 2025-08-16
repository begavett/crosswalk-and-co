# PsyMCA 2025 ---------------------------------------------
# Project       : Crosswalk and Co Simulation Study
# Script        : sim_cwxco.R
# Authors       : Brandon Gavett, Emma Nichols, Nathan Ryder
# Collaborators : Hoda Abdel Magid, Sarah Forrester, Xiaqing Jiang, Katrina Kezios, Adina Zeki Al Hazzouri 
# Created       : 2025-08-12
# Modified      : 2025-08-16
# Version       : 1.0
#
# Description:
#   This script simulates cognitive item response data under different DIF and 
#   anchor-strength scenarios, cocalibrates across groups, estimates regression 
#   relationships with education, and compares results to crosswalks generated 
#   with the `cogxwalkr` package. Parallel processing is used to replicate 
#   simulations efficiently.
#
# Dependencies:
#   - R (>= 4.3)
#   - Packages: dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2,
#               readxl, readr, stringr, simDAG, cogxwalkr, forcats,
#               future, future.apply, progressr
#   - Local source files: helper.R
#   - External objects: models/dem_logr.Rds, data/b_m.xlsx, data/scenarios.csv,
#                       data/RecodedItemName_m.csv
#
# Inputs:
#   - Simulation parameters (sample size, education proportion, regression betas).
#   - Item parameter files for Groups 1 and 2.
#   - Scenario map defining DIF and anchor conditions.
#
# Outputs:
#   - Data frames of simulated factor scores with cocalibration and crosswalks.
#   - Model objects saved as RDS files (cc_models.Rds).
#   - Combined results for comparison across methods.
#
# Notes:
#   - Parallel execution via {future} requires setting a plan and 
#     may need adjustments on cluster/HPC environments.
#   - Simulation seed is optional but recommended for reproducibility.
#   - Memory use may be high for large n_sample × niter.
# ******************************************************************************

# remotes::install_github("jrgant/cogxwalkr")
library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2, readxl, 
       readr, stringr, simDAG, data.table, cogxwalkr, forcats,
       future, future.apply, progressr)

user <- Sys.info()[["user"]]

if (user == "brandon") {
  code_dir <- "~/Dropbox/Projects/crosswalk-and-co/"
} else if (user == "emmanich") {
  code_dir <- "C:/Users/emmanich/code/crosswalk-and-co/"
}

source(paste0(code_dir, "src/helper.R"))
dem_logr <- readRDS(paste0(code_dir, "models/dem_logr.Rds"))

# Setup -----------------------------------------------

scenario_labels <- data.table(scenario = 1:4, 
                              slabel = c("No DIF + Strong anchor", "No DIF + Weak anchor", 
                                         "Weak DIF + Strong anchor", "Strong DIF + Strong anchor"))

RecodedItemName_M <- read_csv(paste0(code_dir, "data/RecodedItemName_m.csv"))
useMemItems <- which(!str_detect(RecodedItemName_M$RecodedItemName, "^rav.+b$"))
ItemNames <- RecodedItemName_M %>%
  slice(useMemItems) %>%
  pull(RecodedItemName)


# Run in Parallel ---------------------------------------------------------

iter_start <- 1
niter <- 50
iter_end <- iter_start + niter - 1
# set.seed(33333)
# results <- replicate(niter, sim_cwxco(n_sample = 5002, prop_high_edu = .30, b = 0, b1 = .2), 
#                      simplify = FALSE)


plan(multisession, workers = 8)

# Enable progress reporting
handlers(global = TRUE)
options(future.globals.maxSize = 1500*1024^2)

# run the simulations in parallel
with_progress({
  p <- progressor(along = iter_start:iter_end)
  
  results <- future_Map(function(i) {
    sim_cwxco(n_sample = 5002, prop_high_edu = .30, b = 0, b1 = .2)
    p()  # Report progress
  }, 
  i = iter_start:iter_end,
  future.seed = TRUE)
})

# shut down parallelization
plan(sequential)
closeAllConnections()

# Cocalibrated models -----------------------------------------------------

if(file.exists("models/cc_models.Rds")) {
  cc_models <- readRDS("models/cc_models.Rds")
  
  cc_models$cc1$mirt_rg %>% plot(type = "infoSE", main = "Scenario 1: No DIF + Strong Anchors")
  cc_models$cc2$mirt_rg %>% plot(type = "infoSE", main = "Scenario 2: No DIF + Weak Anchors")
  cc_models$cc3$mirt_rg %>% plot(type = "infoSE", main = "Scenario 3: Weak DIF + Strong Anchors")
  cc_models$cc4$mirt_rg %>% plot(type = "infoSE", main = "Scenario 4: Strong DIF + Strong Anchors")
}

# source("summarize-results.R)
