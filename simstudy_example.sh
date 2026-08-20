#!/bin/bash
#
#SBATCH -J simstudy
#SBATCH --array=1-1000 # performs job in a loop of defined length 

date
echo ""

cd /home/homedir


## start r process not saving workspace but saving output
R --no-save > /home/homedir/simstudy$SLURM_ARRAY_TASK_ID.out 2>&1 <<EOF

library(tidyverse)
library(readxl)
library(simDAG)
library(mirt)
library(data.table)
library(cogxwalkr)


n_sample <- 5000 # SAMPLE SIZE per group
b1_val <- 0.4  # Effect Size
seedpad <- 000
iter <- $SLURM_ARRAY_TASK_ID + seedpad # Seed
CI <- TRUE  # bootstrap CI's (500 samples) for cogxwalkr method

prop_high_edu <- 0.3
b0 <- 0
dem_prop_cut <- 0.9


code_dir <- "/home/homedir/crosswalk-and-co-main/"
setwd(code_dir)
outdir <- paste0("/home/homedir/results/withci_n", n_sample, "_d", b1_val*10, "/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dem_logr <- readRDS(paste0(code_dir, "models/dem_logr.Rds"))
dem_logr_noedu <- readRDS(paste0(code_dir, "models/dem_logr_noedu.Rds"))
source(paste0(code_dir, "src/helper_sum.R"))  ######################### different for sum


scenario_labels <- data.table(scenario = 1:4, 
                              slabel = c("No DIF + Strong anchor", "No DIF + Weak anchor", 
                                         "Weak DIF + Strong anchor", "Strong DIF + Strong anchor"))



results <- sim_cwxco(iter = iter, N = n_sample, prop_high_edu = prop_high_edu, b0 = b0, b1 = b1_val, 
                     dem_prob_cut = dem_prop_cut, CI = CI)

saveRDS(results, paste0(outdir, "result", $SLURM_ARRAY_TASK_ID, ".RDS"))

EOF

echo "Finished"

