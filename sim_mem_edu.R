library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2, readxl, 
       lm.beta, readr, stringr, simDAG)

source("src/simCog.R")

# Simulation Parameters -----------------------------------------------

n_sample <- 500 # per group
prop_high_edu <- .30
var_edu <- prop_high_edu * (1 - prop_high_edu)

b0 <- 0 # intercept of memory theta regressed on education
b1 <- .2 # slope of memory theta regressed on education

# Error variance to make Var(y) = 1
sigma2 <- 1 - (b1^2 * var_edu)
sigma <- sqrt(sigma2)


## Memory Parameters --------------------------------------------------

RecodedItemName_M <- read_csv("data/RecodedItemName_m.csv")
useMemItems <- which(!str_detect(RecodedItemName_M$RecodedItemName, "^rav.+b$"))
ItemNames <- RecodedItemName_M %>%
  slice(useMemItems) %>%
  pull(RecodedItemName)

Loadings_M <- read_csv("data/a_m.csv") %>%
  slice(useMemItems)


# Group 1 -----------------------------------------------------------------

Thresholds_M_G1 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                              sheet = "b_m") %>%
  select(-a, -RecodedItemName)

item_pars_M_G1 <- traditional2mirt(bind_cols(Loadings_M, Thresholds_M_G1), 
                                   "graded", 
                                   ncat = ncol(Thresholds_M_G1) + 1) 

set.seed(3987)
mem_fscores_G1 <- simCog(itempars = item_pars_M_G1, 
                         itemnames = ItemNames, 
                         edu_prob_1 = .3, 
                         b0 = 0, 
                         b1 = .2, 
                         n_sample = 500, 
                         cogname = "Mem", 
                         groupname = "Group 1")


# Group 2 -----------------------------------------------------------------

Thresholds_M_G2 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                              sheet = "b_m_group2") %>%
  select(ends_with("_dif"), -DIF, -a_dif)

item_pars_M_G2 <- traditional2mirt(bind_cols(Loadings_M, Thresholds_M_G2), 
                                   "graded", 
                                   ncat = ncol(Thresholds_M_G2) + 1)

set.seed(4365)
mem_fscores_G2 <- simCog(itempars = item_pars_M_G2, 
                         itemnames = ItemNames, 
                         edu_prob_1 = .3, 
                         b0 = 0, 
                         b1 = .2, 
                         n_sample = 500, 
                         cogname = "Mem", 
                         groupname = "Group 1")


# Combine -----------------------------------------------------------------

head(mem_fscores_G1)
head(mem_fscores_G2)

mem_fscores <- bind_rows(mem_fscores_G1, mem_fscores_G2)


