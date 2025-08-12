library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2, readxl, 
       lm.beta, readr, stringr, simDAG)

code_dir <- "C:/Users/emmanich/code/crosswalk-and-co/"
dropbox_dir <- "C:/Users/emmanich/P2AGING Dropbox/Emma Nichols/"

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

RecodedItemName_M <- read_csv(paste0(code_dir, "data/RecodedItemName_m.csv"))
useMemItems <- which(!str_detect(RecodedItemName_M$RecodedItemName, "^rav.+b$"))
ItemNames <- RecodedItemName_M %>%
  slice(useMemItems) %>%
  pull(RecodedItemName)

Loadings_M <- read_csv(paste0(code_dir, "data/a_m.csv")) %>%
  slice(useMemItems)

# Group 1 -----------------------------------------------------------------

Thresholds_M_G1 <- read_excel(paste0(dropbox_dir, "PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx"), 
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

Thresholds_M_G2 <- read_excel(paste0(dropbox_dir, "PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx"), 
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
                         groupname = "Group 2")


# Combine -----------------------------------------------------------------

mem_fscores <- as.data.table(bind_rows(mem_fscores_G1, mem_fscores_G2))

# Define Scenarios -----------------------------------------------------------

scenario_map <- fread(paste0(code_dir, "data/scenarios.csv"))
items <- names(scenario_map)[!names(scenario_map) %in% c("scenario", "outcome")]
scenario_map <- melt.data.table(scenario_map, id.vars = c("scenario", "outcome"), variable.name = "item", value.name = "value")

items <- function(s, o){
  as.character(scenario_map[scenario == s & outcome == o & value == 1, item])
}

combos <- as.data.table(expand.grid(s = 1:4, o = 1:2))
item_lists <- purrr::pmap(combos, items)

## make id variable for combos of scenarios and groups
combos[, id := 1:.N]

# Get scenario data -----------------------------------------------------------

get_scenario_data <- function(scenario_num, data = mem_fscores){
  subset <- copy(data)

  ## get items
  items_group1 <- item_lists[[combos[s == scenario_num & o == 1, id]]]
  items_group2 <- item_lists[[combos[s == scenario_num & o == 2, id]]]
  items <- unique(c(items_group1, items_group2))

  ## get subset of variables
  subset <- subset[, c("theta", "edu", "Mem_FS", "Mem_FS_SE", "Group", items), with = FALSE]

  ## set missingness for variables only in one group
  subset[Group == "Group 1", setdiff(items, items_group1) := NA]
  subset[Group == "Group 2", setdiff(items, items_group2) := NA]

  return(subset)
}

scenario_data <- lapply(1:combos[, max(s)], get_scenario_data)
