library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, ggplot2, readxl, 
       lm.beta, readr, stringr, simDAG, data.table)

source("src/simCog.R")
source("src/cocalibrate.R")
dem_logr <- readRDS("models/dem_logr.Rds")

# Simulation Parameters -----------------------------------------------

n_sample <- 1000 # per group
prop_high_edu <- .30

b0 <- 0 # intercept of memory theta regressed on education
b1 <- .2 # slope of memory theta regressed on education

RecodedItemName_M <- read_csv("data/RecodedItemName_m.csv")
useMemItems <- which(!str_detect(RecodedItemName_M$RecodedItemName, "^rav.+b$"))
ItemNames <- RecodedItemName_M %>%
  slice(useMemItems) %>%
  pull(RecodedItemName)


# Group 1 -----------------------------------------------------------------

Loadings_M_G1 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                              sheet = "b_m") %>%
  select(a)

Thresholds_M_G1 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                              sheet = "b_m") %>%
  select(-a, -RecodedItemName)

item_pars_M_G1 <- traditional2mirt(bind_cols(Loadings_M_G1, Thresholds_M_G1), 
                                   "graded", 
                                   ncat = ncol(Thresholds_M_G1) + 1) 

set.seed(3987)
mem_fscores_G1 <- simCog(itempars = item_pars_M_G1, 
                         itemnames = ItemNames, 
                         edu_prob_1 = .3, 
                         b0 = 0, 
                         b1 = .2, 
                         n_sample = n_sample, 
                         cogname = "Mem", 
                         groupname = "Group 1")


# Group 2 -----------------------------------------------------------------

Loadings_M_G2 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                            sheet = "b_m_group2") %>%
  select(a)

Thresholds_M_G2 <- read_excel("~/Dropbox/PsyMCA2025_Sharing/W4.5.4/Memory/b_m.xlsx", 
                              sheet = "b_m_group2") %>%
  select(ends_with("_dif"), -DIF, -a_dif)

item_pars_M_G2 <- traditional2mirt(bind_cols(Loadings_M_G2, Thresholds_M_G2), 
                                   "graded", 
                                   ncat = ncol(Thresholds_M_G2) + 1)

set.seed(4365)
mem_fscores_G2 <- simCog(itempars = item_pars_M_G2, 
                         itemnames = ItemNames, 
                         edu_prob_1 = .3, 
                         b0 = 0, 
                         b1 = .2, 
                         n_sample = n_sample, 
                         cogname = "Mem", 
                         groupname = "Group 2")


# Combine -----------------------------------------------------------------

head(mem_fscores_G1)
head(mem_fscores_G2)

mem_fscores <- bind_rows(mem_fscores_G1, mem_fscores_G2)

mem_fscores <- mem_fscores %>%
  mutate(DemProb = predict(dem_logr, newdata = ., type = "response"),
         Dementia = ifelse(DemProb >= .9, 1, 0))

psych::describe(mem_fscores %>%
                  select(theta, edu, Mem_FS, DemProb, Dementia))

psych::describeBy(mem_fscores %>%
                    select(theta, edu, Mem_FS, DemProb, Dementia),
                  mem_fscores$Group)


# Define Scenarios -----------------------------------------------------------

scenario_map <- fread("data/scenarios.csv")
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
  items_group1 <- item_lists[[combos[s == scenario_num & g == 1, id]]]
  items_group2 <- item_lists[[combos[s == scenario_num & g == 2, id]]]
  items <- unique(c(items_group1, items_group2))
  
  ## get subset of variables
  subset <- subset[, c("theta", "edu", "Mem_FS", "Mem_FS_SE", "Group", "Dementia", items), with = FALSE]
  
  ## set missingness for variables only in one group
  subset[Group == "Group 1", setdiff(items, items_group1) := NA]
  subset[Group == "Group 2", setdiff(items, items_group2) := NA]
  
  return(subset)
}

mem_fscores <- as.data.table(mem_fscores)
scenario_data <- lapply(1:combos[, max(s)], get_scenario_data)



# Scenario 1 --------------------------------------------------------------

scenario_data1 <- scenario_data[[1]]

scenario_data1_rg <- scenario_data1 %>%
  filter(Group == "Group 1") %>%
  data.frame()

scenario_data1_fg <- scenario_data1 %>%
  filter(Group == "Group 2") %>%
  data.frame()

cc1 <- cocalibrate(rg_dat = scenario_data1_rg,
                   fg_dat = scenario_data1_fg,
                   rg_items = item_lists[[combos[s == 1 & g == 1, id]]],
                   fg_items = item_lists[[combos[s == 1 & g == 2, id]]])

cc1_data <- bind_rows(cc1$fscores_rg %>%
                        set_names(c("S1_Mem_FS", "S1_Mem_SE", "Group")), 
                      cc1$fscores_fg %>%
                        set_names(c("S1_Mem_FS", "S1_Mem_SE", "Group")))

psych::describeBy(cc1_data, cc1_data$Group)


mem_fscores <- bind_cols(mem_fscores, cc1_data %>% select(-Group))

psych::describeBy(mem_fscores, mem_fscores$Group)



# Scenario 2 --------------------------------------------------------------

scenario_data2 <- scenario_data[[2]]

scenario_data2_rg <- scenario_data2 %>%
  filter(Group == "Group 1") %>%
  data.frame()

scenario_data2_fg <- scenario_data2 %>%
  filter(Group == "Group 2") %>%
  data.frame()

cc2 <- cocalibrate(rg_dat = scenario_data2_rg,
                   fg_dat = scenario_data2_fg,
                   rg_items = item_lists[[combos[s == 2 & g == 1, id]]],
                   fg_items = item_lists[[combos[s == 2 & g == 2, id]]])

cc2_data <- bind_rows(cc2$fscores_rg %>%
                        set_names(c("S2_Mem_FS", "S2_Mem_SE", "Group")), 
                      cc2$fscores_fg %>%
                        set_names(c("S2_Mem_FS", "S2_Mem_SE", "Group")))

psych::describeBy(cc2_data, cc2_data$Group)


mem_fscores <- bind_cols(mem_fscores, cc2_data %>% select(-Group))

psych::describeBy(mem_fscores, mem_fscores$Group)


# Scenario 3 --------------------------------------------------------------

scenario_data3 <- scenario_data[[3]]

scenario_data3_rg <- scenario_data3 %>%
  filter(Group == "Group 1") %>%
  data.frame()

scenario_data3_fg <- scenario_data3 %>%
  filter(Group == "Group 2") %>%
  data.frame()

cc3 <- cocalibrate(rg_dat = scenario_data3_rg,
                   fg_dat = scenario_data3_fg,
                   rg_items = item_lists[[combos[s == 3 & g == 1, id]]],
                   fg_items = item_lists[[combos[s == 3 & g == 2, id]]])

cc3_data <- bind_rows(cc3$fscores_rg %>%
                        set_names(c("S3_Mem_FS", "S3_Mem_SE", "Group")), 
                      cc3$fscores_fg %>%
                        set_names(c("S3_Mem_FS", "S3_Mem_SE", "Group")))

psych::describeBy(cc3_data, cc3_data$Group)


mem_fscores <- bind_cols(mem_fscores, cc3_data %>% select(-Group))

psych::describeBy(mem_fscores, mem_fscores$Group)


# Scenario 4 --------------------------------------------------------------

scenario_data4 <- scenario_data[[4]]

scenario_data4_rg <- scenario_data4 %>%
  filter(Group == "Group 1") %>%
  data.frame()

scenario_data4_fg <- scenario_data4 %>%
  filter(Group == "Group 2") %>%
  data.frame()

cc4 <- cocalibrate(rg_dat = scenario_data4_rg,
                   fg_dat = scenario_data4_fg,
                   rg_items = item_lists[[combos[s == 4 & g == 1, id]]],
                   fg_items = item_lists[[combos[s == 4 & g == 2, id]]])

cc4_data <- bind_rows(cc4$fscores_rg %>%
                        set_names(c("S4_Mem_FS", "S4_Mem_SE", "Group")), 
                      cc4$fscores_fg %>%
                        set_names(c("S4_Mem_FS", "S4_Mem_SE", "Group")))

psych::describeBy(cc4_data, cc4_data$Group)


mem_fscores <- bind_cols(mem_fscores, cc4_data %>% select(-Group))

psych::describeBy(mem_fscores, mem_fscores$Group)



