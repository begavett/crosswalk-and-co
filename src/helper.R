# cocalibrate -------------------------------------------------------------

#' Cocalibrate item response models across two groups
#'
#' This function performs scale cocalibration between a reference group (rg) 
#' and a focal group (fg) using the `mirt` package. Linking items common to both 
#' groups are used to anchor the scale. The procedure first estimates item 
#' parameters in the reference group, then fixes the linking item parameters in 
#' the focal group to align the scales, while allowing group-level mean and 
#' variance to differ.
#'
#' @param rg_dat Data frame containing item responses for the reference group.
#' @param fg_dat Data frame containing item responses for the focal group.
#' @param rg_items Character vector of item names for the reference group model.
#' @param fg_items Character vector of item names for the focal group model.
#'
#' @return A list containing:
#'   \item{mirt_rg}{Fitted mirt model for the reference group.}
#'   \item{mirt_fg}{Fitted mirt model for the focal group with linking constraints.}
#'   \item{fscores_rg}{Factor scores (and SEs) for the reference group.}
#'   \item{fscores_fg}{Factor scores (and SEs) for the focal group.}
#'   \item{linking_items}{Character vector of items used for linking.}
#'
#' @details
#' The function:
#'   1. Identifies common items across groups to serve as linking items.
#'   2. Estimates a unidimensional graded response model for the reference group.
#'   3. Initializes the focal group model with free parameters, then fixes 
#'      linking item parameters to the reference group estimates.
#'   4. Re-estimates the focal group model with group mean and variance free.
#'   5. Returns factor scores for both groups on the same scale.
#'
#' @examples
#' # Example usage (assuming data prepared):
#' # result <- cocalibrate(rg_dat, fg_dat, rg_items, fg_items)
#' # head(result$fscores_rg)
#' # head(result$fscores_fg)
#'

cocalibrate <- function(rg_dat, 
                        fg_dat, 
                        rg_items,
                        fg_items) {
  
  library(mirt)
  library(dplyr)
  library(magrittr)
  
  linking_items <- intersect(rg_items, fg_items)
  
  if(length(linking_items) == 0) stop ("No linking items")
  
  mirt_rg <- mirt(data = rg_dat %>%
                    select(all_of(rg_items)),
                  model = 1,
                  itemtype = "graded",
                  technical = list(NCYCLES = 5000))
  
  itembank_rg <- data.frame(coef(mirt_rg, simplify = TRUE)$items) %>%
    mutate(item = rownames(.))
  
  linking_itembank_rg <- itembank_rg %>%
    filter(item %in% linking_items)
  
  # multipleGroup estimates parameters using all data,
  # which I don't think is right
  # I think we want to estimate parameters in Group 1 first
  # Otherwise the estimated paramters in a MG model will be based on
  # both G1 and G2 data
  # mirt_mg <- multipleGroup(data = bind_rows(rg_dat, fg_dat) %>%
  #                 select(all_of(union(rg_items, fg_items))),
  #               model = 1,
  #               group = c(rg_dat$Group, fg_dat$Group),
  #               itemtype = "graded",
  #               invariance = c(linking_items, "free_means", "free_vars"), SE = TRUE)
  # 
  # coef(mirt_mg, simplify = TRUE)
  # fscores(mirt_mg, full.scores.SE = TRUE)
  
  
  mirt_fg0 <- mirt(data = fg_dat %>%
                     select(all_of(fg_items)),
                   model = 1,
                   itemtype = "graded",
                   pars = "values",
                   technical = list(NCYCLES = 5000))
  
  mirt_fg0_partable <- mirt_fg0
  
  for(i in linking_items) {
    
    i_pars <- linking_itembank_rg %>%
      t() %>%
      data.frame() %>%
      mutate(par = rownames(.)) %>%
      select(all_of(i), par) %>%
      filter(!is.na(!!sym(i))) %>%
      slice(1:(n()-1)) %>%
      pull(par)
    
    for(p in i_pars) {
      mirt_fg0_partable$value[mirt_fg0_partable$item == i & mirt_fg0_partable$name == p] <- linking_itembank_rg[i, p]
      mirt_fg0_partable$est[mirt_fg0_partable$item == i & mirt_fg0_partable$name == p] <- FALSE
    }
    
    mirt_fg0_partable$est[mirt_fg0_partable$name == "MEAN_1"] <- TRUE
    mirt_fg0_partable$est[mirt_fg0_partable$name == "COV_11"] <- TRUE
  }
  
  mirt_fg <- mirt(data = fg_dat %>%
                    select(all_of(fg_items)),
                  model = 1,
                  itemtype = "graded",
                  pars = mirt_fg0_partable,
                  technical = list(NCYCLES = 5000))
  
  coef(mirt_rg, simplify = TRUE)
  coef(mirt_fg, simplify = TRUE)
  
  fscores_rg <- data.frame(fscores(mirt_rg, full.scores.SE = TRUE)) %>%
    mutate(Group = "Group 1")
  fscores_fg <- data.frame(fscores(mirt_fg, full.scores.SE = TRUE)) %>%
    mutate(Group = "Group 2")
  
  return(list(mirt_rg = mirt_rg, 
              mirt_fg = mirt_fg, 
              fscores_rg = fscores_rg, 
              fscores_fg = fscores_fg,
              linking_items = linking_items))
}


# cocalibrated_regressions ------------------------------------------------
#
# Description:
#   Fits parallel linear regressions in two groups to estimate how education 
#   (edu) predicts memory factor scores for a given scenario, and extracts 
#   the group-specific regression coefficients as a basis for crosswalks.
#
# Arguments:
#   scenario - A numeric or character identifier for the scenario of interest. 
#              Used to dynamically construct the outcome variable name in the 
#              form "S[scenario]_Mem_FS".
#
# Process:
#   1. Builds a regression formula: "S[scenario]_Mem_FS ~ edu".
#   2. Fits the model separately in Group 1 and Group 2 using the dataset 
#      'mem_fscores'.
#   3. Extracts the coefficient for 'edu' in each group using the 
#      'parameters' package and converts results to a data.table.
#   4. Creates a crosswalk-style table showing, for each scenario, the 
#      coefficient linking education to memory scores in each group.
#   5. Merges the results with a 'scenario_labels' lookup table for clarity.
#
# Returns:
#   A data.table with the following columns:
#     - scenario     : Scenario identifier.
#     - crosswalk_to : The target group for crosswalk comparison.
#     - coef         : Regression coefficient for 'edu' in the given group.
#     - (plus any additional labels from 'scenario_labels').
#
# Notes:
#   - Assumes 'mem_fscores' contains the variables:
#       * S[scenario]_Mem_FS (memory factor score for the scenario)
#       * edu (education variable)
#       * Group (with at least "Group 1" and "Group 2")
#   - Relies on 'parameters::model_parameters()' for standardized extraction 
#     of model results.

cocalibrated_regressions <- function(scenario){
  library(parameters)
  library(data.table)
  
  model_formula <- as.formula(paste0("S", scenario, "_Mem_FS ~ edu"))
  
  g1_model <- lm(model_formula, data = mem_fscores[Group == "Group 1"])
  g1_params <- as.data.table(parameters::model_parameters(g1_model))
  
  g2_model <- lm(model_formula, data = mem_fscores[Group == "Group 2"])
  g2_params <- as.data.table(parameters::model_parameters(g2_model))
  
  result <- data.table(scenario = scenario, 
                       crosswalk_to = c("Group 2", "Group 1"), 
                       coef = c(g1_params$Coefficient[g1_params$Parameter == "edu"], 
                                g2_params$Coefficient[g2_params$Parameter == "edu"]))
  result <- merge(result, scenario_labels, by = "scenario")                              
  
  return(result)                               
}


# getcrosswalk ------------------------------------------------------------
#
# Description:
#   Creates a score crosswalk between two cognitive measures, optionally
#   conditioned on dementia status, to enable score comparability across tests.
#
# Arguments:
#   groupdata - A data frame containing the cognitive measures of interest 
#               and a dementia indicator variable.
#   cog1      - A character string giving the column name of the first 
#               cognitive measure.
#   cog2      - A character string giving the column name of the second 
#               cognitive measure.
#
# Process:
#   1. Subsets 'groupdata' to include only the two measures and the Dementia variable.
#   2. Runs the crosswalk procedure using the 'crosswalk' function, 
#      conditioning on Dementia status and looping across levels if specified.
#   3. (Optional) Control parameters such as bootstrap replications, 
#      random seed, and parallel cores can be set in 'control'.
#
# Returns:
#   An object from the 'crosswalk' function containing the estimated 
#   crosswalk between the two measures.
#
# Notes:
#   - The Dementia variable is assumed to exist in 'groupdata' and is 
#     used as a conditioning variable.
#   - The 'crosswalk' function handles the statistical mapping between 
#     the two cognitive measures.

getcrosswalk <- function(groupdata, cog1, cog2)
{
  library(cogxwalkr)
  library(dplyr)
  library(magrittr)
  groupdata <- groupdata %>% select(all_of(c(cog1, cog2)), Dementia)
  crswkout <- crosswalk(cog1, cog2, condition_by = "Dementia",
                        condition_loop = TRUE, data = groupdata)
  ##control = list(nboot = 100, seed = 42, ncores = 4L)
}


# get_scenario_data -------------------------------------------------------

#' Extract scenario-specific item response data
#'
#' This function subsets a dataset of memory factor scores and item responses
#' to return only the variables relevant to a given simulation scenario. 
#' Each scenario corresponds to a predefined combination of items for 
#' Group 1 and Group 2.
#'
#' @param scenario_num Integer. Scenario identifier used to select the 
#'   appropriate item set from `combos` and `item_lists`.
#' @param data Data table containing factor scores, group labels, and item 
#'   responses. Defaults to `mem_fscores`.
#'
#' @return A data.table containing:
#'   \item{theta}{True latent ability values.}
#'   \item{edu}{Education variable.}
#'   \item{Mem_FS}{Estimated memory factor score.}
#'   \item{Mem_FS_SE}{Standard error of memory factor score.}
#'   \item{Group}{Group identifier (e.g., "Group 1" / "Group 2").}
#'   \item{Dementia}{Dementia status variable.}
#'   \item{[items]}{All items relevant to the specified scenario. Items not 
#'                  belonging to a given group can optionally be set to `NA`.}
#'
#' @details
#' The function:
#'   1. Looks up the set of items for Group 1 and Group 2 associated with the 
#'      given `scenario_num`.
#'   2. Combines the group-specific items into a single unique item list.
#'   3. Subsets the dataset to include only key variables and scenario items.
#'   4. (Optional) Sets missing values for items not administered to a group 
#'      (currently commented out).
#'
#' @examples
#' \dontrun{
#' # Extract data for scenario 5
#' scenario_data <- get_scenario_data(5, data = mem_fscores)
#' head(scenario_data)
#' }
#'

get_scenario_data <- function(scenario_num, data = mem_fscores){
  subset <- copy(data)
  
  ## get items
  items_group1 <- item_lists[[combos[s == scenario_num & g == 1, id]]]
  items_group2 <- item_lists[[combos[s == scenario_num & g == 2, id]]]
  items <- unique(c(items_group1, items_group2))
  
  ## get subset of variables
  subset <- subset[, c("theta", "edu", "Mem_FS", "Mem_FS_SE", "Group", "Dementia", items), with = FALSE]
  
  ## set missingness for variables only in one group
  # subset[Group == "Group 1", setdiff(items, items_group1) := NA]
  # subset[Group == "Group 2", setdiff(items, items_group2) := NA]
  
  return(subset)
}

# getscore ----------------------------------------------------------------
#
# Description:
#   Fits a unidimensional graded response IRT model to a set of items and 
#   appends the resulting factor scores to the original dataset.
#
# Arguments:
#   scendat  - A data frame containing the item responses (and any other variables).
#   itemnms  - A character vector of column names in 'scendat' that represent 
#              the item responses to be modeled.
#   scorenm  - A character string specifying the name to assign to the 
#              newly created factor score column.
#
# Process:
#   1. Selects the item columns from the dataset.
#   2. Fits a 1-factor graded response model using the mirt package.
#   3. Extracts factor scores (no standard errors).
#   4. Renames the score column as specified in 'scorenm'.
#   5. Binds the score back to the original dataset.
#
# Returns:
#   A data frame identical to 'scendat' but with an added column containing 
#   the estimated factor scores.
#
# Notes:
#   - By default, the model is unidimensional with graded response items.
#   - The function uses 5000 EM cycles to help ensure convergence.

getscore <- function(scendat, itemnms, scorenm)
{
  library(dplyr)
  library(magrittr)
  library(mirt)
  data <- scendat %>% select(all_of(itemnms)) 
  mirtout <- mirt(data, model = 1, itemtype = "graded",
                  technical = list(NCYCLES = 5000))
  score <- fscores(mirtout, full.scores.SE = FALSE)
  colnames(score) <- scorenm
  scendat <- cbind.data.frame(scendat, score)
  ## left_join(cbind.data.frame(score, ID = data$ID), by = "ID")
  return(scendat)
}


# sim_cwxco ---------------------------------------------------------------

#' Simulate and compare cocalibration vs. crosswalk approaches
#'
#' This function generates item response data for two groups, 
#' applies cocalibration under multiple DIF/anchor scenarios, 
#' and compares regression estimates with those obtained via crosswalks 
#' using the `cogxwalkr` package.
#'
#' @param n_sample Integer. Number of participants per group (default = 5002).
#' @param prop_high_edu Numeric in [0,1]. Proportion of participants 
#'   with high education in the simulation (default = 0.30).
#' @param b Numeric. Intercept parameter for the latent regression 
#'   on education (default = 0).
#' @param b1 Numeric. Slope parameter for the latent regression 
#'   on education (default = 0.20).
#'
#' @details
#' The function:
#'   1. Reads group-specific item parameters from Excel files.
#'   2. Simulates cognitive scores for Group 1 and Group 2 using `simCog()`.
#'   3. Defines DIF/anchor scenarios based on an external `scenarios.csv`.
#'   4. Runs cocalibration with `cocalibrate()` across 4 scenarios.
#'   5. Fits regression models on cocalibrated scores (`cocalibrated_regressions()`).
#'   6. Computes crosswalks for each scenario (`getcrosswalk()`).
#'   7. Combines results into a single data frame for comparison.
#'
#' @return A data frame with regression coefficients across scenarios, 
#'   methods (`Cocalibration` vs. `cogxwalkr`), and groups.
#'
#' @examples
#' \dontrun{
#' results <- sim_cwxco(n_sample = 1000, prop_high_edu = 0.25, b = 0, b1 = 0.3)
#' head(results)
#' }
#'
#' @note
#'   - Saves cocalibration model objects to `models/cc_models.Rds` if not present.
#'   - Requires external files in `data/` (item parameters, scenario map).
#'   - Memory usage may be high for large sample sizes.
#'
#' @seealso [simCog()], [cocalibrate()], [getcrosswalk()], [cocalibrated_regressions()]
#'
#'

sim_cwxco <- function(n_sample = 5002, prop_high_edu = .30, b = 0, b1 = .2) {
  
  
  # Group 1 -----------------------------------------------------------------
  
  Loadings_M_G1 <- read_excel(paste0(code_dir, "data/b_m.xlsx"), 
                              sheet = "b_m") %>%
    select(a)
  
  Thresholds_M_G1 <- read_excel(paste0(code_dir, "data/b_m.xlsx"), 
                                sheet = "b_m") %>%
    select(-a, -RecodedItemName)
  
  item_pars_M_G1 <- traditional2mirt(bind_cols(Loadings_M_G1, Thresholds_M_G1), 
                                     "graded", 
                                     ncat = ncol(Thresholds_M_G1) + 1) 
  
  
  mem_fscores_G1 <- simCog(itempars = item_pars_M_G1, 
                           itemnames = ItemNames, 
                           prop_high_edu = .3, 
                           b0 = 0, 
                           b1 = .2, 
                           n_sample = n_sample, 
                           cogname = "Mem", 
                           groupname = "Group 1")
  
  
  # Group 2 -----------------------------------------------------------------
  
  Loadings_M_G2 <- read_excel(paste0(code_dir, "data/b_m.xlsx"), 
                              sheet = "b_m_group2") %>%
    select(a = a_dif)
  
  Thresholds_M_G2 <- read_excel(paste0(code_dir, "data/b_m.xlsx"), 
                                sheet = "b_m_group2") %>%
    select(ends_with("_dif"), -DIF, -a_dif)
  
  item_pars_M_G2 <- traditional2mirt(bind_cols(Loadings_M_G2, Thresholds_M_G2), 
                                     "graded", 
                                     ncat = ncol(Thresholds_M_G2) + 1)
  
  mem_fscores_G2 <- simCog(itempars = item_pars_M_G2, 
                           itemnames = ItemNames, 
                           edu_prob_1 = .3, 
                           b0 = 0, 
                           b1 = .2, 
                           n_sample = n_sample, 
                           cogname = "Mem", 
                           groupname = "Group 2")
  
  
  # Combine -----------------------------------------------------------------
  
  # head(mem_fscores_G1)
  # head(mem_fscores_G2)
  
  mem_fscores <- bind_rows(mem_fscores_G1, mem_fscores_G2)
  
  mem_fscores <- mem_fscores %>%
    mutate(DemProb = predict(dem_logr, newdata = ., type = "response"),
           Dementia = ifelse(DemProb >= .9, 1, 0))
  
  # psych::describe(mem_fscores %>%
  #                   select(theta, edu, Mem_FS, DemProb, Dementia))
  # 
  # psych::describeBy(mem_fscores %>%
  #                     select(theta, edu, Mem_FS, DemProb, Dementia),
  #                   mem_fscores$Group)
  
  
  # Define Scenarios -----------------------------------------------------------
  
  scenario_map <- fread(paste0(code_dir, "data/scenarios.csv"))
  items <- names(scenario_map)[!names(scenario_map) %in% c("scenario", "outcome")]
  scenario_map <- melt.data.table(scenario_map, id.vars = c("scenario", "outcome"), 
                                  variable.name = "item", value.name = "value")
  
  items <- function(s, g){
    as.character(scenario_map[scenario == s & outcome == g & value == 1, item])
  }
  
  combos <- as.data.table(expand.grid(s = 1:4, g = 1:2))
  item_lists <- purrr::pmap(combos, items)
  
  ## make id variable for combos of scenarios and groups
  combos[, id := 1:.N]
  
  # Get scenario data -----------------------------------------------------------
  
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
  
  # psych::describeBy(cc1_data, cc1_data$Group)
  
  
  mem_fscores <- bind_cols(mem_fscores, cc1_data %>% select(-Group))
  
  # psych::describeBy(mem_fscores, mem_fscores$Group)
  
  
  
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
  
  # psych::describeBy(cc2_data, cc2_data$Group)
  
  
  mem_fscores <- bind_cols(mem_fscores, cc2_data %>% select(-Group))
  
  # psych::describeBy(mem_fscores, mem_fscores$Group)
  
  
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
  
  # psych::describeBy(cc3_data, cc3_data$Group)
  
  
  mem_fscores <- bind_cols(mem_fscores, cc3_data %>% select(-Group))
  
  # psych::describeBy(mem_fscores, mem_fscores$Group)
  
  
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
  
  # psych::describeBy(cc4_data, cc4_data$Group)
  
  
  mem_fscores <- bind_cols(mem_fscores, cc4_data %>% select(-Group))
  
  # psych::describeBy(mem_fscores, mem_fscores$Group)
  
  
  
  # Run regression models on cocalibrated scores -------------------------------------
  
  
  
  cocalibration_results <- rbindlist(lapply(1:combos[, max(s)], cocalibrated_regressions))
  
  
  # x <- ggplot(mem_fscores, aes(x = edu, y = S2_Mem_FS)) + 
  #   geom_jitter() + 
  #   theme_bw()
  # 
  # x
  
  # Crosswalk ---------------------------------------------------------------
  
  
  crswkoutlist <- outdfs <- datasets <- vector("list", 4)
  for(scen in 1:4)
  {
    itemnmslist <- list(item_lists[[combos[s == scen & g == 1, id]]],
                        item_lists[[combos[s == scen & g == 2, id]]])
    groupdatalist <- list(scenario_data[[scen]] %>% filter(Group == "Group 1") %>%
                            getscore(itemnmslist[[1]], "score1") %>%
                            getscore(itemnmslist[[2]], "score2"),
                          scenario_data[[scen]] %>% filter(Group == "Group 2") %>%
                            getscore(itemnmslist[[1]], "score1") %>%
                            getscore(itemnmslist[[2]], "score2"))
    datasets[[scen]] <- groupdatalist
    crswkoutlist[[scen]] <- list(getcrosswalk(groupdatalist[[1]], "score2", "score1"),
                                 getcrosswalk(groupdatalist[[2]], "score1", "score2"))
    edufits <- list(lm(score2 ~ edu, groupdatalist[[2]]),
                    lm(score1 ~ edu, groupdatalist[[1]]))
    coefs = sapply(1:2, function(i)
    {
      coef(crswkoutlist[[scen]][[i]]$fit)*coef(edufits[[i]])["edu"]
    })
    outdfs[[scen]] <- cbind.data.frame(scenario = scen, crosswalk_to = c("Group 2", "Group 1"),
                                       coefs = coefs)
  }
  allout <- do.call(rbind.data.frame, outdfs)
  
  cocalibration_results
  
  cwxco <- cocalibration_results %>%
    mutate(Method = "Cocalibration") %>%
    bind_rows(allout %>%
                mutate(slabel = cocalibration_results$slabel,
                       Method = "cogxwalkr") %>%
                rename(coef = coefs))
  
  cc_models <- list(cc1 = cc1, cc2 = cc2, cc3 = cc3, cc4 = cc4)
  if(!file.exists("models/cc_models.Rds")) saveRDS(cc_models, "models/cc_models.Rds")
  
  return(cwxco)
}
