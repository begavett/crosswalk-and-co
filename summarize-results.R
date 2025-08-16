# ------------------------------------------------------------------------------
# Simulation Results Visualization Script
#
# Description:
#   This script loads simulation results from RDS files and generates plots
#   to evaluate bias in regression coefficient estimates under different
#   crosswalking methods, sample sizes, and scenarios.
#
# Data Inputs:
#   - allres5002.RDS : Results from simulations with sample size N = 5,002
#   - allres1002.RDS : Results from simulations with sample size N = 1,002
#
# Process:
#   1. Load results and create derived variables:
#        * bias  = estimated coefficient - 0.2 (true effect size)
#        * bias2 = estimated coefficient - truecoefs (scenario-specific target)
#        * relabel "Crosswalk" method to "cogxwalkr" for clarity
#   2. For each dataset (N = 5002 and N = 1002):
#        * Plot distribution of bias vs. truth
#        * Plot distribution of bias vs. target coefficient
#   3. For N = 1002 specifically:
#        * Compute and plot mean bias ± 1 SE across scenarios, by method
#
# Outputs:
#   - Boxplots of bias by scenario, method, and target group
#   - Barplots with error bars showing mean bias across replications
#
# Key Variables:
#   - coef       : estimated regression coefficient
#   - truecoefs  : scenario-specific target coefficient
#   - rep        : replication ID
#   - slabel     : scenario label
#   - Method     : estimation method (e.g., cogxwalkr, etc.)
#   - crosswalk_to : target group for crosswalk
#
# Notes:
#   - Bias is evaluated relative to both the "true" value (0.2) 
#     and the scenario-specific target coefficient.
#   - Iteration counts (niter) are automatically included in plot subtitles.
#   - Uses ggplot2 for visualization and viridis color scale for accessibility.
# ------------------------------------------------------------------------------


library(pacman)
p_load(magrittr, dplyr, forcats, ggplot2, tidyr)

# N = 5002 ----------------------------------------------------------------
allres5002 <- readRDS("~/Dropbox/Projects/crosswalk-and-co/data/allres5002.RDS")

results_df <- allres5002 %>%
  mutate(bias = coef - .2,
         bias2 = coef - truecoefs,
         slabel = fct_inorder(slabel),
         Method = case_when(Method == "Crosswalk" ~ "cogxwalkr",
                            TRUE ~ Method))

niter <- length(unique(results_df$rep))


## Bias vs. Truth ----------------------------------------------------------

results_df %>%
  ggplot(aes(x = crosswalk_to, y = bias, fill = Method)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~slabel) +
  ylab("Bias") +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Bias (vs. Truth) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 5,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  scale_y_continuous(limits = c(-.3, .3), n.breaks = 6) +
  geom_hline(yintercept = 0, colour = "red", lty = 2)

## Bias vs. Target ----------------------------------------------------------

results_df %>%
  ggplot(aes(x = crosswalk_to, y = bias2, fill = Method)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~slabel) +
  ylab("Bias") +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Bias (vs. Target Coefficient) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 5,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  scale_y_continuous(limits = c(-.3, .3), n.breaks = 6) +
  geom_hline(yintercept = 0, colour = "red", lty = 2)


# N = 1002 ----------------------------------------------------------------

allres1002 <- readRDS("~/Dropbox/Projects/crosswalk-and-co/data/allres1002.RDS")

results_df <- allres1002 %>%
  mutate(bias = coef - .2,
         bias2 = coef - truecoefs,
         slabel = fct_inorder(slabel),
         Method = case_when(Method == "Crosswalk" ~ "cogxwalkr",
                            TRUE ~ Method))

niter <- length(unique(results_df$rep))

## Bias vs. Truth ----------------------------------------------------------

results_df %>%
  ggplot(aes(x = crosswalk_to, y = bias, fill = Method)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~slabel) +
  ylab("Bias") +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Bias (vs. Truth) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 1,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  scale_y_continuous(limits = c(-.3, .3), n.breaks = 6) +
  geom_hline(yintercept = 0, colour = "red", lty = 2)

## Bias vs. Target ----------------------------------------------------------

results_df %>%
  ggplot(aes(x = crosswalk_to, y = bias2, fill = Method)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~slabel) +
  ylab("Bias") +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Bias (vs. Target Coefficient) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 1,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  scale_y_continuous(limits = c(-.3, .3), n.breaks = 6) +
  geom_hline(yintercept = 0, colour = "red", lty = 2)


## Mean Bias vs. Truth ----------------------------------------------------------

allres1002 %>%
  mutate(bias = coef - .2, 
         bias2 = coef - truecoefs) %>%
  select(crosswalk_to, slabel, Method, n_sample, bias, bias2) %>%
  pivot_longer(cols = c(bias, bias2)) %>%
  mutate(name = case_when(name == "bias" ~ "Bias vs. True",
                          name == "bias2" ~ "Bias vs. Target"),
         slabel = fct_inorder(slabel)) %>%
  filter(name == "Bias vs. True") %>%
  ggplot(aes(x = crosswalk_to, y = value, fill = Method)) +
  stat_summary(fun = mean,
               geom = "bar",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               position = position_dodge(width = 0.9),
               width = 0.2) +
  facet_wrap(~ slabel) +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean Bias (vs. True) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 1,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  ylab("Mean Bias (± 1 SE)") +
  geom_hline(yintercept = 0, colour = "red", lty = 2)

## Mean Bias vs. Target ----------------------------------------------------------

allres1002 %>%
  mutate(bias = coef - .2, 
         bias2 = coef - truecoefs) %>%
  select(crosswalk_to, slabel, Method, n_sample, bias, bias2) %>%
  pivot_longer(cols = c(bias, bias2)) %>%
  mutate(name = case_when(name == "bias" ~ "Bias vs. True",
                          name == "bias2" ~ "Bias vs. Target"),
         slabel = fct_inorder(slabel)) %>%
  filter(name == "Bias vs. Target") %>%
  ggplot(aes(x = crosswalk_to, y = value, fill = Method)) +
  stat_summary(fun = mean,
               geom = "bar",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               position = position_dodge(width = 0.9),
               width = 0.2) +
  facet_wrap(~ slabel) +
  xlab("Target Group") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Mean Bias (vs. Target) by Scenario, Method, and Target Group",
          subtitle = paste0("Results from ", niter, " iterations, each iteration N = 1,002")) +
  scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
  ylab("Mean Bias (± 1 SE)") +
  geom_hline(yintercept = 0, colour = "red", lty = 2)
