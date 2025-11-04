###############################################################################
# Simulation Results Summary and Plotting Script
#
# Purpose:
#   This script loads simulation results from RDS files, computes bias metrics
#   relative to either the true coefficient or a target coefficient, and 
#   generates publication-ready plots (boxplots or barplots) comparing methods
#   across scenarios and target groups.
#
# Libraries:
#   - pacman: for easy package loading
#   - magrittr: for piping (%>%)
#   - dplyr: for data manipulation
#   - forcats: for factor level manipulation
#   - ggplot2: for plotting
#   - tidyr: for data reshaping
#   - patchwork: for combining ggplots
#
# Main Function:
#   summarize_sims(outfilename, ref_grp = 1, plot_type = c("boxplot", "barplot"),
#                  vs = c("truth", "target"))
#
# Arguments:
#   outfilename : character
#       Name of the RDS file containing simulation results.
#   ref_grp     : integer (1 or 2)
#       Reference group used for cocalibration.
#   plot_type   : character
#       Type of plot to generate. Options: "boxplot" (individual iteration
#       results) or "barplot" (mean ± SE across iterations).
#   vs          : character
#       Comparison reference. Options: "truth" (compare to true coefficient)
#       or "target" (compare to target coefficient).
#
# Output:
#   - Generates a ggplot showing bias by scenario, method, and target group.
#   - Returns a summarized data frame with mean values across iterations.
#
# Usage:
#   # Boxplot of bias vs. true coefficient with Group 1 as reference
#   summarize_sims("allresn5000_d1.Rds", ref_grp = 1, plot_type = "boxplot", vs = "truth")
#
# Notes:
#   - The script assumes RDS files are located in `file.path(code_dir, "data")`.
#   - Bias is computed as difference between estimated coefficients and
#     either the true coefficient or the target coefficient.
#   - Plots use viridis color scales and facet_wrap by scenario label.
#   - Horizontal red dashed line indicates zero bias.
##########################################################

library(pacman)
p_load(magrittr, dplyr, forcats, ggplot2, tidyr, patchwork)

user <- Sys.info()[["user"]]

if (user == "brandon") {
  code_dir <- "~/Dropbox/Projects/crosswalk-and-co/"
} else if (user == "emmanich") {
  code_dir <- "C:/Users/emmanich/code/crosswalk-and-co/"
}


summarize_sims <- function(outfilename, 
                           ref_grp = 1, 
                           plot_type = c("boxplot", "barplot"),
                           vs = c("truth", "target")) {
  
  res <- readRDS(file.path(code_dir, "data", outfilename))
  beta1 <- first(res$b1)
  plot_type <- match.arg(plot_type)
  vs <- match.arg(vs)
  
  results_df <- res %>%
    filter(cc_rg == paste0("Group ", ref_grp) | Method == "cogxwalkr") %>%
    mutate(bias = coef - b1,
           bias2 = coef - truecoefs,
           slabel = fct_inorder(slabel),
           Method = case_when(Method == "Crosswalk" ~ "cogxwalkr",
                              TRUE ~ Method))
  
  niter <- length(unique(results_df$rep))
  
  
  if (plot_type == "boxplot") {
    
    if(vs == "truth") {
      ## Bias vs. Truth ----------------------------------------------------------
      
      p <- results_df %>%
        ggplot(aes(x = crosswalk_to, y = bias, fill = Method)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel) +
        ylab("Bias") +
        xlab("Target Group") +
        theme_bw() +
        theme(legend.position = "bottom") +
        ggtitle(paste0("Bias (vs. Truth) by Scenario, Method, and Target Group (\u03B2 = ", beta1, ")"),
                subtitle = paste0("Results from ", niter, " iterations, each iteration N = ", first(results_df$n_sample))) +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_y_continuous(limits = c(min(min(results_df$bias), min(results_df$bias2)), max(max(results_df$bias), max(results_df$bias2))), n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2) +
        plot_annotation(caption = paste0("Group ", ref_grp, " was the reference group for cocalibration"))
      
    } else if (vs == "target") {
      
      ## Bias vs. Target ----------------------------------------------------------
      
      p <- results_df %>%
        ggplot(aes(x = crosswalk_to, y = bias2, fill = Method)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel) +
        ylab("Bias") +
        xlab("Target Group") +
        theme_bw() +
        theme(legend.position = "bottom") +
        ggtitle(paste0("Bias (vs. Target Coefficient) by Scenario, Method, and Target Group (\u03B2 = ", beta1, ")"),
                subtitle = paste0("Results from ", niter, " iterations, each iteration N = ", first(results_df$n_sample))) +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_y_continuous(limits = c(min(min(results_df$bias), min(results_df$bias2)), max(max(results_df$bias), max(results_df$bias2))), n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2) +
        plot_annotation(caption = paste0("Group ", ref_grp, " was the reference group for cocalibration"))
    }
  } else if (plot_type == "barplot") {
    
    if (vs == "truth") {
      
      ## Mean Bias vs. Truth ----------------------------------------------------------
      
      p <- results_df %>%
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
        ggtitle(paste0("Mean Bias (vs. True) by Scenario, Method, and Target Group (\u03B2 = ", beta1, ")"),
                subtitle = paste0("Results from ", niter, " iterations, each iteration N = ", first(results_df$n_sample))) +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        ylab("Mean Bias (± 1 SE)") +
        geom_hline(yintercept = 0, colour = "red", lty = 2) +
        plot_annotation(caption = paste0("Group ", ref_grp, " was the reference group for cocalibration"))
      
    } else if (vs == "target") {
      
      ## Mean Bias vs. Target ----------------------------------------------------------
      
      p <- results_df %>%
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
        ggtitle(paste0("Mean Bias (vs. Target) by Scenario, Method, and Target Group (\u03B2 = ", beta1, ")"),
                subtitle = paste0("Results from ", niter, " iterations, each iteration N = ", first(results_df$n_sample))) +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        ylab("Mean Bias (± 1 SE)") +
        geom_hline(yintercept = 0, colour = "red", lty = 2) +
        plot_annotation(caption = paste0("Group ", ref_grp, " was the reference group for cocalibration"))
    }
  }
  plot(p)
  return(results_df %>%
           select(-rep) %>%
           group_by(scenario, crosswalk_to, slabel, cc_rg, Method, n_sample, b1) %>%
           summarise(across(everything(), mean),
                     .groups = "drop"))
}





summarize_sims(outfilename = "allresn5000_d1.Rds",
               ref_grp = 1, 
               plot_type = "boxplot",
               vs = "target")
