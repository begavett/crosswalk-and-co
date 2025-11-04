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
p_load(magrittr, dplyr, forcats, ggplot2, tidyr, data.table, patchwork)

user <- Sys.info()[["user"]]

if (user == "brandon") {
  code_dir <- "~/Dropbox/Projects/crosswalk-and-co/"
} else if (user == "emmanich") {
  code_dir <- "C:/Users/emmanich/code/crosswalk-and-co/"
}

# LOAD FILES -----------------------------------------------------------

files <- list.files(paste0(code_dir, "models/"), pattern = "allresn", full.names = TRUE)
alldata <- rbindlist(lapply(files, function(file) as.data.table(readr::read_rds(file)))) 

# FORMAT DATA -----------------------------------------------------------

## for cocalibration - scaling should be from the version that you are crosswalking to - so only 
## keep the version that where the reference matches the version you are crosswalking to 
alldata <- alldata[!(Method == "Cocalibration" & !crosswalk_to == cc_rg)]

## define bias and other variables 
alldata[, `:=` (bias_sim = coef-b1, bias_truth = coef-truecoefs, naive_diff = (naivecoefs-truecoefs) - (coef-truecoefs))]
alldata[Method == "cogxwalkr", Method := "Cogxwalkr"][, Method := factor(Method, levels = c("Cocalibration", "Cogxwalkr"))]
alldata[, `:=` (slabel = fct_inorder(slabel), eslabel = factor(paste0("ES ", b1), levels = c("ES 0.1", "ES 0.2")), 
                crosswalk_to = factor(gsub("Group ", "G", crosswalk_to), levels = c("G1", "G2")))]
alldata[, x_factor := fct_cross(crosswalk_to, eslabel, sep = " - ")]

## remove the extra simulation sensitivty (put aside)
maindata <- copy(alldata[rep <=1000])
simsens_dt <- copy(alldata[rep >1000])

# COMPARISONS TO TRUE REGRESSION COEFFICIENTS ----------------------------

## summarizing over group 1 and group 2 - no real difference - also can then show different effect sizes

get_trueplot <- function(it){
    data <- copy(maindata[n_sample == it])
    lims <- c(maindata[, min(bias_sim)], maindata[, max(bias_sim)])
    p <- ggplot(data, aes(x = x_factor, y = bias_sim, fill = Method)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab("Bias") +
        xlab("Target Group") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

trueplots <- lapply(c(500, 1000, 5000), function(n) get_trueplot(n))

trueplot <- wrap_plots(trueplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

# COMPARISONS TO TARGET COEFFICIENTS ---------------------------------------

## summarizing over group 1 and group 2 - no real difference - also can then show different effect sizes

get_targetplot <- function(it){
    data <- copy(maindata[n_sample == it])
    lims <- c(maindata[, min(bias_truth)], maindata[, max(bias_truth)])
    p <- ggplot(data, aes(x = x_factor, y = bias_truth, fill = Method)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab("Bias") +
        xlab("Target Group") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

targetplots <- lapply(c(500, 1000, 5000), function(n) get_targetplot(n))

targetplot <- wrap_plots(targetplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

# BARPLOT COMPARISON --------------------------------------------------------------

## for a given iteration number and outcome 

get_barplot <- function(yvar, it){
    data <- copy(maindata[n_sample == it])
    lims <- c(maindata[, min(get(yvar))], maindata[, max(get(yvar))])
    p <- ggplot(data, aes(x = x_factor, y = get(yvar), fill = Method)) +
        stat_summary(fun = mean,
                     geom = "bar",
                     position = position_dodge(width = 0.9)) +
        stat_summary(fun.data = mean_se,
                     geom = "errorbar",
                     position = position_dodge(width = 0.9),
                     width = 0.2) +
        facet_wrap(~ slabel, nrow = 1) +
        xlab("Target Group") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        ylab("Mean Bias (± 1 SE)") +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

barplots1000 <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y, 1000))

barplot1000 <- wrap_plots(barplots1000) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

barplots500 <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y, 500))
barplot500 <- wrap_plots(barplots500) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

barplots5000 <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y, 5000))
barplot5000 <- wrap_plots(barplots5000) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

# COMPARISON WITH MORE ITERATIONS -------------------------------------------------

## FINISH THIS TOMORROW

