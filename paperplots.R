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

date <- format(Sys.Date(), "%Y_%m_%d")

# LOAD FILES -----------------------------------------------------------

files <- list.files(paste0(code_dir, "models/"), pattern = "allresn[0-9]*_d[0-9]_sum", full.names = TRUE)
alldata <- rbindlist(lapply(files, function(file) as.data.table(readr::read_rds(file)))) 

# FORMAT DATA -----------------------------------------------------------

## for cocalibration - scaling should be from the version that you are crosswalking to - so only 
## keep the version that where the reference matches the version you are crosswalking to 
alldata <- alldata[!(Method == "Cocalibration" & !crosswalk_to == cc_rg)]

## define bias and other variables 
alldata[, `:=` (bias_sim = coef-b1, bias_truth = coef-truecoefs, naive_diff = abs(naivecoefs-truecoefs) - abs(coef-truecoefs), 
                naivetrue_diff = naivecoefs - truecoefs)]
alldata[Method == "cogxwalkr", Method := "Cogxwalkr"][, Method := factor(Method, levels = c("Cocalibration", "Cogxwalkr"))]
alldata[, `:=` (slabel = fct_inorder(slabel), eslabel = factor(paste0("ES ", b1), levels = c("ES 0.2", "ES 0.4")), 
                crosswalk_to = factor(gsub("Group ", "G", crosswalk_to), levels = c("G1", "G2")), 
                N_label = factor(paste0("N=", n_sample), levels = c("N=500", "N=1000", "N=5000")))]
alldata[, x_factor := fct_cross(eslabel, N_label, sep = "\n")]

maindata <- copy(alldata)

# COMPARISONS TO TRUE REGRESSION COEFFICIENTS ----------------------------

get_boxplot <- function(yvar = "bias_sim", dataset = maindata, cw_to = "G2", ylab = "Bias"){
    data <- copy(dataset[crosswalk_to == cw_to])
    lims <- c(dataset[, min(get(yvar))], dataset[, max(get(yvar))])
    p <- ggplot(data, aes(x = x_factor, y = get(yvar), fill = Method)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab(ylab) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

trueplot <- get_boxplot(yvar = "bias_sim", ylab = "Bias \n(estimate vs. sim parameter)")

ggsave(paste0(code_dir, "plots/bias_vs_simcoef_", date, ".jpg"), trueplot, width = 14, height = 6)

# COMPARISONS TO TARGET COEFFICIENTS ---------------------------------------

targetplot <- get_boxplot(yvar = "bias_truth", ylab = "Bias \n(estimate vs. target coef)")

ggsave(paste0(code_dir, "plots/bias_vs_targetcoef_", date, ".jpg"), targetplot, width = 14, height = 6)

# BARPLOT COMPARISON --------------------------------------------------------------

## for a given iteration number and outcome 

get_barplot <- function(yvar, dataset = maindata, cw_to = "G2"){
    data <- copy(dataset[crosswalk_to == cw_to])
    lims <- c(dataset[, min(get(yvar))], dataset[, max(get(yvar))])
    p <- ggplot(data, aes(x = x_factor, y = get(yvar), fill = Method)) +
        stat_summary(fun = mean,
                     geom = "bar",
                     position = position_dodge(width = 0.9)) +
        stat_summary(fun.data = mean_se,
                     geom = "errorbar",
                     position = position_dodge(width = 0.9),
                     width = 0.2) +
        facet_wrap(~ slabel, nrow = 1) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        ylab("Mean Bias (± 1 SE)") +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

barplots <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y))

barplot <- wrap_plots(barplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(code_dir, "plots/barplots_", date, ".jpg"), barplot, width = 14, height = 8)

# CHECK OTHER DIRECTION -----------------------------------------------------------

barplots_otherdirection <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y, cw_to = "G1"))

barplot_otherdirection <- wrap_plots(barplots_otherdirection) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(code_dir, "plots/barplots_otherdirection_", date, ".jpg"), barplot_otherdirection, width = 14, height = 8)

# COMPARISON WITH MORE ITERATIONS -------------------------------------------------

iteration_dt <- rbind(maindata[b1 == 0.2 & n_sample == 1000, .(slabel, crosswalk_to, Method, bias_sim, bias_truth, version = "1000 iterations")], 
                      rbind(maindata[b1 == 0.2 & n_sample == 1000, .(slabel, crosswalk_to, Method, bias_sim, bias_truth, version = "2000 iterations")], 
                            simsens_dt[b1 == 0.2 & n_sample == 1000, .(slabel, crosswalk_to, Method, bias_sim, bias_truth, version = "2000 iterations")]))
iteration_dt[, method := factor(fcase(Method == "Cocalibration", "Cocal.", 
                                        Method == "Cogxwalkr", "Cogx."), levels = c("Cocal.", "Cogx."))]
iteration_dt[, x_factor := fct_cross(crosswalk_to, method, sep = " - ")]


get_iterplot <- function(yvar){
    data <- copy(iteration_dt)
    lims <- c(iteration_dt[, min(get(yvar))], iteration_dt[, max(get(yvar))])
    p <- ggplot(data, aes(x = x_factor, y = get(yvar), fill = as.factor(version))) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab("Bias") +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_fill_viridis_d(option = "plasma", begin = .2, end = .8, name = "") +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

iterplots <- lapply(c("bias_sim", "bias_truth"), function(y) get_iterplot(y))

iterplot <- wrap_plots(iterplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

# NAIVE COMPARISON -----------------------------------------------------------

naiveplots_box <- lapply(c(500, 1000, 5000), function(n) get_boxplot(n, yvar = "naive_diff") + ylab("Naive Bias - Estimated Bias"))

naiveplot_box <- wrap_plots(naiveplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

naiveplots_bar <- lapply(c(500, 1000, 5000), function(n) get_barplot("naive_diff", n) + ylab("Naive Bias - Estimated Bias"))

naiveplot_bar <- wrap_plots(naiveplots_bar) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")
## reason why we see low means for cocalibration: 
maindata[, .(min = min(naive_diff), max = max(naive_diff)), by = c("x_factor", "Method")]
ggplot(maindata[Method == "Cocalibration" & n_sample == 1000 & x_factor == "G1 - ES 0.2"], aes(x = coef)) + geom_histogram(fill = "white", color = "black") + theme_bw()

## boxplot for just cogxwalkr 
naiveplots_cogx_box <- lapply(c(500, 1000, 5000), function(n) get_boxplot(n, yvar = "naive_diff", dataset = maindata[Method == "Cogxwalkr"]) + ylab("Naive Bias - Estimated Bias"))

naiveplot_cogx_box <- wrap_plots(naiveplots_cogx_box) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

## barplot for just cogxwalkr
naiveplots_cogx_bar <- lapply(c(500, 1000, 5000), function(n) get_barplot("naive_diff", n, dataset = maindata[Method == "Cogxwalkr"]) + ylab("Naive Bias - Estimated Bias"))

naiveplot_cogx_bar <- wrap_plots(naiveplots_cogx_bar) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

## difference between naive and true effects (is this not really that big?)
naivetrueplots_box <- lapply(c(500, 1000, 5000), function(n) get_boxplot(n, yvar = "naivetrue_diff", dataset = maindata[Method == "Cogxwalkr"]) + ylab("Naive - True Coefficient") + theme(legend.position = "none"))

naivetrueplot_box <- wrap_plots(naivetrueplots_box) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".")

