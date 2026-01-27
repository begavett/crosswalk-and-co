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
p_load(magrittr, dplyr, forcats, ggplot2, tidyr, data.table, patchwork, ggtext)

user <- Sys.info()[["user"]]

if (user == "brandon") {
  code_dir <- "~/Dropbox/Projects/crosswalk-and-co/"
} else if (user == "emmanich") {
  if (Sys.info()[["sysname"]] == "Windows") {
    code_dir <- "C:/Users/emmanich/code/crosswalk-and-co/"
  } else {
    code_dir <- "/Users/emmanich/code/crosswalk-and-co/"
  }
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
alldata[, fill_factor := fct_cross(Method, crosswalk_to, sep = " - ")]

## make absolute value versions of bias  
alldata[, `:=` (abs_bias_sim = abs(bias_sim), abs_bias_truth = abs(bias_truth))]

maindata <- copy(alldata)

# COMPARISONS TO TRUE REGRESSION COEFFICIENTS ----------------------------

get_boxplot <- function(yvar = "bias_sim", dataset = maindata, cw_to = "G2", ylab = "Bias"){
    if (!is.na(cw_to)){
        dataset <- copy(dataset[crosswalk_to == cw_to])
    }
    lims <- c(dataset[, min(get(yvar))], dataset[, max(get(yvar))])
    p <- ggplot(dataset, aes(x = x_factor, y = get(yvar), fill = fill_factor)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab(ylab) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        #scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_fill_manual(name = "", values = c("#4287f5", "#4287f5", "#f06311", "#f06311")) +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

trueplot <- get_boxplot(yvar = "bias_sim", cw_to = NA, ylab = "Absolute Bias \n(estimate vs. sim parameter)")

ggsave(paste0(code_dir, "plots/bias_vs_simcoef_", date, ".jpg"), trueplot, width = 14, height = 6)

# COMPARISONS TO TARGET COEFFICIENTS ---------------------------------------

targetplot <- get_boxplot(yvar = "bias_truth", cw_to = NA, ylab = "Absolute Bias \n(estimate vs. target coef)")

ggsave(paste0(code_dir, "plots/bias_vs_targetcoef_", date, ".jpg"), targetplot, width = 14, height = 6)

targetplot_aaic <- targetplot + 
    ggtitle("**Figure 2.** Bias of cogxwalkr crosswalk and co-calibration methods across simulation replications.") +  
    labs(caption = "Note: ES = effect size.") + 
    theme(plot.caption.position = "plot", plot.title.position = "plot", 
        plot.title = element_textbox_simple(),
        plot.caption = element_textbox_simple(width = unit(1, "npc"), margin = margin(t = 10))) 

ggsave(paste0(code_dir, "plots/bias_vs_targetcoef_aaic_", date, ".jpg"), targetplot_aaic, width = 12, height = 5)

# BARPLOT COMPARISON --------------------------------------------------------------

## for a given iteration number and outcome 

get_meanplot <- function(yvar, dataset = maindata, cw_to = "G2"){
    if (!is.na(cw_to)){
        dataset <- copy(dataset[crosswalk_to == cw_to])
    }
    lims <- c(dataset[, min(get(yvar))], dataset[, max(get(yvar))])
    p <- ggplot(dataset, aes(x = x_factor, y = get(yvar), color = fill_factor)) +
        stat_summary(fun = mean,
                     geom = "point",
                     position = position_dodge(width = 0.9)) +
        stat_summary(fun.data = mean_cl_normal,
                     geom = "errorbar",
                     position = position_dodge(width = 0.9),
                     width = 0.2) +
        facet_wrap(~ slabel, nrow = 1) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        #scale_color_viridis_d(name = "", option = "turbo", begin = .2, end = .8) +
        scale_color_manual(name = "", values = c("#4287f5", "#4287f5", "#f06311", "#f06311")) +
        ylab("Mean Bias (95% CI)") +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

meanplots <- lapply(c("bias_sim", "bias_truth"), function(y) get_meanplot(y, cw_to = NA))

meanplot <- wrap_plots(meanplots) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(code_dir, "plots/meanplots_", date, ".jpg"), meanplots, width = 14, height = 8)

barplot_aaic <- meanplots[[2]] + 
    ggtitle("**Figure 3.** Mean bias of cogxwalkr crosswalk and co-calibration.") +  
    theme(plot.caption.position = "plot", plot.title.position = "plot", 
        plot.title = element_textbox_simple(),
        plot.caption = element_textbox_simple(width = unit(1, "npc"), margin = margin(t = 10))) 

ggsave(paste0(code_dir, "plots/barplot_aaic_", date, ".jpg"), barplot_aaic, width = 12, height = 5)

# CHECK OTHER DIRECTION -----------------------------------------------------------

barplots_otherdirection <- lapply(c("bias_sim", "bias_truth"), function(y) get_barplot(y, cw_to = "G1"))

barplot_otherdirection <- wrap_plots(barplots_otherdirection) + 
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(code_dir, "plots/barplots_otherdirection_", date, ".jpg"), barplot_otherdirection, width = 14, height = 8)

# COMPARISON WITH MORE ITERATIONS -------------------------------------------------

it1000_dt <- copy(maindata[rep <= 1000]); it2000_dt <- copy(maindata)

iteration_dt <- rbind(it1000_dt[, version := "1000 iterations"], it2000_dt[, version := "2000 iterations"])


get_iterplot <- function(yvar){
    data <- copy(iteration_dt)
    lims <- c(iteration_dt[, min(get(yvar))], iteration_dt[, max(get(yvar))])
    p <- ggplot(data, aes(x = x_factor, y = get(yvar), fill = as.factor(version))) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel+Method, nrow = 2) +
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

ggsave(paste0(code_dir, "plots/iterplot_sim_", date, ".jpg"), iterplots[[1]], width = 12, height = 8)


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

