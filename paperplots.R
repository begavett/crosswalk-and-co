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
p_load(magrittr, dplyr, forcats, ggplot2, tidyr, data.table, patchwork, ggtext, binom)

user <- Sys.info()[["user"]]

if (user == "brandon") {
    code_dir <- "~/Dropbox/Projects/crosswalk-and-co/"
} else if (user == "emmanich") {
    if (Sys.info()[["sysname"]] == "Windows") {
        code_dir <- "C:/Users/emmanich.HP/code/crosswalk-and-co/"
    } else {
        code_dir <- "/Users/emmanich/code/crosswalk-and-co/"
    }
}

date <- format(Sys.Date(), "%Y_%m_%d")

# LOAD FILES -----------------------------------------------------------

files <- list.files(paste0(code_dir, "models/"), pattern = "allresn[0-9]*_d[0-9]_ci_aug2026", full.names = TRUE)
alldata <- rbindlist(lapply(files, function(file) as.data.table(readr::read_rds(file))))

# FORMAT DATA -----------------------------------------------------------

format_data <- function(data, CI = FALSE) {
    ## for cocalibration - scaling should be from the version that you are crosswalking to - so only
    ## keep the version that where the reference matches the version you are crosswalking to
    data <- data[!(Method == "Cocalibration" & !crosswalk_to == cc_rg)]

    ## define bias and other variables
    data[, `:=`(
        bias_sim = coef - b1, bias_truth = coef - truecoefs, naive_diff = abs(naivecoefs - truecoefs) - abs(coef - truecoefs),
        naivetrue_diff = naivecoefs - truecoefs
    )]
    data[Method == "cogxwalkr", Method := "ES Crosswalking"][, Method := factor(Method, levels = c("Cocalibration", "ES Crosswalking"))]
    data[, `:=`(
        slabel = fct_inorder(slabel), eslabel = factor(paste0("ES ", b1), levels = c("ES 0.2", "ES 0.4")),
        crosswalk_to = factor(gsub("Group ", "O", crosswalk_to), levels = c("O1", "O2")),
        N_label = factor(paste0("N=", n_sample), levels = c("N=500", "N=1000", "N=5000"))
    )]
    data[, x_factor := fct_cross(eslabel, N_label, sep = "\n")]
    data[, fill_factor := fct_cross(Method, crosswalk_to, sep = " - ")]

    ## make absolute value versions of bias
    data[, `:=`(abs_bias_sim = abs(bias_sim), abs_bias_truth = abs(bias_truth))]

    ## education conditions
    data[, dementia_generation_normal := as.numeric(demwithedu == FALSE | is.na(demwithedu))]
    data[, dementia_generation := fcase(
        demwithedu == TRUE, "With education",
        demwithedu == FALSE, "Standard",
        is.na(demwithedu), NA_character_
    )]
    data[, dementia_generation := fct_cross(dementia_generation, crosswalk_to, sep = " - ")]

    ## get percent bias
    data[, `:=`(pct_bias_sim = bias_sim / b1, pct_bias_truth = bias_truth / b1)]

    ## CI coverage and width if exists
    if (CI == TRUE) {
        data[, CIcoverage := as.numeric(truecoefs > lwr & truecoefs < upr)]
        data[, CIcoverage_sim := as.numeric(b1 > lwr & b1 < upr)]
        data[, CIwidth := upr - lwr]
    }

    return(data)
}

maindata <- format_data(alldata, CI = TRUE)

# COMPARISONS TO TRUE REGRESSION COEFFICIENTS ----------------------------

get_boxplot <- function(yvar = "bias_sim", dataset = maindata, cw_to = "O2", ylab = "Bias") {
    if (!is.na(cw_to)) {
        dataset <- copy(dataset[crosswalk_to == cw_to])
    }
    dataset <- dataset[dementia_generation_normal == TRUE]
    lims <- c(dataset[, min(get(yvar))], dataset[, max(get(yvar))])
    p <- ggplot(dataset, aes(x = x_factor, y = get(yvar), fill = fill_factor)) +
        geom_boxplot(notch = TRUE) +
        facet_wrap(~slabel, nrow = 1) +
        ylab(ylab) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        # scale_fill_viridis_d(option = "turbo", begin = .2, end = .8) +
        scale_fill_manual(name = "", values = c("#5e7cad", "#1b3c71", "#e78147", "#d25309")) +
        scale_y_continuous(limits = lims, n.breaks = 6) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

trueplot <- get_boxplot(yvar = "bias_sim", cw_to = NA, ylab = "Absolute Bias \n(estimate vs. sim parameter)")

# COMPARISONS TO TARGET COEFFICIENTS ---------------------------------------

targetplot <- get_boxplot(yvar = "bias_truth", cw_to = NA, ylab = "Absolute Bias \n(estimate vs. target coef)")

# JOINT FIGURE 1 --------------------------------------------------------------

fig1 <- trueplot / targetplot +
    plot_annotation(tag_levels = "A", tag_suffix = ".") +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(code_dir, "plots/Figure2_distributions_", date, ".pdf"), fig1, width = 14, height = 9)

# MEAN COMPARISON --------------------------------------------------------------

## for a given iteration number and outcome

#### observing super-population
### SD of bias of estimates - divided by square root of sample size

## add percent coverage for the confidence intervals

get_meanplot <- function(yvar, dataset = maindata, cw_to = "O2", color_var = "fill_factor", percent_labels = TRUE,
                         ylabel = "Mean Bias (95% CI)") {
    if (!is.na(cw_to)) {
        dataset <- copy(dataset[crosswalk_to == cw_to])
    }
    if (color_var == "fill_factor") {
        dataset <- dataset[dementia_generation_normal == TRUE]
        cols <- c("#5e7cad", "#1b3c71", "#e78147", "#d25309")
    } else {
        dataset <- dataset[!is.na(get(color_var))]
        cols <- c("#e78147", "#d25309", "#a1f164", "#458c0f")
    }
    # Precompute groupwise mean, sd, n, se and t-based 95% CI for the chosen yvar
    summary_dt <- dataset[, .(
        mean = mean(get(yvar), na.rm = TRUE),
        sd = sd(get(yvar), na.rm = TRUE),
        n = .N
    ), by = c("slabel", "x_factor", color_var)]

    summary_dt[, se := sd / sqrt(n)]
    # Use t critical value with df = n-1 (guard df >= 1)
    summary_dt[, `:=`(
        tcrit = qt(0.975, df = pmax(n - 1, 1)),
        lower = mean - qt(0.975, df = pmax(n - 1, 1)) * se,
        upper = mean + qt(0.975, df = pmax(n - 1, 1)) * se
    )]

    p <- ggplot(summary_dt, aes(x = x_factor, y = mean, color = get(color_var), shape = get(color_var))) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
        facet_wrap(~slabel, nrow = 1) +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_shape_manual(name = "", values = c(16, 16, 17, 17)) +
        scale_y_continuous(
            labels = if (percent_labels) scales::label_percent() else ggplot2::waiver(),
            breaks = scales::breaks_pretty(n = 5)
        ) +
        # scale_color_viridis_d(name = "", option = "turbo", begin = .2, end = .8) +
        scale_color_manual(name = "", values = cols) +
        ylab(ylabel) +
        geom_hline(yintercept = 0, colour = "red", lty = 2)
    return(p)
}

meanplots <- lapply(c("pct_bias_sim", "pct_bias_truth"), function(y) get_meanplot(y, cw_to = NA))

ggsave(paste0(code_dir, "plots/Figure3_meanplots_", date, ".pdf"), meanplots[[2]], width = 14, height = 5)
ggsave(paste0(code_dir, "plots/FigureS1_meanplots_sim_", date, ".pdf"), meanplots[[1]], width = 14, height = 5)


# EDUCATION COMPARISON -----------------------------------------------------

education_meanplot_comparison <- get_meanplot(yvar = "pct_bias_truth", cw_to = NA, color_var = "dementia_generation")

# CI COVERAGE ---------------------------------------------------------------------

get_coverageplot <- function(yvar = "CIcoverage", dataset = maindata, cw_to = "O2") {
    if (!is.na(cw_to)) {
        dataset <- copy(dataset[crosswalk_to == cw_to])
    }
    dataset <- dataset[dementia_generation_normal == TRUE]

    group_vars <- c("slabel", "x_factor", "fill_factor")
    coverage_dt <- dataset[,
        {
            covered <- get(yvar)

            .(
                num_coverage = sum(covered, na.rm = TRUE),
                num_sims = sum(!is.na(covered)),
                num_missing = sum(is.na(covered))
            )
        },
        by = group_vars
    ]

    coverage_dt[, pct_coverage := num_coverage / num_sims]

    coverage_ci <- binom.confint(
        x = coverage_dt$num_coverage,
        n = coverage_dt$num_sims,
        conf.level = 0.95,
        methods = "wilson"
    )

    coverage_dt[, `:=`(
        lower = coverage_ci$lower,
        upper = coverage_ci$upper
    )]

    p <- ggplot(coverage_dt, aes(
        x = x_factor, y = pct_coverage, ymin = lower,
        ymax = upper, color = fill_factor, shape = fill_factor
    )) +
        geom_point(position = position_dodge(width = 0.4)) +
        geom_errorbar(width = 0, position = position_dodge(width = 0.4)) +
        geom_hline(yintercept = 0.95, colour = "red", lty = 2) +
        facet_wrap(~slabel, nrow = 1) +
        ylab("Coverage") +
        xlab("Condition") +
        theme_bw() +
        theme(legend.position = "bottom") +
        scale_y_continuous(
            limits = c(0.01, 1),
            breaks = c(0, seq(0.75, 1, by = 0.05)),
            labels = scales::label_percent()
        ) +
        ggbreak::scale_y_break(
            c(0, 0.75)
        ) +
        scale_shape_manual(name = "", values = c(16, 16, 17, 17)) +
        scale_color_manual(name = "", values = c("#5e7cad", "#1b3c71", "#e78147", "#d25309"))
    return(p)
}

coverageplot <- get_coverageplot(cw_to = NA)
coverageplot_sim <- get_coverageplot(cw_to = NA, yvar = "CIcoverage_sim")

ggsave(paste0(code_dir, "plots/Figure4_coverageplot_", date, ".pdf"), coverageplot, width = 14, height = 5)

# CI WIDTH ------------------------------------------------------------------

maindata[, .(ciwidth = mean(CIwidth), cicoverage = mean(CIcoverage)), by = c("slabel", "x_factor", "fill_factor")]
widthplot <- get_meanplot(yvar = "CIwidth", cw_to = NA, percent_labels = FALSE, ylabel = "Mean CI Width")
ggsave(paste0(code_dir, "plots/FigureS2_CIwidth_", date, ".pdf"), widthplot, width = 14, height = 5)
