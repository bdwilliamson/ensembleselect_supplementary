# make plots based on the simulation output

# ------------------------------------------------------------------------------
# load required libraries, set up code and output directories
# source required functions, set up args
# ------------------------------------------------------------------------------
library("ggplot2")
library("dplyr")
library("tidyr")
library("cowplot")
theme_set(theme_cowplot())
library("argparse")
library("data.table")

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    code_dir <- "sim_code/"
    raw_output_dir <- "sim_output/"
    compiled_output_dir <- "sim_output/compiled_results/"
    plots_dir <- "plots/"
} else {
    code_dir <- "./"
    raw_output_dir <- "../sim_output/"
    compiled_output_dir <- "../sim_output/compiled_results/"
    plots_dir <- "../plots/"
}
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
}
source(paste0(code_dir, "utils.R"))

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "binomial-probit-linear-normal-nested",
                    help = "the name of the simulation")
parser$add_argument("--b", type = "double", default = 100,
                    help = "number of bootstrap replicates")
parser$add_argument("--m", type = "double", default = 1,
                    help = "number of MI replicates")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                     help = "total number of Monte-Carlo replicates")
parser$add_argument("--all-miss", type = "double", default = 1,
                    help = "should we include all missing data proportions in plot?")
args <- parser$parse_args()

if (!grepl("probit", args$sim_name)) {
    plots_dir <- paste0(plots_dir, "logistic/")
    if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
    }
}

read_func <- function(x, type = "select") {
    if (all(is.na(x))) {
        NA
    } else {
        these_results <- readRDS(x)
        if (!any(grepl("est", names(x))) & type == "vim") {
            these_results %>% mutate(est = rank)
        } else {
            these_results
        }
    }
}
# ------------------------------------------------------------------------------
# read in compiled results for the given sim name, b, m
# (note that we're compiling all)
# ------------------------------------------------------------------------------
all_files <- list.files(compiled_output_dir)
files_to_load <- paste0(compiled_output_dir,
                        all_files[grepl(args$sim_name, all_files) &
                                      grepl("select", all_files)])
all_results <- lapply(
    as.list(files_to_load),
    read_func
)
non_na_results <- unlist(lapply(all_results, function(x) !all(is.na(x))))
output_tib <- as_tibble(rbindlist(all_results[non_na_results]))

# ------------------------------------------------------------------------------
# average results within a given estimator, extra layer, missing percentage
# ------------------------------------------------------------------------------
true_perf <- output_tib %>%
    filter(missing_perc == 0) %>%
    group_by(p) %>%
    summarize(true_perf = mean(true_perf))
output_tib_fixed_trueperf <- output_tib %>%
    select(-true_perf) %>%
    dplyr::left_join(true_perf, by = "p")
if (!grepl("probit", args$sim_name)) {
    output_tib_fixed_trueperf <- output_tib_fixed_trueperf %>%
        filter(missing_perc %in% c(0, 0.1, 0.3)) %>%
        mutate(sensitivity = get_sensitivity(selected_vars, true_active_set = 1:6),
               specificity = get_specificity(selected_vars, true_active_set = 1:6, p = p))
}
performance_tib_init <- output_tib_fixed_trueperf %>%
  rename(sensitivity_init = sensitivity, specificity_init = specificity) %>%
  mutate(set_size_init = get_selected_set_size(selected_vars),
         bias_init = sqrt(n) * (perf - true_perf),
         mse_init = n * (perf - true_perf) ^ 2,
         est_fct = factor(case_when(
           estimator_type == "lasso" & extra_layer == "none" ~ 1,
           estimator_type == "SL" & extra_layer == "none" ~ 2,
           estimator_type == "lasso" & extra_layer == "SS" ~ 3,
           estimator_type == "SL" & extra_layer == "SS" ~ 4,
           estimator_type == "lasso" & (extra_layer == "KF" | extra_layer == "knockoffs") ~ 5,
           estimator_type == "base-SL" & extra_layer == "none" ~ 6,
           estimator_type == "SPVIM-gFWER" & extra_layer == "none" ~ 7,
           estimator_type == "SPVIM-PFP" & extra_layer == "none" ~ 8,
           estimator_type == "SPVIM-FDR" & extra_layer == "none" ~ 9,
         ), levels = 1:9, labels = c("lasso", "SL", "lasso + SS",
                                     "SL + SS", "lasso + KF", "SL (benchmark)",
                                     "SPVIM + gFWER", "SPVIM + PFP", "SPVIM + FDR"),
         ordered = TRUE))
performance_tib <- performance_tib_init %>%
    group_by(est_fct, missing_perc, n, p) %>%
    summarize(set_size = mean(set_size_init), mse = mean(mse_init),
              mse_mcse = sqrt(sum((mse_init - mse)^2)) / args$nreps_total,
              bias = mean(bias_init), bias_mcse = sqrt(sum((bias_init - bias) ^ 2)) / args$nreps_total,
              variance = var(perf),
              sensitivity = mean(sensitivity_init),
              specificity = mean(specificity_init),
              sensitivity_mcse = sqrt(sum((sensitivity_init - sensitivity)^2)) /
                  args$nreps_total,
              specificity_mcse = sqrt(sum((specificity_init - specificity)^2)) /
                  args$nreps_total,
            .groups = "drop") %>%
    mutate(missing = paste0(round(missing_perc * 100), "%"), variance = n * variance)
plot_tib <- performance_tib %>%
  # filter(!grepl("SPVIM", est_fct)) %>%
  mutate(n_fct = factor(n))
# ------------------------------------------------------------------------------
# create plots/tables
# ------------------------------------------------------------------------------
miss_percs <- unique(plot_tib$missing_perc)
point_size <- 1.5
axis_text_size <- 8.5
lgnd_text_size <- 8.5
title_text_size <- 10
label_size <- title_text_size
fig_width <- 8.5
fig_height <- 4
shapes <- c(15:18, 3, 4, 8, 5, 6)
unique_p <- unique(plot_tib$p)
image_types <- c("png", "tiff")

# determine which to plot; for main MS figures, don't want to show all missing proportions
# or all SPVIM combos
if (args$all_miss == 1) {
  miss_vec <- miss_percs
} else {
  miss_vec <- miss_percs[c(1, 3)]
}

# --------------------------
# plots of bias and variance
# --------------------------
min_bias <- min(min(plot_tib %>% filter(!grepl("FDR", est_fct) & !grepl("PFP", est_fct)) %>%
                      pull(bias)), 0)
max_bias <- max(max(plot_tib %>% pull(bias)), 0)
bias_lim <- max(abs(min_bias), abs(max_bias))
bias_plot <- plot_tib %>%
  filter(missing_perc %in% miss_vec) %>%
  ggplot(aes(x = n_fct, y = bias, shape = est_fct)) +
  geom_point(position = position_dodge(0.8), size = point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab(expression(paste(sqrt(n), " x empirical ", bias[n]))) +
  xlab("n") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
  scale_shape_manual(values = shapes) +
  labs(shape = "Estimator") +
  guides(shape = guide_legend(nrow = 1)) +
  ylim(c((-1) * bias_lim - 0.05, bias_lim + 0.05)) +
  facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))

variance_plot <- plot_tib %>%
  filter(missing_perc %in% miss_vec) %>%
  ggplot(aes(x = n_fct, y = variance, shape = est_fct)) +
  geom_point(position = position_dodge(0.8), size = point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab(expression(paste(n, " x empirical ", variance[n]))) +
  xlab("n") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
  scale_shape_manual(values = shapes) +
  facet_grid(~ missing_perc + p) +
  labs(shape = "Estimator") +
  guides(shape = guide_legend(nrow = 1)) +
  facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))
shared_legend <- get_legend(bias_plot +
                              theme(legend.spacing.x = unit(0.1, "in")) +
                              guides(shape = guide_legend(nrow = 2))
)
bias_var_title <- ggdraw() + draw_label(bquote("Prediction performance: AUC"), size = title_text_size)
bias_var_plot <- plot_grid(bias_var_title,
                           plot_grid(
                             bias_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                             variance_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                             nrow = 1, ncol = 2, labels = "AUTO"
                           ),
                           shared_legend, ncol = 1, nrow = 3,
                           rel_heights = c(.075, 1, .1))

for (image_type in image_types) {
  ggsave(filename = paste0(plots_dir, args$sim_name, "_pred-perf.", image_type),
         plot = bias_var_plot, width = fig_width, height = fig_height,
         device = image_type, units = "in", dpi = 300)
}

# separate plots for the supplement
if (!args$all_miss) {
  min_bias <- min(min(plot_tib %>% filter(!grepl("FDR", est_fct) & !grepl("PFP", est_fct)) %>%
                        pull(bias)), 0)
  max_bias <- max(max(plot_tib %>% pull(bias)), 0)
  bias_lim <- max(abs(min_bias), abs(max_bias))
  bias_plot <- plot_tib %>%
    filter(!(missing_perc %in% miss_vec)) %>%
    ggplot(aes(x = n_fct, y = bias, shape = est_fct)) +
    geom_point(position = position_dodge(0.8), size = point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylab(expression(paste(sqrt(n), " x empirical ", bias[n]))) +
    xlab("n") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
    scale_shape_manual(values = shapes) +
    labs(shape = "Estimator") +
    guides(shape = guide_legend(nrow = 1)) +
    ylim(c((-1) * bias_lim - 0.05, bias_lim + 0.05)) +
    facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
    theme(legend.direction = "horizontal", legend.position = "bottom",
          title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
          axis.text = element_text(size = axis_text_size),
          legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
          strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
          plot.margin = unit(c(0, 2, 0, 0), "cm"))

  variance_plot <- plot_tib %>%
    filter(!(missing_perc %in% miss_vec)) %>%
    ggplot(aes(x = n_fct, y = variance, shape = est_fct)) +
    geom_point(position = position_dodge(0.8), size = point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ylab(expression(paste(n, " x empirical ", variance[n]))) +
    xlab("n") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
    scale_shape_manual(values = shapes) +
    facet_grid(~ missing_perc + p) +
    labs(shape = "Estimator") +
    guides(shape = guide_legend(nrow = 1)) +
    facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
    theme(legend.direction = "horizontal", legend.position = "bottom",
          title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
          axis.text = element_text(size = axis_text_size),
          legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
          strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
          plot.margin = unit(c(0, 2, 0, 0), "cm"))
  shared_legend <- get_legend(bias_plot +
                                theme(legend.spacing.x = unit(0.1, "in")) +
                                guides(shape = guide_legend(nrow = 2))
  )
  bias_var_title <- ggdraw() + draw_label(bquote("Prediction performance: AUC"), size = title_text_size)
  bias_var_plot <- plot_grid(bias_var_title,
                             plot_grid(
                               bias_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                               variance_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                               nrow = 1, ncol = 2, labels = "AUTO"
                             ),
                             shared_legend, ncol = 1, nrow = 3,
                             rel_heights = c(.075, 1, .1))

  for (image_type in image_types) {
    ggsave(filename = paste0(plots_dir, args$sim_name, "_pred-perf_supp.", image_type),
           plot = bias_var_plot, width = fig_width, height = fig_height,
           device = image_type, units = "in", dpi = 300)
  }
}

# --------------------------------
# plots of sensitivity/specificity
# --------------------------------
sens_plot <- plot_tib %>%
  filter(missing_perc %in% miss_vec) %>%
  ggplot(aes(x = n_fct, y = sensitivity, shape = est_fct)) +
  geom_point(position = position_dodge(0.8), size = point_size) +
  ylab("Empirical sensitivity") +
  xlab("n") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
  scale_shape_manual(values = shapes) +
  labs(shape = "Estimator") +
  guides(shape = guide_legend(nrow = 1)) +
  facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))

spec_plot <- plot_tib %>%
  filter(missing_perc %in% miss_vec) %>%
  ggplot(aes(x = n_fct, y = specificity, shape = est_fct)) +
  geom_point(position = position_dodge(0.8), size = point_size) +
  ylab("Empirical specificity") +
  xlab("n") +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
  scale_shape_manual(values = shapes) +
  facet_grid(~ missing_perc + p) +
  labs(shape = "Estimator") +
  guides(shape = guide_legend(nrow = 1)) +
  facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))
select_legend <- get_legend(sens_plot +
                              theme(legend.spacing.x = unit(0.1, "in")) +
                              guides(shape = guide_legend(nrow = 2))
)
select_title <- ggdraw() + draw_label(bquote("Variable selection performance"), size = title_text_size)
select_plot <- plot_grid(select_title,
                         plot_grid(
                           sens_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                           spec_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                           nrow = 1, ncol = 2, labels = "AUTO"
                         ),
                         select_legend, ncol = 1, nrow = 3,
                         rel_heights = c(.075, 1, .1))

for (image_type in image_types) {
  ggsave(filename = paste0(plots_dir, args$sim_name, "_select-perf.", image_type),
         plot = select_plot, width = fig_width, height = fig_height,
         device = image_type, units = "in", dpi = 300)
}

# separate plots for the supplement
if (!args$all_miss) {
  sens_plot <- plot_tib %>%
    filter(!(missing_perc %in% miss_vec)) %>%
    ggplot(aes(x = n_fct, y = sensitivity, shape = est_fct)) +
    geom_point(position = position_dodge(0.8), size = point_size) +
    ylab("Empirical sensitivity") +
    xlab("n") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
    scale_shape_manual(values = shapes) +
    labs(shape = "Estimator") +
    guides(shape = guide_legend(nrow = 1)) +
    facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
    theme(legend.direction = "horizontal", legend.position = "bottom",
          title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
          axis.text = element_text(size = axis_text_size),
          legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
          strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
          plot.margin = unit(c(0, 2, 0, 0), "cm"))

  spec_plot <- plot_tib %>%
    filter(!(missing_perc %in% miss_vec)) %>%
    ggplot(aes(x = n_fct, y = specificity, shape = est_fct)) +
    geom_point(position = position_dodge(0.8), size = point_size) +
    ylab("Empirical specificity") +
    xlab("n") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), color = "gray85") +
    scale_shape_manual(values = shapes) +
    facet_grid(~ missing_perc + p) +
    labs(shape = "Estimator") +
    guides(shape = guide_legend(nrow = 1)) +
    facet_grid(cols = vars(p), rows = vars(missing), labeller = label_both, scales = "free") +
    theme(legend.direction = "horizontal", legend.position = "bottom",
          title = element_text(size = title_text_size), text = element_text(size = axis_text_size),
          axis.text = element_text(size = axis_text_size),
          legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
          strip.background = element_blank(), panel.grid.major.y = element_line(color = "gray85"),
          plot.margin = unit(c(0, 2, 0, 0), "cm"))
  select_legend <- get_legend(sens_plot +
                                theme(legend.spacing.x = unit(0.1, "in")) +
                                guides(shape = guide_legend(nrow = 2))
  )
  select_title <- ggdraw() + draw_label(bquote("Variable selection performance"), size = title_text_size)
  select_plot <- plot_grid(select_title,
                           plot_grid(
                             sens_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                             spec_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                             nrow = 1, ncol = 2, labels = "AUTO"
                           ),
                           select_legend, ncol = 1, nrow = 3,
                           rel_heights = c(.075, 1, .1))

  for (image_type in image_types) {
    ggsave(filename = paste0(plots_dir, args$sim_name, "_select-perf_supp.", image_type),
           plot = select_plot, width = fig_width, height = fig_height,
           device = image_type, units = "in", dpi = 300)
  }
}

# table with average estimated AUC and true AUC for each dgm and estimator
aucs <- output_tib %>%
    group_by(estimator_type, extra_layer, missing_type, missing_perc, n, p) %>%
    summarize(est_auc = mean(perf), true_auc = mean(true_perf), .groups = "drop") %>%
    mutate(round_auc = round(est_auc, 3), round_true_auc = round(true_auc, 3)) %>%
    filter(n == 500 | n == 3000, missing_perc == 0 | missing_perc == 0.5) %>%
    ungroup()
readr::write_csv(aucs %>% filter(p == unique_p[1]) %>% arrange(missing_perc) %>% select(n, missing_perc, estimator_type, extra_layer, est_auc, round_auc, round_true_auc), path = paste0(plots_dir, args$sim_name, "_auc_tib_p_30.csv"))
readr::write_csv(aucs %>% filter(p == unique_p[2]) %>% arrange(missing_perc) %>% select(n, missing_perc, estimator_type, extra_layer, est_auc, round_auc, round_true_auc), path = paste0(plots_dir, args$sim_name, "_auc_tib_p_500.csv"))

print_tab <- aucs %>%
    filter(!(estimator_type == "SL" & extra_layer == "knockoffs"),
           extra_layer != "knockoffs", extra_layer != "SS") %>%
    arrange(p, missing_perc) %>%
    select(n, p, missing_perc, estimator_type, extra_layer, est_auc, round_auc, round_true_auc) %>%
    mutate(Estimator = paste(estimator_type, extra_layer, sep = " + ")) %>%
    select(-estimator_type, -extra_layer, -est_auc, -round_true_auc) %>%
    rename(`% missing` = missing_perc, `Estimated AUC` = round_auc) %>%
    pivot_wider(names_from = n, values_from = `Estimated AUC`) %>%
    rename(`Estimated AUC (n = 500)` = `500`, `Estimated AUC (n = 3000)` = `3000`)
get_setting <- function(sim_name) {
    if (grepl("linear-normal", sim_name) & !grepl("non", sim_name)) {
        return("A")
    } else if (grepl("nonlinear-normal", sim_name)) {
        return("B")
    } else if (grepl("linear-nonnormal", sim_name) & !grepl("nonlinear", sim_name)) {
        return("C")
    } else {
        return("D")
    }
}
knitr::kable(print_tab,
             label = paste0(gsub("-nested", "", gsub("binomial-", "", args$sim_name)), "-aucs"),
    caption = paste0("Setting ", get_setting(args$sim_name), ": AUC (estimated on test data) of the estimated panel (estimated on training data) averaged over Monte Carlo replications for each combination of $n \\in \\{500, 4000\\}$, $p$, percentage of missing data $\\in \\{0, 0.5\\}$, and estimator. The true AUC of the optimal panel is ", round(mean(aucs$true_auc), 3), "."),
    format = "latex", booktabs = TRUE, linesep = "") %>%
    cat(., file = paste0(plots_dir, args$sim_name, "_auc_tib.tex"))
