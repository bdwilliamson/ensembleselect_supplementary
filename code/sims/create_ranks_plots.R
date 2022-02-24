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
parser$add_argument("--sim-name", default = "ranks-binomial-linear-normal",
                    help = "the name of the simulation")
parser$add_argument("--b", type = "double", default = 20,
                    help = "number of bootstrap replicates")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "total number of Monte-Carlo replicates")
args <- parser$parse_args()

# -------------------------------------------------------------------------------------
# read in compiled results for the given sim name, b, m (note that we're compiling all)
# -------------------------------------------------------------------------------------
all_files <- list.files(compiled_output_dir)
if (args$sim_name == "ranks-binomial-linear-normal") {
  files_to_load <- paste0(compiled_output_dir, all_files[grepl(args$sim_name, all_files)])
} else {
  files_to_load <- paste0(compiled_output_dir, all_files[grepl(args$sim_name, all_files) & grepl(args$b, all_files)])
}

all_results <- lapply(as.list(files_to_load), function(x) {
  if(all(is.na(x))) {
    NA
  } else {
    readRDS(x)
  }})
non_na_results <- unlist(lapply(all_results, function(x) !all(is.na(x))))
output_tib <- as_tibble(rbindlist(all_results[non_na_results]))

# -------------------------------------------------------------------------------------
# average results within a given estimator, extra layer, missing percentage
# -------------------------------------------------------------------------------------
selection_types <- c("sl_none", "sl_rank", "sl_screen")
performance_tib <- output_tib %>%
  group_by(selection_type, extra_layer, n, p) %>%
  mutate(set_size = get_selected_set_size(selected_vars),
         sensitivity_init = get_sensitivity(selected_vars, true_active_set = 1:6),
         specificity_init = get_specificity(selected_vars, true_active_set = 1:6, p = p),
         mse_init = n * (perf - true_perf) ^ 2) %>%
  summarize(set_size = mean(set_size), mse = mean(mse_init),
            mse_mcse = sqrt(sum((mse_init - mse)^2)) / args$nreps_total,
            sensitivity = mean(sensitivity_init), specificity = mean(specificity_init),
            sensitivity_mcse = sqrt(sum((sensitivity_init - sensitivity)^2)) / args$nreps_total,
            specificity_mcse = sqrt(sum((specificity_init - specificity)^2)) / args$nreps_total,
            .groups = "drop") %>%
  filter(selection_type %in% selection_types)
# -------------------------------------------------------------------------------------
# create plots/tables
# -------------------------------------------------------------------------------------
est_levels <- unique(paste(performance_tib$selection_type, performance_tib$extra_layer, sep = "; "))
est_labels <- make_est_levels(performance_tib)
num_ns <- length(unique(performance_tib$n))
point_vals <- (0:25)[1:length(est_levels)]
point_size <- 2
axis_text_size <- 20
lgnd_text_size <- 14
title_text_size <- 30
label_size <- title_text_size
fig_width <- 43
fig_height <- 25
dodge_width <- 50
legend_ords_init <- c(1, 2, 5, 3, 7, 6, 4, 8)
legend_ords <- rep(legend_ords_init, each = num_ns)
point_val_vec <- rep(point_vals, each = num_ns)
if (grepl("weaklinear", args$sim_name)) {
  minus_vec <- (5 * num_ns + 1):(5 * num_ns + 2)
} else if (grepl("nonlinear", args$sim_name)) {
  minus_vec <- (5 * num_ns + 1):(5 * num_ns + 2)
} else {
  minus_vec <- 9 * num_ns + 1
}
legend_ords <- legend_ords[-minus_vec]
point_val_vec <- point_val_vec[-minus_vec]
performance_tib2 <- performance_tib %>%
  ungroup() %>%
  mutate(nice_est = factor(paste(selection_type, extra_layer, sep = "; "),
                           levels = est_levels, labels = est_labels),
         legend_ord = legend_ords,
         point_val = point_val_vec)


# ------------------------
# plots of MSE
# ------------------------
these_point_vals <- unique(performance_tib2$point_val)
mse_plot <- performance_tib2 %>%
  ggplot(aes(x = factor(n), y = mse,
             shape = reorder(nice_est, legend_ord))) +
  geom_point(position = position_dodge(0.5),
             size = 2 * point_size) +
  scale_shape_manual(values = these_point_vals) +
  geom_errorbar(aes(ymin = mse - 1.96 * mse_mcse, ymax = mse + 1.96 * mse_mcse),
                position = position_dodge(0.5),
                size = point_size / 2) +
  ylab(expression(paste(n, " x empirical ", MSE[n]))) +
  xlab("n") +
  labs(shape = "Selection type; extra layer") +
  guides(shape = guide_legend(nrow = 2)) +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))

# a second plot, zoomed in on the small stuff
mse_filter <- 20
mse_plot_filtered <- performance_tib2 %>%
  ggplot(aes(x = factor(n), y = mse,
             shape = reorder(nice_est, legend_ord))) +
  geom_point(position = position_dodge(0.5),
             size = 2 * point_size) +
  scale_shape_manual(values = these_point_vals) +
  geom_errorbar(aes(ymin = mse - 1.96 * mse_mcse, ymax = mse + 1.96 * mse_mcse),
                position = position_dodge(0.5),
                size = point_size / 2) +
  ggtitle(paste0("Excluding MSE > ", mse_filter)) +
  ylab(expression(paste(n, " x empirical ", MSE[n]))) +
  ylim(c(0, 20)) +
  xlab("n") +
  labs(shape = "Selection type; extra layer") +
  scale_shape_manual(values = these_point_vals) +
  guides(shape = guide_legend(nrow = 2)) +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))

mse_plt_row <- plot_grid(mse_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                         mse_plot_filtered + theme(legend.position = "none"),
                         align = "vh", axis = "t", labels = "AUTO", label_size = label_size,
                         nrow = 1)
mse_lgnd <- get_legend(mse_plot +
                         theme(legend.spacing.x = unit(1.5, "cm"))
)
mse_title <- ggdraw() + draw_label(bquote(bold(paste("Empirical MSE scaled by ", n))), size = title_text_size)
mse_plt <- plot_grid(mse_title,
                     mse_plt_row,
                     mse_lgnd,
                     ncol = 1, nrow = 3,
                     # axis = "t", align = "vh",
                     rel_heights = c(.075, 1, .1))
ggsave(filename = paste0(plots_dir, args$sim_name, "_mse_plots.png"),
       plot = mse_plt,
       device = "png",
       width = fig_width, height = fig_height, units = "cm", dpi = 300)


# --------------------------------
# plots of sensitivity/specificity
# --------------------------------
y_lim <- c(min(min(performance_tib$specificity - 1.96 * performance_tib$specificity_mcse),
               min(performance_tib$sensitivity - 1.96 * performance_tib$sensitivity_mcse)),
           max(max(performance_tib$specificity + 1.96 * performance_tib$specificity_mcse),
               max(performance_tib$sensitivity + 1.96 * performance_tib$sensitivity_mcse)))
these_point_vals <- unique(performance_tib2$point_val)

sens_plot <- performance_tib2 %>%
  ggplot(aes(x = factor(n), y = sensitivity,
             group = factor(paste(selection_type, extra_layer, n, sep = "; ")),
             shape = reorder(nice_est, legend_ord))) +
  geom_point(position = position_dodge(0.5),
             size = 2 * point_size) +
  scale_shape_manual(values = these_point_vals) +
  geom_errorbar(aes(ymin = sensitivity - 1.96 * sensitivity_mcse,
                    ymax = sensitivity + 1.96 * sensitivity_mcse),
                position = position_dodge(0.5),
                size = point_size / 2) +
  ylab("Empirical sensitivity") +
  xlab("n") +
  labs(shape = "Selection type; extra layer") +
  guides(shape = guide_legend(nrow = 2)) +
  theme(legend.direction = "horizontal", legend.position = "bottom",
        legend.text = element_text(size = lgnd_text_size), legend.title = element_text(size = lgnd_text_size),
        plot.margin = unit(c(0, 2, 0, 0), "cm"))

spec_plot <- performance_tib2 %>%
  ggplot(aes(x = factor(n), y = specificity,
             group = factor(paste(selection_type, extra_layer, n, sep = "; ")),
             shape = reorder(nice_est, legend_ord))) +
  geom_point(position = position_dodge(0.5),
             size = 2 * point_size) +
  scale_shape_manual(values = these_point_vals) +
  geom_errorbar(aes(ymin = specificity - 1.96 * specificity_mcse,
                    ymax = specificity + 1.96 * specificity_mcse),
                position = position_dodge(0.5),
                size = point_size / 2) +
  ylab("Empirical specificity") +
  xlab("n") +
  labs(shape = "Selection type; extra layer") +
  guides(shape = guide_legend(nrow = 2)) +
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm"))

sens_spec_rows <- plot_grid(sens_plot + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
                            spec_plot + theme(legend.position = "none"))
sens_spec_lgnd <- get_legend(sens_plot +
                               theme(legend.spacing.x = unit(1.5, "cm")))
sens_title <- ggdraw() + draw_label(bquote(bold(paste("Empirical sensitivity"))), size = title_text_size, hjust = 0.35)
spec_title <- ggdraw() + draw_label(bquote(bold(paste("Empirical specificity"))), size = title_text_size, hjust = 0.35)
# create and save plots, one for each p (args$sim_name handles linear predictor and x distribution)
sens_spec_plt <- plot_grid(plot_grid(sens_title, spec_title),
                           sens_spec_rows,
                           sens_spec_lgnd,
                           ncol = 1, nrow = 3,
                           # axis = "t", align = "vh",
                           rel_heights = c(.05, 1, .1))
ggsave(filename = paste0(plots_dir, args$sim_name, "_sens-spec_plots.png"),
       plot = sens_spec_plt,
       device = "png",
       width = fig_width, height = fig_height + 2, units = "cm", dpi = 300)

# table with average estimated AUC and true AUC for each dgm and estimator
aucs <- output_tib %>%
  filter(selection_type %in% selection_types) %>% 
  group_by(selection_type, extra_layer, n, p) %>%
  summarize(est_auc = mean(perf), true_auc = mean(true_perf), .groups = "drop") %>%
  mutate(round_auc = round(est_auc, 3), round_true_auc = round(true_auc, 3)) %>%
  filter(n == 200 | n == 500 | n == 4000) %>%
  ungroup() %>%
  mutate(nice_est = factor(paste(selection_type, extra_layer, sep = "; "),
                           levels = est_levels, labels = est_labels))

print_tab <- aucs %>%
  arrange(p) %>%
  select(n, p, nice_est, est_auc, round_auc, round_true_auc) %>%
  rename(Estimator = nice_est) %>%
  select(-est_auc, -round_true_auc) %>%
  rename(`Estimated AUC` = round_auc) %>%
  pivot_wider(names_from = n, values_from = `Estimated AUC`) %>%
  rename(`Estimated AUC (n = 200)` = `200`, `Estimated AUC (n = 500)` = `500`,
         `Estimated AUC (n = 4000)` = `4000`)
get_setting <- function(sim_name) {
  if (grepl("linear-normal", sim_name) & !grepl("non", sim_name) & !grepl("weak", sim_name)) {
    return("A")
  } else if (grepl("nonlinear-normal", sim_name)) {
    return("B")
  } else if (grepl("linear-nonnormal", sim_name) & !grepl("nonlinear", sim_name)) {
    return("C")
  } else if (grepl("weak", sim_name)){
    return("E")
  } else {
    return("D")
  }
}
knitr::kable(print_tab,
             label = paste0(gsub("-nested", "", gsub("binomial-", "", args$sim_name)), "-aucs"),
             caption = paste0("AUC (estimated on test data) of the estimated panel (estimated on training data) averaged over Monte Carlo replications for each combination of $n \\in \\{200, 500, 4000\\}$ and estimator in setting ", get_setting(args$sim_name),". The true AUC of the optimal panel is ", round(mean(aucs$true_auc), 3), "."),
             format = "latex", booktabs = TRUE, linesep = "") %>%
  cat(., file = paste0(plots_dir, args$sim_name, "_auc_tib.tex"))
