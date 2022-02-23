# create ROC curves for SL and SS

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
library("flevr")

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

# ------------------------------------------------------------------------------
# ROC curves for different rank cutoffs for SL, any SS, and SPVIM
# ------------------------------------------------------------------------------
# read in the variable importance results; good for debugging
estimators <- c("lasso_SS", "SL_none", "SL_SS", "SPVIM_none")
vim_files_to_load <- paste0(compiled_output_dir,
                            all_files[grepl(args$sim_name, all_files) &
                                        grepl("vim", all_files) &
                                        apply(do.call(rbind, lapply(as.list(estimators),
                                                                    function(x) grepl(x, all_files) & !grepl("base", all_files))),
                                              2, any)])
all_vim_results <- lapply(as.list(vim_files_to_load), read_func)
non_na_vim_results <- unlist(lapply(all_vim_results, function(x) !all(is.na(x))))
vim_output_tib <- rbindlist(all_vim_results[non_na_vim_results], fill = TRUE)

# get sensitivity and specificity
rank_thresholds_30 <- seq(1, 25, by = 1) # for SL
rank_thresholds_500 <- seq(5, 400, by = 25) # for SL
pi_thresholds <- seq(0.7, 1, by = 0.01) # for SS
ns <- c(200, 500, 1500, 3000)
all_specificities_30 <- unlist(lapply(as.list(ns), function(n) set_specificity(n, ns = ns, p = 30)))
all_specificities_500 <- unlist(lapply(as.list(ns), function(n) set_specificity(n, ns = ns, p = 500)))
initial_gfwer_ks_30 <- ceiling((1 - all_specificities_30) * (30 - (7 - 1)))
initial_gfwer_ks_500 <- ceiling((1 - all_specificities_500) * (500 - (7 - 1)))
initial_pfp_qs_30 <- initial_gfwer_ks_30 / ((30 - 6) / 30 * (sqrt(ns) / sqrt(200)) + 6)
initial_pfp_qs_500 <- initial_gfwer_ks_500 / ((500 - 6) / 500 * (sqrt(ns) / sqrt(200)) + 6)
# at 30, min was 3 and max was 5
gfwer_ks_30 <- seq(1, 25, by = 1)
# at 500, min was 11 and max was 45
gfwer_ks_500 <- seq(10, 250, by = 10)
pfp_qs <- seq(0.2, 1, by = 0.05)
alpha <- 0.05
true_active_set <- 1:6
true_null_set_30 <- 7:30
true_null_set_500 <- 7:500
sl_30 <- as_tibble(vim_output_tib[(estimator_type == "SL" & extra_layer != "SS") & (p == 30), ])
sl_500 <- as_tibble(vim_output_tib[(estimator_type == "SL"  & extra_layer != "SS") & (p == 500), ])
ss_30 <- as_tibble(vim_output_tib[(extra_layer == "SS") & (p == 30), ])
ss_500 <- as_tibble(vim_output_tib[(extra_layer == "SS") & (p == 500), ])
spvim_30 <- as_tibble(vim_output_tib[(estimator_type == "SPVIM") & (p == 30), ])
spvim_500 <- as_tibble(vim_output_tib[(estimator_type == "SPVIM") & (p == 500), ])

# sensitivity and specificity for the base Super Learner
sl_sens_spec_30 <- rank_thresholds_30 %>%
  purrr::map_dfr(function(x) {
    sl_30 %>%
      mutate(threshold = x,
             tpr_contribution = ifelse(rank <= threshold & feature_num %in% true_active_set, 1, 0),
             fpr_contribution = ifelse(rank <= threshold & feature_num %in% 7:30, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:30), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop")
  })
sl_sens_spec_500 <- rank_thresholds_500 %>%
  purrr::map_dfr(function(x) {
    sl_500 %>%
      mutate(threshold = x,
             tpr_contribution = ifelse(rank <= threshold & feature_num %in% true_active_set, 1, 0),
             fpr_contribution = ifelse(rank <= threshold & feature_num %in% 7:500, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:500), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop")
  })
# sensitivity and specificity for stability selection-based methods (lasso, SL)
ss_sens_spec_30 <- pi_thresholds %>%
  purrr::map_dfr(function(x) {
    ss_30 %>%
      mutate(threshold = x,
             tpr_contribution = ifelse(est >= threshold & feature_num %in% true_active_set, 1, 0),
             fpr_contribution = ifelse(est >= threshold & feature_num %in% 7:30, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:30), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop")
  })
ss_sens_spec_500 <- pi_thresholds %>%
  purrr::map_dfr(function(x) {
    ss_500 %>%
      mutate(threshold = x,
             tpr_contribution = ifelse(est >= threshold & feature_num %in% true_active_set, 1, 0),
             fpr_contribution = ifelse(est >= threshold & feature_num %in% 7:500, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:500), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop")
  })
# sensitivity and specificity for SPVIM
# first, compute adjusted p-values for each
spvim <- bind_rows(spvim_30, spvim_500)
spvim$adj_p <- unlist(spvim %>%
                        group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p, missing_perc, n, mc_id) %>%
                        group_map(~ flevr::get_base_set(p_values = .x$p_value, test_statistics = .x$p_value,
                                                        alpha = alpha, method = "Holm")$p_values))
spvim$initially_rejected <- spvim$adj_p <= alpha
spvim_30 <- spvim %>% filter(p == 30)
spvim_500 <- spvim %>% filter(p == 500)
# gFWER(k)
spvim_gfwer_sens_spec_30 <- gfwer_ks_30 %>%
  purrr::map_dfr(function(x) {
    spvim_30 %>%
      mutate(threshold = x, augmented_set = unlist(
               spvim_30 %>%
                 group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p, missing_perc, n, mc_id) %>%
                 group_map(~ flevr::get_augmented_set(p_values = .x$p_value, num_rejected = sum(.x$initially_rejected),
                                                      alpha = alpha, quantity = "gFWER", k = x)$set)
             ), selected = initially_rejected | augmented_set,
             tpr_contribution = ifelse(selected & feature_num %in% true_active_set, 1, 0),
             fpr_contribution = ifelse(selected & feature_num %in% 7:30, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:30), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop") %>%
      mutate(mn_sens = pmin(mn_sens, 1), mn_spec = pmin(mn_spec, 1))
  })
spvim_gfwer_sens_spec_500 <- gfwer_ks_500 %>%
  purrr::map_dfr(function(x) {
    spvim_500 %>%
      mutate(threshold = x, augmented_set = unlist(
        spvim_500 %>%
          group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p, missing_perc, n, mc_id) %>%
          group_map(~ flevr::get_augmented_set(p_values = .x$p_value, num_rejected = sum(.x$initially_rejected),
                                               alpha = alpha, quantity = "gFWER", k = x)$set)
      ), selected = initially_rejected | augmented_set,
      tpr_contribution = ifelse(selected & feature_num %in% true_active_set, 1, 0),
      fpr_contribution = ifelse(selected & feature_num %in% 7:500, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:500), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop") %>%
      mutate(mn_sens = pmin(mn_sens, 1), mn_spec = pmin(mn_spec, 1))
  })
# PFP(q)
spvim_pfp_sens_spec_30 <- pfp_qs %>%
  purrr::map_dfr(function(x) {
    spvim_30 %>%
      mutate(threshold = x, augmented_set = unlist(
        spvim_30 %>%
          group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p, missing_perc, n, mc_id) %>%
          group_map(~ flevr::get_augmented_set(p_values = .x$p_value, num_rejected = sum(.x$initially_rejected),
                                               alpha = alpha, quantity = "PFP", q = x)$set)
      ), selected = initially_rejected | augmented_set,
      tpr_contribution = ifelse(selected & feature_num %in% true_active_set, 1, 0),
      fpr_contribution = ifelse(selected & feature_num %in% 7:30, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:30), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop") %>%
      mutate(mn_sens = pmin(mn_sens, 1), mn_spec = pmin(mn_spec, 1))
  })
spvim_pfp_sens_spec_500 <- pfp_qs %>%
  purrr::map_dfr(function(x) {
    spvim_500 %>%
      mutate(threshold = x, augmented_set = unlist(
        spvim_500 %>%
          group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p, missing_perc, n, mc_id) %>%
          group_map(~ flevr::get_augmented_set(p_values = .x$p_value, num_rejected = sum(.x$initially_rejected),
                                               alpha = alpha, quantity = "PFP", q = x)$set)
      ), selected = initially_rejected | augmented_set,
      tpr_contribution = ifelse(selected & feature_num %in% true_active_set, 1, 0),
      fpr_contribution = ifelse(selected & feature_num %in% 7:500, 1, 0)) %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold, mc_id) %>%
      summarize(sens = sum(tpr_contribution) / length(true_active_set),
                spec = 1 - sum(fpr_contribution) / length(7:500), .groups = "drop") %>%
      group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer,
               p, missing_perc, n, threshold) %>%
      summarize(mn_sens = mean(sens), mn_spec = mean(spec), .groups = "drop") %>%
      mutate(mn_sens = pmin(mn_sens, 1), mn_spec = pmin(mn_spec, 1))
  })

# pull it all together
est_sens_spec <- sl_sens_spec_30 %>%
  bind_rows(sl_sens_spec_500) %>%
  bind_rows(ss_sens_spec_30) %>%
  bind_rows(ss_sens_spec_500) %>%
  bind_rows(spvim_gfwer_sens_spec_30 %>% mutate(extra_layer = "gFWER")) %>%
  bind_rows(spvim_gfwer_sens_spec_500 %>% mutate(extra_layer = "gFWER")) %>%
  bind_rows(spvim_pfp_sens_spec_30 %>% mutate(extra_layer = "PFP")) %>%
  bind_rows(spvim_pfp_sens_spec_500 %>% mutate(extra_layer = "PFP"))
# plot sensitivity vs 1 - specificity
est_sens_spec_plot_tib <- est_sens_spec %>%
  filter(missing_perc %in% c(0, 0.4)) %>%
  mutate(fpr = 1 - mn_spec, tpr = mn_sens,
         n_miss = paste0(n, "_", missing_perc),
         missing = factor(missing_perc, labels = c("0%", "40%"), ordered = TRUE),
         est_fct = factor(case_when(
           estimator_type == "lasso" & extra_layer == "SS" ~ 1,
           estimator_type == "SL" & extra_layer == "SS" ~ 2,
           estimator_type == "SL" & extra_layer == "none" ~ 3,
           estimator_type == "SPVIM" & extra_layer == "gFWER" ~ 4,
           estimator_type == "SPVIM" & extra_layer == "PFP" ~ 5
         ), labels = c("lasso + SS", "SL + SS", "SL", "SPVIM + gFWER",
                       "SPVIM + PFP")))

point_size <- 2
point_size <- 1.5
axis_text_size <- 8.5
lgnd_text_size <- 8.5
title_text_size <- 10
label_size <- title_text_size
fig_width <- 20
fig_height <- 12

est_sens_spec_plot_0 <- est_sens_spec_plot_tib %>%
  filter(missing == "0%") %>%
  ggplot(aes(x = fpr, y = tpr, linetype = est_fct, color = est_fct)) +
  geom_line() +
  xlim(c(0, 1)) +
  scale_color_viridis_d(begin = 0, end = 0.5) +
  scale_linetype_discrete() +
  labs(y = "Sensitivity", x = "1 - Specificity", linetype = "Est. procedure", color = "Est. procedure") +
  facet_grid(rows = vars(p), cols = vars(n), scales = "free", space = "fixed",
             labeller = label_both) +
  theme(strip.background = element_blank(), legend.direction = "horizontal",
        legend.position = "bottom", text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size),
        legend.title = element_text(size = lgnd_text_size))

est_sens_spec_plot_40 <- est_sens_spec_plot_tib %>%
  filter(missing == "40%") %>%
  ggplot(aes(x = fpr, y = tpr, linetype = est_fct, color = est_fct)) +
  geom_line() +
  xlim(c(0, 1)) +
  scale_color_viridis_d(begin = 0, end = 0.5) +
  scale_linetype_discrete() +
  labs(y = "Sensitivity", x = "1 - Specificity", linetype = "Est. procedure", color = "Est. procedure") +
  facet_grid(rows = vars(p), cols = vars(n), scales = "free", space = "fixed",
             labeller = label_both) +
  theme(strip.background = element_blank(), legend.direction = "horizontal",
        legend.position = "bottom", text = element_text(size = axis_text_size),
        axis.text = element_text(size = axis_text_size),
        legend.text = element_text(size = lgnd_text_size),
        legend.title = element_text(size = lgnd_text_size))

est_sens_spec_plot <- plot_grid(est_sens_spec_plot_0 + theme(legend.position = "none"),
                                est_sens_spec_plot_40 + theme(legend.position = "none"),
                                labels = "AUTO", nrow = 2)
lgnd <- get_legend(est_sens_spec_plot_40 +
                     theme(legend.spacing.x = unit(0.1, "in")) +
                     guides(linetype = guide_legend(nrow = 2)))
ggsave(
  filename = paste0(plots_dir, args$sim_name, "_roc-curves.png"),
  plot = plot_grid(est_sens_spec_plot, lgnd, ncol = 1, nrow = 2, rel_heights = c(1, .1)),
  device = "png",
  width = fig_width, height = fig_height, units = "cm", dpi = 300,
)

# compute max sensitivity for fixed specificity = 0.85, n = 200, missing = 0
# print out to a table
setting_txt <- "A"
if (grepl("nonlinear", args$sim_name)) {
    if (grepl("nonnormal", args$sim_name)) {
        setting_txt <- "D"
    } else {
        setting_txt <- "C"
    }
} else if (grepl("nonnormal", args$sim_name)) {
    setting_txt <- "B"
}
max_sens <- est_sens_spec %>%
  filter(n == 200, missing_perc == 0, !grepl("PFP", extra_layer)) %>%
  group_by(family, linear, x_dist, missing_type, estimator_type, extra_layer, p) %>%
    filter(abs(mn_spec - 0.85) == min(abs(mn_spec - 0.85))) %>% 
  ungroup()
knitr::kable(max_sens %>% 
               select(-family, -linear, -x_dist, -missing_type, -n, -missing_perc) %>% 
               arrange(p), 
             format = "latex", 
             col.names = c("Estimator", "Extra layer", "p", "Tuning parameter value", "Mean sensitivity", "Mean specificity"), 
             caption = paste0("Mean sensitivity for mean specificity closest to 0.85 ",
             "at $n = 200$ with no missing data, for each $p$ and estimator in setting ", 
             setting_txt, "."), label = paste0("max-spec-", tolower(setting_txt)), 
             row.names = FALSE, digits = 3) %>%
    kableExtra::kable_styling(latex_options = c("scale_down")) %>%
    cat(., file = paste0(plots_dir, args$sim_name, "_max-spec.tex"))
