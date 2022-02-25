# create data analysis plots

# set up -----------------------------------------------------------------------
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())
library("here")
library("data.table")

# read in the results ----------------------------------------------------------
outcomes <- c("mucinous", "high_malignancy")
est_types <- c("lasso_glm", "lasso-KF_glm", "lasso-SS_glm",
               "baseSL_SL", "SL_SL", "SL-SS_SL", "SPVIM_SL", "SPVIM-RR_SL")
algos <- c("1", "2")
file_table <- expand.grid(result = c("perf", "vim"), est = est_types, algo = algos) %>%
    filter(!(est == "SPVIM-RR_SL" & algo == "2"))
perf <- as_tibble(rbindlist(lapply(as.list(apply(file_table %>%
                                                     filter(result == "perf"),
                                                 1, paste, collapse = "_")), function(l) {
    readRDS(here::here("results", "data_analysis", paste0(l, ".rds")))
})))

# create plots of prediction performance ---------------------------------------
unique_procedures <- c("lasso", "SL", "lasso + SS", "SL + SS", "lasso + KF",
                       "SL (benchmark)", "SPVIM + gFWER", "SPVIM + PFP", "SPVIM + FDR",
                       paste0("SPVIM-RR + ", c("gFWER", "PFP", "FDR")))
unique_levels <- c(6, 1, 3, 5, 2, 4, 7:9, 10:12)
perf_summary <- perf %>%
    group_by(outcome, covariates, biomarkers, selection_procedure, pred_alg, algo, n) %>%
    summarize(perf = mean(est, na.rm = TRUE), var_init = mean(var, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(se = sqrt(var_init / n), cil = perf - qnorm(.975) * se,
           ciu = pmin(perf + qnorm(.975) * se, 1),
           est_proc = factor(case_when(
               selection_procedure == "lasso" ~ 2,
               selection_procedure == "SL" ~ 5,
               selection_procedure == "lasso-SS" ~ 3,
               selection_procedure == "SL-SS" ~ 6,
               selection_procedure == "lasso-KF" ~ 4,
               selection_procedure == "baseSL" ~ 1,
               selection_procedure == "SPVIM + gFWER" ~ 7,
               selection_procedure == "SPVIM + PFP" ~ 8,
               selection_procedure == "SPVIM + FDR" ~ 9,
               selection_procedure == "SPVIM-RR + gFWER" ~ 10,
               selection_procedure == "SPVIM-RR + PFP" ~ 11,
               selection_procedure == "SPVIM-RR + FDR" ~ 12
           ), levels = 1:length(unique_levels), labels = unique_procedures[unique_levels]),
           variables = factor(case_when(
               biomarkers == 0 ~ "Clinical variables only",
               covariates == 0 ~ "Biomarkers only",
               biomarkers == 1 & covariates == 1 ~ "Biomarkers and clinical variables"
           )),
           ci_txt = paste0(round(perf, 3), " [", round(cil, 3), ", ", round(ciu, 3), "]"))

for (vars in unique(as.character(perf_summary$variables))) {
    for (i in seq_len(2)) {
        for (j in seq_len(2)) {
            this_outcome <- outcomes[i]
            this_alg <- algos[j]
            this_perf <- perf_summary %>%
                filter(outcome == this_outcome, algo == this_alg,
                       variables == vars)
            this_plot <- this_perf %>%
                ggplot(aes(x = perf, y = est_proc, label = ci_txt)) +
                geom_point() +
                geom_text(aes(x = 0), color = "blue", size = 5, hjust = 0) +
                geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
                geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
                xlim(c(0, 1)) +
                xlab("CV-AUC") +
                ylab("Variable selection procedure")
            ggsave(filename = paste0(here::here("plots"), "/data-analysis_",
                                     this_outcome, "_alg_", this_alg,
                                     "_covariates_", this_perf$covariates[1],
                                     "_biomarkers_", this_perf$biomarkers[1], ".png"),
                   plot = this_plot,
                   width = 10, height = 5, units = "in", dpi = 300)
        }
    }
}
# get selected sets ------------------------------------------------------------
biomarker_labs <- c("CEA", "CEA mucinous call",
                    "ACTB", "Molecules (M) score",
                    "M neoplasia call", "Telomerase (T) score",
                    "T neoplasia call", "AREG (A) score",
                    "A mucinous call", "Glucose (G) score",
                    "G mucinous call", "A and G mucinous call",
                    "Fluorescence (F) score", "F mucinous call",
                    "DNA mucinous call", "DNA neoplasia call (v1)",
                    "DNA neoplasia call (v2)", "MUC3AC score",
                    "MUC5AC score", "Ab score",
                    "Ab neoplasia call")
selection_types <- c(paste0("lasso", c("", "-KF", "-SS")),
                     paste0("SL", c("", "-SS")),
                     paste0("SPVIM", c("", "-RR")))
selected <- as_tibble(
    rbindlist(lapply(as.list(selection_types), function(l) {
        mucinous <- readRDS(here::here("results", "data_analysis",
                           paste0("output_", l, "/impute_10/", "selected-vars_mucinous_0_1_final.rds")))
        high_malignancy <- readRDS(here::here("results", "data_analysis",
                                              paste0("output_", l, "/impute_10/", "selected-vars_high_malignancy_0_1_final.rds")))
        bind_rows(mucinous %>% mutate(outcome = "mucinous"), high_malignancy %>% mutate(outcome = "high_malignancy"))
    }))
) %>%
    filter(selected != "institution") %>%
    mutate(procedure = case_when(
        procedure == "lasso" ~ "lasso",
        procedure == "SL" ~ "SL",
        procedure == "lasso-SS" ~ "lasso + SS",
        procedure == "SL-SS" ~ "SL + SS",
        procedure == "lasso-KF" ~ "lasso + KF",
        procedure == "SPVIM + gFWER" ~ "SPVIM + gFWER",
        procedure == "SPVIM + PFP" ~ "SPVIM + PFP",
        procedure == "SPVIM + FDR" ~ "SPVIM + FDR",
        procedure == "SPVIM-RR + gFWER" ~ "SPVIM-RR + gFWER",
        procedure == "SPVIM-RR + PFP" ~ "SPVIM-RR + PFP",
        procedure == "SPVIM-RR + FDR" ~ "SPVIM-RR + FDR"
    ),
    yesno = "Yes") %>%
    pivot_wider(names_from = procedure, values_from = yesno) %>%
    mutate(across(where(is.character), ~ ifelse(is.na(.x), "No", .x)),
           `Number of procedures` = rowSums(across(-selected, ~ .x == "Yes")),
           selected = factor(case_when(
               selected == "cea" ~ 1,
               selected == "cea_call" ~ 2,
               selected == "lab1_actb" ~ 3,
               selected == "lab1_molecules_score" ~ 4,
               selected == "lab1_molecules_neoplasia_call" ~ 5,
               selected == "lab1_telomerase_score" ~ 6,
               selected == "lab1_telomerase_neoplasia_call" ~ 7,
               selected == "lab4_areg_score" ~ 8,
               selected == "lab4_areg_mucinous_call" ~ 9,
               selected == "lab4_glucose_score" ~ 10,
               selected == "lab4_glucose_mucinous_call" ~ 11,
               selected == "lab4_combined_mucinous_call" ~ 12,
               selected == "lab2_fluorescence_score" ~ 13,
               selected == "lab2_fluorescence_mucinous_call" ~ 14,
               selected == "lab5_mucinous_call" ~ 15,
               selected == "lab5_pancrea_seq_v1_neoplasia_call" ~ 16,
               selected == "lab5_pancrea_seq_v2_neoplasia_call" ~ 17,
               selected == "lab3_muc3ac_score" ~ 18,
               selected == "lab3_muc5ac_score" ~ 19,
               selected == "lab6_m_ab_das_score" ~ 20,
               selected == "lab6_m_ab_das_neoplasia_call" ~ 21
           ), levels = 1:21, labels = biomarker_labs)) %>%
    rename(Biomarker = selected)

for (i in seq_len(2)) {
    this_outcome <- outcomes[i]
    outcome_txt <- ifelse(this_outcome == "mucinous", "is mucinous", "has high malignancy potential")
    def_txt <- "Table~\\ref{tab:biomarker_definitions}"
    these_selected <- selected %>%
        filter(outcome == this_outcome, !is.na(selected)) %>%
        select(-outcome) %>%
        arrange(Biomarker)
    knitr::kable(these_selected, format = "latex",
                 caption = paste0("Biomarkers selected by each selection procedure for ",
                                  "predicting whether a cyst ", outcome_txt,
                                  " on the full imputed dataset. Full definitions of ",
                                  "each variable are provided in ", def_txt, "."),
                 label = this_outcome, escape = FALSE) %>%
        kableExtra::kable_styling() %>%
        kableExtra::column_spec(column = 1, width = "10em") %>%
        kableExtra::column_spec(column = 3, width = "3em") %>%
        kableExtra::column_spec(column = 4, width = "3em") %>%
        kableExtra::column_spec(column = 6, width = "3em") %>%
        kableExtra::column_spec(column = 7, width = "3em") %>%
        kableExtra::column_spec(column = 8, width = "3em") %>%
        kableExtra::column_spec(column = 9, width = "5em") %>%
        cat(file = paste0(here::here("plots"), "/data-analysis_", this_outcome, "_selected.tex"))
}
