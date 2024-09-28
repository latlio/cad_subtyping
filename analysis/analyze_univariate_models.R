# Ideally run this on Minerva
# Author: Lathan Liou
# Contact: lathan.liou@icahn.mssm.edu
# Inputs: output of analyze_prset.R

library(optparse)

option_list <- list(
  make_option(c("--original_ukb_path"), type = "character", dest = "original_ukb_path"),
  make_option(c("--pheno_df_path"), type = "character", dest = "pheno_df_path"),
  make_option(c("--results_path"), type = "character", dest = "results_path"),
  make_option(c("--train_prop"), type = "numeric", dest = "train_prop"),
  make_option(c("--run_name"), type = "character", dest = "run_name"),
  make_option(c("--prsice_run_name"), type = "character", dest = "prsice_run_name"),
  make_option(c("--prset_run_name"), type = "character", dest = "prset_run_name"),
  make_option(c("--outcome_var"), type = "character", dest = "outcome_var"),
  make_option(c("--seed"), type = "numeric", dest = "seed"),
  make_option(c("--threshold"), type = "numeric", dest = "threshold")
)

argv <- commandArgs(trailingOnly = TRUE)
argv <- parse_args(OptionParser(option_list = option_list))

library(tidyverse)
library(data.table)
library(here)
library(janitor)
library(glmnet)
library(conflicted)
library(glue)
library(pROC)
library(caret)
library(PredictABEL)
library(ggsci)
library(ggrepel)
library(boot)
library(patchwork)
library(qs)
library(lmtest)
library(doParallel)
# nodes <- detectCores() - 2
# cl <- makeCluster(nodes)
# registerDoParallel(cl)
# library(glmnetUtils)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
# conflicted::conflict_prefer("cv.glmnet", "glmnetUtils")

# source("cad_subtype/R/common_prset_fcts.R")
source("common_prset_fcts.R")

build_and_train_bin <- function(formula,
                                .train_data,
                                .test_data,
                                .outcome_var,
                                .threshold,
                                .results_path,
                                .run_name,
                                file_name = "") {
  model <- glm(as.formula(formula),
               data = .train_data,
               family = binomial)
  train_pred <- predict(model, new = .train_data %>% select(-all_of(.outcome_var)),
                        type = "response")
  test_pred <- predict(model, new = .test_data %>% select(-all_of(.outcome_var)),
                       type = "response")
  # test_pred_class <- ifelse(test_pred > .threshold, 1, 0)
  
  #get AUC
  train_auc <- pROC::auc(.train_data %>% pull(!!enquo(.outcome_var)), train_pred)
  test_auc <- pROC::auc(.test_data %>% pull(!!enquo(.outcome_var)), test_pred)
  
  #get roc
  train_roc <- roc(.train_data %>% 
                     pull(!!enquo(.outcome_var)), 
                   train_pred)
  
  test_roc <- roc(.test_data %>% 
                    pull(!!enquo(.outcome_var)), 
                  test_pred)
  
  # get Nagelkerke's R2 of train
  # r2 <- compute_nagelkerkeR2(model)
  
  test_discrim_df <- evaluate_single_classifier(tibble(outcome = .test_data %>% 
                                                         pull(!!enquo(.outcome_var)),
                                                       pred_prob = test_pred) %>%
                                                  mutate(predicted_label = ifelse(pred_prob > .threshold, 1, 0)))
  
  message("Bootstrapping...\n")
  
  prs_boot_ci <- ci.auc(.test_data %>% 
                          pull(!!enquo(.outcome_var)), 
                        test_pred,
                        boot.n = 5000,
                        parallel = TRUE) 
  qsave(prs_boot_ci, paste0(.results_path, "/", .run_name, file_name,
                            "_prs_only_boot.qs"))
  
  out <- list(model = model,
              train_auc = train_auc,
              test_auc = test_auc,
              train_roc = train_roc,
              test_roc = test_roc,
              test_discrim_df = test_discrim_df,
              test_pred = test_pred)
  # train_r2 = r2)
  
  return(out)
}

boot_auc <- function(data, indices, formula, .test_data, .outcome_var) {
  d <- data[indices,] #allows boot to select sample
  model <- glm(formula, data=d)
  
  test_pred <- predict(model, new = d %>% select(-all_of(.outcome_var)))
  
  test_auc <- pROC::auc(d %>% pull(!!enquo(.outcome_var)), test_pred)
  
  return(test_auc)
}

boot_lasso_auc <- function(data, indices, fit, .outcome_var, .lambda, 
                           .pred_vars) {
  
  d <- data[indices,] #allows boot to select sample
  
  test_pred <- predict(fit, 
                       newx = d %>% 
                         select(all_of(.pred_vars)) %>% 
                         as.matrix(),
                       type = "response",
                       s = .lambda)
  
  # test_auc <- pROC::auc(d %>% pull(!!enquo(.outcome_var)), test_pred)
  test_auc <- pROC::auc(d$ldl_bin, test_pred)
  
  return(test_auc)
}

test_auc_difference <- function(.results_path,
                                .run_name) {
  genome_test <- map(seq(length(univariate_models)),
                     ~ pROC::roc.test(univariate_models[[.x]]$test_roc, 
                                      genome_models[[.x]]$test_roc,
                                      boot.n = 5000,
                                      parallel = TRUE))
  
  prset_test <- map(seq(length(univariate_models)),
                    ~ pROC::roc.test(univariate_models[[.x]]$test_roc, 
                                     prset_models[[.x]]$test_roc,
                                     boot.n = 5000,
                                     parallel = TRUE))
  
  genome_delta <- map(seq(length(genome_test)), 
                      ~ genome_test[[.x]]$estimate[2] - genome_test[[.x]]$estimate[1]) %>%
    unlist()
  
  genome_delta_p <- map(seq(length(genome_test)), 
                        ~ genome_test[[.x]]$p.value) %>%
    unlist()
  
  prset_delta <- map(seq(length(prset_test)), 
                     ~ prset_test[[.x]]$estimate[2] - prset_test[[.x]]$estimate[1]) %>%
    unlist()
  
  prset_delta_p <- map(seq(length(prset_test)), 
                       ~ prset_test[[.x]]$p.value) %>%
    unlist()
  
  out <- tibble(
    variable = formulas$univariate,
    genome_delta = genome_delta,
    genome_delta_p = genome_delta_p,
    prset_delta = prset_delta,
    prset_delta_p = prset_delta_p
  )
  
  write_csv(out, paste0(.results_path, "/", .run_name,
                        "_auc_differences.csv"))
  return(out)
}

build_and_train_cont <- function(formula,
                                 .train_data,
                                 .test_data,
                                 .outcome_var) {
  model <- lm(formula, .train_data)
  train_pred <- predict(model, new = .train_data %>% select(-all_of(.outcome_var)))
  in_sample_cor <- cor.test(train_pred, .train_data[[.outcome_var]])
  
  test_pred <- predict(model, new = .test_data %>% select(-pheno))
  out_sample_cor <- cor.test(test_pred, .test_data[[.outcome_var]])
  
  out <- list(model = model,
              performance = tibble(
                train_cor = in_sample_cor$estimate^2,
                train_cor_p = in_sample_cor$p.value,
                test_cor = out_sample_cor$estimate^2,
                test_cor_p = out_sample_cor$p.value
              ),
              test_pred = test_pred)
  return(out)
}

build_and_lasso_bin <- function(pred_vars,
                                .train_data,
                                .test_data,
                                .outcome_var,
                                .threshold,
                                .results_path,
                                .run_name,
                                file_name = "",
                                .alpha = 1,
                                .lambda_type = "lambda.1se") {
  
  message("Training lasso model...\n")
  
  train_mat <- as.matrix(.train_data %>% select(all_of(pred_vars)))
  model <- cv.glmnet(y = as.matrix(.train_data[[.outcome_var]]),
                     x = train_mat,
                     alpha = .alpha,
                     family = binomial,
                     standardize = TRUE)
  # penalty.factor = c(rep(0, 18), rep(1, ncol(train_mat) - 18)))
  train_pred <- predict(model, newx = train_mat, type = "response", s = .lambda_type)
  test_pred <- predict(model, newx = as.matrix(.test_data %>% 
                                                 select(all_of(pred_vars))),
                       type = "response", s = .lambda_type)
  
  message("Evaluating model performance...\n")
  
  #get AUC
  train_auc <- pROC::auc(.train_data %>% pull(!!enquo(.outcome_var)), train_pred)
  test_auc <- pROC::auc(.test_data %>% pull(!!enquo(.outcome_var)), test_pred)
  print(test_auc)
  
  #get roc
  train_roc <- roc(.train_data %>% 
                     pull(!!enquo(.outcome_var)), 
                   train_pred)
  test_roc <- roc(.test_data %>% 
                    pull(!!enquo(.outcome_var)), 
                  test_pred)
  
  # get Nagelkerke's R2 (?) 
  # r2 <- 1 - model$cvm/var(as.matrix(.train_data[[.outcome_var]]))
  
  test_discrim_df <- evaluate_single_classifier(tibble(outcome = .test_data %>% 
                                                         pull(!!enquo(.outcome_var)),
                                                       pred_prob = test_pred) %>%
                                                  mutate(predicted_label = ifelse(pred_prob > .threshold, 1, 0)))
  
  message("Bootstrapping...\n")
  
  # bootlasso <- HDCI::bootLasso(x = as.matrix(.test_data %>% 
  #                                 select(all_of(pred_vars))),
  #                 y = .test_data %>% 
  #                   pull(!!enquo(.outcome_var)),
  #                 B = 5000,
  #                 nfolds = 5,
  #                 parallel = TRUE,
  #                 standardize = FALSE,
  #                 parallel.boot = TRUE,
  #                 ncores.boot = parallel::detectCores() - 2)
  # prset_boot <- boot(data = .test_data,
  #                    statistic = boot_lasso_auc,
  #                    stype = "i",
  #                    R = 1000,
  #                    parallel = "multicore",
  #                    ncpus = parallel::detectCores() - 2,
  #                    fit = model,
  #                    .outcome_var = .outcome_var,
  #                    .lambda = .lambda_type,
  #                    .pred_vars = pred_vars)
  # prset_boot_ci <- boot.ci(prset_boot,
  #                          type = "norm")
  prset_boot_ci <- ci.auc(.test_data %>% 
                            pull(!!enquo(.outcome_var)), 
                          test_pred,
                          boot.n = 5000,
                          parallel = TRUE) 
  qsave(prset_boot_ci, paste0(.results_path, "/", .run_name, file_name,
                              "_prset_only_boot_lasso.qs"))
  
  out <- list(model = model,
              train_auc = train_auc,
              test_auc = test_auc,
              train_roc = train_roc,
              test_roc = test_roc,
              test_discrim_df = test_discrim_df,
              test_pred = test_pred)
  # train_r2 = r2)
  return(out)
}

build_and_lasso_cont <- function(pred_vars,
                                 .train_data,
                                 .test_data,
                                 .outcome_var,
                                 .alpha = 1,
                                 .lambda_type = "lambda.1se") {
  train_mat <- as.matrix(.train_data %>% select(all_of(pred_vars)))
  model <- cv.glmnet(y = as.matrix(.train_data[[.outcome_var]]),
                     x = train_mat,
                     alpha = .alpha,
                     standardize = FALSE)
  # penalty.factor = c(rep(0, 18), rep(1, ncol(train_mat) - 18)))
  train_pred <- predict(model, newx = train_mat, type = "response", s = .lambda_type)
  test_pred <- predict(model, newx = as.matrix(.test_data %>% 
                                                 select(all_of(pred_vars))),
                       type = "response", s = .lambda_type)
  
  in_sample_cor <- cor.test(train_pred, .train_data[[.outcome_var]])
  out_sample_cor <- cor.test(test_pred, .test_data[[.outcome_var]])
  
  out <- list(model = model,
              performance = tibble(
                train_cor = in_sample_cor$estimate^2,
                train_cor_p = in_sample_cor$p.value,
                test_cor = out_sample_cor$estimate^2,
                test_cor_p = out_sample_cor$p.value
              ),
              test_pred = test_pred)
  return(out)
}

run_analysis <- function(original_ukb_path = "clean_original_ukb.csv",
                         pheno_df_path = "ldl_bin_mod1.csv",
                         results_path = "cad_subtype/result_data",
                         train_prop = 0.8,
                         run_name = "ukb_cad_ldlsub_bin",
                         prsice_run_name = "prsice_ukb_cad_ldlsub_bin",
                         prset_run_name = "prset_ukb_cad_ldlsub_bin",
                         outcome_var = "pheno",
                         seed = 47,
                         threshold = 0.5,
                         env_sensitivity = FALSE) {
  browser()
  # Loading and Preparing Data ----
  message("Preparing data ... \n")
  pheno_df <- fread(pheno_df_path)
  original_ukb_df <- fread(original_ukb_path)
  prsice_best <- fread(here(results_path, "/", paste0(prsice_run_name, ".best"))) %>%
    janitor::clean_names() %>%
    select(-in_regression)
  
  # standardize <- function(x, n) {
  #   (x - mean(x)) / (sqrt(sum(x^2) / n))
  # }
  # load in lasso-selected pathways
  
  prset_best <- fread(here(results_path, "/", paste0(prset_run_name, ".best"))) %>%
    janitor::clean_names() %>%
    select(-c(in_regression, base)) %>% 
    mutate(across(!fid & !iid, \(x) standardize(x, n = n())))
  # mutate(across(!fid & !iid, ~scale(.x)[,1]))
  
  prset_summary <- fread(here(results_path, "/", paste0(prset_run_name, ".summary"))) %>%
    janitor::clean_names() %>%
    filter(competitive_p < 0.05) %>%
    mutate(set = tolower(set),
           set = str_replace(set, ":", "_"))
  
  select_prset_best <- prset_best %>%
    select(fid, iid, all_of(prset_summary$set))
  
  # prset_lasso <- qread(here(results_path, "/", paste0(prset_run_name, "_results.qs")))
  # 
  # lasso_selected_pathways <- save_model_coef(prset_lasso$model) %>%
  #   filter(go_id != "(Intercept)") %>%
  #   mutate(go_id = stringr::str_to_lower(go_id),
  #          go_id = str_replace_all(go_id, ":", "_")) %>%
  #   pull(go_id)
  # 
  # select_prset_best <- prset_best %>%
  #   select(fid, iid, all_of(lasso_selected_pathways))
  
  # 8/1/2024: Adding Lpa as univariate per ATVB reviewer comment
  
  univariate <- c("ascvd", "apoa_std", "apob_std", "trigly_std", "ldl_std", "hdl_std",
                  "crp_std", "lpa_std")
  outcomes <- c("cad_bin", "unstable_cad_bin", "occlusive_cad_bin", "ldl_bin",
                "lpa_bin", "stemi_cad_bin")
  prset_colnames <- colnames(select_prset_best %>% select(-c(fid, iid)))
  
  combined_df <- pheno_df %>%
    left_join(original_ukb_df %>%
                select(fid, iid, chol, hdl_chol, sys_bp_auto, smoke, diabetes,
                       apoa_std, apob_std, trigly_std, ldl_std, hdl_std,
                       crp_std, lpa_std, bp_med_bin),
              by = c("fid", "iid")) %>% 
    drop_na(age, chol, hdl_chol, sys_bp_auto, smoke, diabetes) %>% 
    mutate(ascvd = case_when(
      # untreated male
      sex == 0 & bp_med_bin == 0 ~ 1 - 0.9144^exp((12.344*log(age) + 11.853*log(chol) - 2.664*log(age)*log(chol) - 7.99*log(hdl_chol) +
                                                     1.769*log(age)*log(hdl_chol) + 1.764*log(sys_bp_auto) + 7.837*smoke - 1.795*log(age)*smoke + 0.658*diabetes) - 61.18),
      # treated male
      sex == 0 & bp_med_bin == 1 ~ 1 - 0.9144^exp((12.344*log(age) + 11.853*log(chol) - 2.664*log(age)*log(chol) - 7.99*log(hdl_chol) +
                                                     1.769*log(age)*log(hdl_chol) + 1.797*log(sys_bp_auto) + 7.837*smoke - 1.795*log(age)*smoke + 0.658*diabetes) - 61.18),
      # untreated female
      sex == 1 & bp_med_bin == 0 ~ 1 - 0.9665^exp((-29.799*log(age) + 4.884*log(age)^2 + 13.540*log(chol) - 3.114*log(age)*log(chol) - 13.578*log(hdl_chol) +
                                                     3.149*log(age)*log(hdl_chol) + 1.957*log(sys_bp_auto) + 7.574*smoke - 1.665*log(age)*smoke + 0.661*diabetes) + 29.18),
      # treated female
      sex == 1 & bp_med_bin == 1 ~ 1 - 0.9665^exp((-29.799*log(age) + 4.884*log(age)^2 + 13.540*log(chol) - 3.114*log(age)*log(chol) - 13.578*log(hdl_chol) +
                                                     3.149*log(age)*log(hdl_chol) + 2.019*log(sys_bp_auto) + 7.574*smoke - 1.665*log(age)*smoke + 0.661*diabetes) + 29.18)
    ),
    ascvd_cad_bin = ifelse(ascvd > 0.075, 1, 0)) %>%
    left_join(prsice_best) %>%
    left_join(select_prset_best) %>%
    drop_na(c(all_of(univariate), prs, colnames(select_prset_best))) %>%
    mutate(prs_std = scale(prs)[,1])
  # ses, diet, pa, pm25))
  
  # Supplementary Exploratory Plots ----
  ## Density of PRS by sex
  message("Making exploratory plots ... \n")
  prs_density_sex_plot <- ggplot(combined_df %>%
                                   mutate(sex = as.factor(case_when(
                                     sex == 0 ~ glue("Male (N = {sum(combined_df$sex == 0)})"),
                                     sex == 1 ~ glue("Female (N = {sum(combined_df$sex == 1)})"),
                                   ))), aes(x = prs_std, fill = sex)) + 
    geom_density(alpha = 0.4) + 
    theme_bw() +
    labs(x = "PRS (standardized)",
         y = "Density",
         fill = "Sex") + 
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.8,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_density_sex_plot.png"), 
         prs_density_sex_plot,
         width = 10,
         height = 9)
  
  ## Scatter plot of PRS by age colored by sex
  prs_scatter_age_sex_plot <- ggplot(combined_df %>%
                                       mutate(sex = as.factor(case_when(
                                         sex == 0 ~ glue("Male (N = {sum(combined_df$sex == 0)})"),
                                         sex == 1 ~ glue("Female (N = {sum(combined_df$sex == 1)})"),
                                       ))), aes(x = age, 
                                                y = prs_std,
                                                color = sex)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "Age",
         y = "PRS (standardized)",
         color = "Sex") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.2,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_age_sex_plot.png"), 
         prs_scatter_age_sex_plot,
         width = 10,
         height = 9)
  
  ## Scatter plot of PRS by univariate
  prs_scatter_ascvd_plot <- ggplot(combined_df, aes(x = ascvd, 
                                                    y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "ASCVD",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_ascvd_plot.png"), 
         prs_scatter_ascvd_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$ascvd, combined_df$prs_std)
  
  prs_scatter_apoa_plot <- ggplot(combined_df, aes(x = apoa_std, 
                                                   y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "ApoA (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_apoa_plot.png"), 
         prs_scatter_apoa_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$apoa_std, combined_df$prs_std)
  
  prs_scatter_apob_plot <- ggplot(combined_df, aes(x = apob_std, 
                                                   y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "ApoB (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_apob_plot.png"), 
         prs_scatter_apob_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$apob_std, combined_df$prs_std)
  
  prs_scatter_ldl_plot <- ggplot(combined_df, aes(x = ldl_std, 
                                                  y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "LDL (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_ldl_plot.png"), 
         prs_scatter_ldl_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$ldl_std, combined_df$prs_std)
  
  prs_scatter_hdl_plot <- ggplot(combined_df, aes(x = hdl_std, 
                                                  y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "HDL (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_hdl_plot.png"), 
         prs_scatter_hdl_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$hdl_std, combined_df$prs_std)
  
  prs_scatter_tri_plot <- ggplot(combined_df, aes(x = trigly_std, 
                                                  y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "Triglycerides (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_tri_plot.png"), 
         prs_scatter_tri_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$trigly_std, combined_df$prs_std)
  
  prs_scatter_crp_plot <- ggplot(combined_df, aes(x = crp_std, 
                                                  y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "C-Reactive Protein (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_crp_plot.png"), 
         prs_scatter_crp_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$crp_std, combined_df$prs_std)
  
  prs_scatter_lpa_plot <- ggplot(combined_df, aes(x = lpa_std, 
                                                  y = prs_std)) + 
    geom_point() + 
    theme_bw() +
    labs(x = "Lp(a) (standardized)",
         y = "PRS (standardized)") +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position=c(0.1,0.9))
  ggsave(paste0(results_path, "/", run_name, "_prs_scatter_lpa_plot.png"), 
         prs_scatter_lpa_plot,
         width = 10,
         height = 9)
  cor.test(combined_df$lpa_std, combined_df$prs_std)
  
  # Train models ----
  set.seed(seed)
  train_df <- combined_df %>%
    sample_frac(train_prop)
  test_df <- combined_df %>% 
    dplyr::setdiff(train_df)
  
  formulas <- expand.grid(outcomes = outcomes, univariate = univariate) %>%
    mutate(univariate_formula = glue('{outcomes} ~ {univariate} + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'),
           genome_formula = glue("{outcomes} ~ {univariate} + prs_std + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15"),
           pathway_formula = glue("{outcomes} ~ {univariate} + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + {paste(prset_colnames, collapse = '+')}"))
  
  full_adj_formulas <- tibble(outcomes = outcomes) %>%
    mutate(full_adj_genome_formula = glue("{outcomes} ~ prs_std + {paste(univariate, collapse = '+')}"),
           full_adj_pathway_formula = glue("{outcomes} ~ {paste(prset_colnames, collapse = '+')} + {paste(univariate, collapse = '+')}"))
  
  env_formulas <- expand.grid(outcomes = outcomes, univariate = univariate) %>%
    mutate(form1 = glue('{outcomes} ~ {univariate} + age + age2 + sex + centre + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'),
           form2 = glue('{outcomes} ~ {univariate} + ses + age + age2 + sex + centre + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'),
           form3 = glue('{outcomes} ~ {univariate} + diet + pa + age + age2 + sex + centre + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'),
           form4 = glue('{outcomes} ~ {univariate} + ses + diet + pa + age + age2 + sex + centre + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'),
           form5 = glue('{outcomes} ~ {univariate} + pm25 + ses + diet + pa + age + age2 + sex + centre + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15'))
  
  switch(outcome_var,
         pheno = {
           formulas <- tibble(outcomes = "pheno",
                              univariate = univariate) %>%
             mutate(univariate_formula = glue("{outcomes} ~ {univariate} + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15"),
                    genome_formula = glue("{outcomes} ~ {univariate} + prs_std + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15"),
                    pathway_formula = glue("{outcomes} ~ {univariate} + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15 + {paste(prset_colnames, collapse = '+')}"))
           
           full_adj_formulas <- tibble(outcomes = "pheno",
                                       univariate = univariate) %>%
             mutate(full_adj_genome_formula = glue("{outcomes} ~ prs_std + {paste(univariate, collapse = '+')}"),
                    full_adj_pathway_formula = glue("{outcomes} ~ {paste(prset_colnames, collapse = '+')} + {paste(univariate, collapse = '+')}"))
           
           if(str_detect(run_name, "ldlsub")) {
             formulas <- formulas %>% 
               filter(!univariate %in% c("ldl_std", "apob_std"))
           } else if(str_detect(run_name, "ascvd_cad")) {
             formulas <- formulas %>% 
               filter(univariate != "ascvd") %>%
               mutate(univariate_formula = str_replace(univariate_formula, fixed("age + age2 + sex + "), ""),
                      genome_formula = str_replace(genome_formula, fixed("age + age2 + sex + "), ""),
                      pathway_formula = str_replace(pathway_formula, fixed("age + age2 + sex + "), ""))
           } else {
             formulas
           }
         },
         cad_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var)
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var)
         }, 
         ldl_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var,
                    !univariate %in% c("ldl_std", "apob_std"))
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var,
                    !univariate %in% c("ldl_std", "apob_std"))
         },
         lpa_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var)
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var)
         },
         stemi_cad_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var)
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var)
         },
         unstable_cad_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var)
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var)
         },
         occlusive_cad_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var)
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var)
         },
         ascvd_cad_bin = {
           formulas <- formulas %>%
             filter(outcomes == outcome_var,
                    univariate != "ascvd") %>%
             mutate(univariate_formula = str_replace(univariate_formula, fixed("age + age2 + sex + "), ""),
                    genome_formula = str_replace(genome_formula, fixed("age + age2 + sex + "), ""),
                    pathway_formula = str_replace(pathway_formula, fixed("age + age2 + sex + "), ""))
           
           env_formulas <- env_formulas %>%
             filter(outcomes == outcome_var,
                    univariate != "ascvd")
         })
  
  if(str_detect(pheno_df_path, "ldl")) {
    plot_title <- "LDL Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := as.factor(case_when(
        !!sym(outcome_var) == 0 ~ "Normal LDL",
        !!sym(outcome_var) == 1 ~ "High LDL",
      )))
    subtype_pie_colors <- c(
      "High LDL" = "#F8766D",
      "Normal LDL" = "#00BFC4"
    )
  } else if(str_detect(pheno_df_path, "lpa")) {
    plot_title <- "Lpa Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "Normal Lpa",
        !!sym(outcome_var) == 1 ~ "High Lpa",
      ))
    subtype_pie_colors <- c(
      "High Lpa" = "#F8766D",
      "Normal Lpa" = "#00BFC4"
    )
  } else if(str_detect(pheno_df_path, "stemi")) {
    plot_title <- "STEMI CAD Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "NSTEMI",
        !!sym(outcome_var) == 1 ~ "STEMI",
      ))
    subtype_pie_colors <- c(
      "STEMI" = "#F8766D",
      "NSTEMI" = "#00BFC4"
    )
  } else if(str_detect(pheno_df_path, "unstable")) {
    plot_title <- "Unstable CAD Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "Stable",
        !!sym(outcome_var) == 1 ~ "Unstable",
      ))
    subtype_pie_colors <- c(
      "Unstable" = "#F8766D",
      "Stable" = "#00BFC4"
    )
  } else if(str_detect(pheno_df_path, "occlusive")) {
    plot_title <- "Occlusive CAD Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "Non-occlusive",
        !!sym(outcome_var) == 1 ~ "Occlusive",
      ))
    subtype_pie_colors <- c(
      "Occlusive" = "#F8766D",
      "Non-occlusive" = "#00BFC4"
    )
  } else if(str_detect(pheno_df_path, "ascvd")) {
    plot_title <- "ASCVD Risk Score Subtype"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "Low ASCVD Score",
        !!sym(outcome_var) == 1 ~ "High ASCVD Score",
      ))
    subtype_pie_colors <- c(
      "High ASCVD Score" = "#F8766D",
      "Low ASCVD Score" = "#00BFC4"
    )
  } else {
    plot_title <- "CAD Case/Control"
    subtype_df <- combined_df %>%
      group_by(!!sym(outcome_var)) %>%
      summarize(n = n(),
                prop = n()/nrow(.) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5*prop) %>%
      mutate("{outcome_var}" := case_when(
        !!sym(outcome_var) == 0 ~ "No CAD",
        !!sym(outcome_var) == 1 ~ "CAD",
      ))
    subtype_pie_colors <- c(
      "CAD" = "#F8766D",
      "No CAD" = "#00BFC4"
    )
  }
  
  if(str_detect(run_name, "_bin")) {
    
    if(env_sensitivity) {
      #clinical variable only
      models1 <- map(env_formulas$form1, 
                     ~build_and_train_bin(formula = .x,
                                          .train_data = train_df,
                                          .test_data = test_df,
                                          .outcome_var = outcome_var,
                                          .threshold = threshold))
      
      #clinical variable + ses
      models2 <-  map(env_formulas$form2, 
                      ~build_and_train_bin(formula = .x,
                                           .train_data = train_df,
                                           .test_data = test_df,
                                           .outcome_var = outcome_var,
                                           .threshold = threshold))
      
      #clinical variable + lifestyle (diet, pa)
      models3 <- map(env_formulas$form3, 
                     ~build_and_train_bin(formula = .x,
                                          .train_data = train_df,
                                          .test_data = test_df,
                                          .outcome_var = outcome_var,
                                          .threshold = threshold))
      
      #clinical variable + ses + lifestyle
      models4 <-  map(env_formulas$form4, 
                      ~build_and_train_bin(formula = .x,
                                           .train_data = train_df,
                                           .test_data = test_df,
                                           .outcome_var = outcome_var,
                                           .threshold = threshold))
      #pm25
      models5 <- map(env_formulas$form5, 
                     ~build_and_train_bin(formula = .x,
                                          .train_data = train_df,
                                          .test_data = test_df,
                                          .outcome_var = outcome_var,
                                          .threshold = threshold))
      
      mod1_aucs <- tibble(auc = map(models1, function(x) {
        x$test_auc
      }) %>%
        flatten_dbl(),
      variable = env_formulas$univariate,
      model = "1: clinical variable")
      
      mod2_aucs <- tibble(auc = map(models2, function(x) {
        x$test_auc
      }) %>%
        flatten_dbl(),
      variable = env_formulas$univariate,
      model = "2: clinical variable + ses")
      
      mod3_aucs <- tibble(auc = map(models3, function(x) {
        x$test_auc
      }) %>%
        flatten_dbl(),
      variable = env_formulas$univariate,
      model = "3: clinical variable + lifestyle")
      
      mod4_aucs <- tibble(auc = map(models4, function(x) {
        x$test_auc
      }) %>%
        flatten_dbl(),
      variable = env_formulas$univariate,
      model = "4: clinical variable + ses + lifestyle")
      
      mod5_aucs <- tibble(auc = map(models5, function(x) {
        x$test_auc
      }) %>%
        flatten_dbl(),
      variable = env_formulas$univariate,
      model = "5: clinical variable + pm25 + all")
      
      all_env_aucs <- bind_rows(mod1_aucs,
                                mod2_aucs,
                                mod3_aucs,
                                mod4_aucs,
                                mod5_aucs) %>%
        mutate(variable = factor(variable, levels=c("ascvd",
                                                    "apoa_std",
                                                    "apob_std",
                                                    "ldl_std",
                                                    "hdl_std",
                                                    "trigly_std",
                                                    "crp_std")))
      all_env_auc_plot <- ggplot(all_env_aucs, aes(x = variable, y = auc, fill = model)) +
        geom_col(position = position_dodge2(width = .9, preserve = "single")) +
        # geom_text(aes(x = variable, y = auc, 
        #               label = glue("{format(round(auc, 2), nsmall = 2)}")),
        #           position = position_dodge2(width = .9, preserve = "single"),
        #           size = 3.5,
        #           stat = "identity",
        #           vjust = 0) +
        theme_bw() +
        theme(plot.title = element_text(size = 20),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size = 16),
              panel.grid = element_blank(), 
              # legend.position ="bottom",
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
              axis.text.y = element_text(size = 16)) +
        coord_cartesian(ylim = c(0.5, 1)) +
        ggsci::scale_fill_npg() +
        labs(x = "Variable",
             y = "Test AUC",
             fill = "Model",
             # title = gsub("cad_subtype/output_data/(.*?)_bin_mod1\\.csv", 
             #              '\\1', pheno_df_path),
             title = plot_title)
    }
    # Running univariate models ----
    message("Running PRS-only model")
    if(outcome_var == "ascvd_cad_bin") {
      prs_formula <- glue("{outcome_var} ~ prs_std + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15")
    } else {
      prs_formula <- glue("{outcome_var} ~ prs_std + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15")
    }
    prs_only_model <- build_and_train_bin(formula = prs_formula,
                                          .train_data = train_df,
                                          .test_data = test_df,
                                          .outcome_var = outcome_var,
                                          .threshold = threshold,
                                          .results_path = results_path,
                                          .run_name = run_name)
    #prs_only_model <- qread(paste0(results_path, "/", run_name,"_prs_only_bin_results.qs"))
    qsave(prs_only_model, paste0(results_path, "/", run_name,
                                 "_prs_only_bin_results.qs"))
    prs_only_boot_ci <- qread(paste0(results_path, "/", run_name,
                                     "_prs_only_boot.qs"))
    
    # prs_only_boot <- boot(data = test_df,
    #                       statistic = boot_auc,
    #                       stype = "i",
    #                       R = 5000,
    #                       parallel = "multicore",
    #                       ncpus = parallel::detectCores() - 2,
    #                       formula = prs_formula,
    #                       .outcome_var = outcome_var)
    # prs_only_boot_ci <- boot.ci(prs_only_boot,
    #                             type = "norm")
    # qsave(prs_only_boot_ci, paste0(results_path, "/", run_name,
    #                                  "_prs_only_boot.qs"))
    
    message("Running PRSet-only model")
    if(outcome_var == "ascvd_cad_bin") {
      prset_pred_vars <- c(paste0("pc", 1:15), prset_colnames)
    } else {
      prset_pred_vars <- c("age", "age2", "sex",
                           paste0("pc", 1:15), prset_colnames)
    }
    prset_only_model <- build_and_lasso_bin(pred_vars = prset_pred_vars,
                                            .train_data = train_df,
                                            .test_data = test_df,
                                            .outcome_var = outcome_var,
                                            .threshold = threshold,
                                            .results_path = results_path,
                                            .run_name = run_name)
    #prset_only_model <- qread(paste0(results_path, "/", run_name, "_prset_only_bin_results.qs"))
    qsave(prset_only_model, paste0(results_path, "/", run_name,
                                   "_prset_only_bin_results.qs"))
    prset_only_boot_ci <- qread(paste0(results_path, "/", run_name, "_prset_only_boot_lasso.qs"))
    
    # Lasso CV plot
    png(paste0(results_path, "/", run_name, "_prset_only_model_cv_plot.png"),
        width = "720",
        height = "648")
    plot(prset_only_model$model)
    dev.off()
    
    # Pathway Weight plot ----
    prset_only_coefs_df <- prset_only_model$model %>%
      coef(s = "lambda.1se") %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column() %>%
      setNames(c("pathway_id", "coef")) %>% 
      filter(!pathway_id %in% c("(Intercept)",
                                "age",
                                "age2",
                                "sex",
                                paste0("pc", 1:15))) %>%
      arrange(desc(abs(coef))) %>%
      mutate(index = row_number())
    
    prset_only_weights_plot <- ggplot(prset_only_coefs_df,
                                      aes(x = index,
                                          y = coef)) + 
      geom_point() +
      geom_col(width = 0.1) + 
      geom_hline(yintercept = 0) +
      geom_text_repel(data = prset_only_coefs_df %>%
                        slice_head(n = 10), 
                      aes(label = pathway_id)) +
      theme_bw() +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20)) +
      labs(x = "Index",
           y = "Pathway Coefficient")
    
    ggsave(paste0(results_path, "/", run_name, "_prset_only_weights_plot.png"), 
           prset_only_weights_plot,
           width = 10,
           height = 9)
    
    message("Running univariate models...\n")
    # Running univariate models ----
    univariate_models <- map(seq(formulas$univariate_formula), 
                             ~build_and_train_bin(formula = formulas$univariate_formula[.x],
                                                  .train_data = train_df,
                                                  .test_data = test_df,
                                                  .outcome_var = outcome_var,
                                                  .threshold = threshold,
                                                  .results_path = results_path,
                                                  .run_name = run_name,
                                                  file_name = paste0("_",
                                                                     formulas$univariate[.x],
                                                                     "_univariate")))
    
    # Running genome models ----
    message("Running genome models...\n")
    genome_models <- map(seq(formulas$genome_formula),
                         ~build_and_train_bin(formula = formulas$genome_formula[.x],
                                              .train_data = train_df,
                                              .test_data = test_df,
                                              .outcome_var = outcome_var,
                                              .threshold = threshold,
                                              .results_path = results_path,
                                              .run_name = run_name,
                                              file_name = paste0("_",
                                                                 formulas$univariate[.x],
                                                                 "_genome")))
    
    # Running pathway models ----
    message("Running prset models...")
    prset_models <- map(seq(nrow(formulas)),
                        ~build_and_lasso_bin(pred_vars = c(prset_pred_vars,
                                                           formulas$univariate %>%
                                                             as.vector() %>%
                                                             .[.x]),
                                             .train_data = train_df,
                                             .test_data = test_df,
                                             .outcome_var = outcome_var,
                                             .threshold = threshold,
                                             .results_path = results_path,
                                             .run_name = run_name,
                                             file_name = paste0("_",
                                                                formulas$univariate[.x],
                                                                "_prset")))
    # Calibration curves ----
    # Plotting calibration curves for prs, prset, 
    # clinical risk factor, clinical risk factor + prs,
    # clinical risk factor + prset
    convert_list_to_tibble <- function(model,
                                       .colnames) {
      test_preds <- model %>%
        map_dfc(~ as.data.frame(.x$test_pred))
      
      colnames(test_preds) <- .colnames
      
      test_preds
    }
    
    message("Calculating prediction stats ... \n")
    
    pred_df <- tibble(
      truth = test_df[[outcome_var]],
      prs_only = prs_only_model$test_pred,
      prset_only = prset_only_model$test_pred
    ) %>%
      # bind_cols(convert_list_to_tibble(univariate_models,
      #                                  formulas$univariate)) %>%
      bind_cols(convert_list_to_tibble(genome_models,
                                       paste0(formulas$univariate, "_prs"))) %>%
      bind_cols(convert_list_to_tibble(prset_models,
                                       paste0(formulas$univariate, "_prset"))) %>%
      pivot_longer(-truth, names_to = "model", values_to = "prediction") %>%
      mutate(prediction = as.numeric(prediction)) #coerce to vector
    
    #calibration curve
    cal_curve_prs <- pred_df %>%
      filter(model %in% c("prs_only",
                          "ascvd_prs",
                          "apoa_std_prs",
                          "apob_std_prs",
                          "ldl_std_prs",
                          "hdl_std_prs",
                          "trigly_std_prs",
                          "crp_std_prs")) %>% 
      group_by(model) %>%
      mutate(bin = ntile(prediction, 10)) %>% 
      # Bin prediction into 10ths
      group_by(bin,
               model) %>%
      mutate(n = n(), # Get ests and CIs
             bin_pred = mean(prediction),
             bin_prob = mean(as.numeric(truth)), 
             se = sqrt((bin_prob * (1 - bin_prob)) / n), 
             ul = bin_prob + 1.96 * se, 
             ll = bin_prob - 1.96 * se) %>% 
      ungroup()
    
    cal_curve_prs_plot <- ggplot(cal_curve_prs %>% 
                                   select(bin_pred,
                                          bin_prob,
                                          ll,
                                          ul,
                                          model) %>%
                                   distinct(),
                                 aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul, color = model)) +
      geom_abline() +
      geom_point() +
      geom_errorbar() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      
      # geom_point(position=position_dodge(width=0.1)) +
      # geom_errorbar(position=position_dodge(width=0.1)) +
      # geom_smooth(aes(x = prediction, y = as.numeric(truth)),
      #             method = "loess",
      #             se = FALSE) +
      labs(x = "Mean Predicted Risk",
           y = "Observed Event Rate",
           color = "Model") +
      theme_minimal() +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.title = element_text(size = 20), 
            legend.text = element_text(size = 18),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = c(0.15,0.8))
    
    cal_curve_prs_density_plot <- ggplot(cal_curve_prs, 
                                         aes(x = prediction,
                                             fill = model)) +
      geom_density(alpha = 0.4) +
      scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.1)) +
      xlab("Predicted Probability") +
      ylab("Density") +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.title = element_text(size = 20), 
            legend.text = element_text(size = 18),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none")
    
    cal_curve_prs_final_plot <- cal_curve_prs_plot / cal_curve_prs_density_plot + plot_layout(heights =  c(4,1))
    
    ggsave(paste0(results_path, "/", run_name, "_cal_curve_prs_plot.png"), 
           cal_curve_prs_final_plot,
           width = 10,
           height = 9)
    
    cal_curve_prset <- pred_df %>%
      filter(model %in% c("prset_only",
                          "ascvd_prset",
                          "apoa_std_prset",
                          "apob_std_prset",
                          "ldl_std_prset",
                          "hdl_std_prset",
                          "trigly_std_prset",
                          "crp_std_prset")) %>% 
      group_by(model) %>%
      mutate(bin = ntile(prediction, 10)) %>% 
      # Bin prediction into 10ths
      group_by(bin,
               model) %>%
      mutate(n = n(), # Get ests and CIs
             bin_pred = mean(prediction),
             bin_prob = mean(as.numeric(truth)), 
             se = sqrt((bin_prob * (1 - bin_prob)) / n), 
             ul = bin_prob + 1.96 * se, 
             ll = bin_prob - 1.96 * se) %>% 
      ungroup()
    
    cal_curve_prset_plot <- ggplot(cal_curve_prset %>% 
                                     select(bin_pred,
                                            bin_prob,
                                            ll,
                                            ul,
                                            model) %>%
                                     distinct(),
                                   aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul, color = model)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      geom_abline() +
      geom_point() +
      geom_errorbar() +
      # geom_point(position=position_dodge(width=0.1)) +
      # geom_errorbar(position=position_dodge(width=0.1)) +
      # geom_smooth(aes(x = prediction, y = as.numeric(truth)),
      #             method = "loess",
      #             se = FALSE) +
      labs(x = "Mean Predicted Risk",
           y = "Observed Event Rate",
           color = "Model") +
      theme_minimal() +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.title = element_text(size = 20), 
            legend.text = element_text(size = 18),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = c(0.15,0.8))
    
    cal_curve_prset_density_plot <- ggplot(cal_curve_prset, 
                                           aes(x = prediction,
                                               fill = model)) +
      geom_density(alpha = 0.4) +
      # geom_histogram(bins = 100, alpha = 0.4) + 
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      xlab("Predicted Probability") +
      ylab("Density") +
      theme(axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.title = element_text(size = 20), 
            legend.text = element_text(size = 18),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none")
    
    cal_curve_prset_final_plot <- cal_curve_prset_plot / cal_curve_prset_density_plot + plot_layout(heights =  c(4,1))
    
    ggsave(paste0(results_path, "/", run_name, "_cal_curve_prset_plot.png"), 
           cal_curve_prset_final_plot,
           width = 10,
           height = 9)
    
    # Calculating NRI ----
    map(seq(nrow(formulas)), function(x) {
      #variable vs. variable + prs
      sink(paste0(results_path, "/", run_name, "_",
                  formulas$univariate[[x]], "_nri_univariate_vs_prs.txt"))
      reclassification(data = data.frame(train_df), cOutcome = match(outcome_var, names(train_df)),
                       predrisk1 = predRisk(glm(as.formula(formulas$univariate_formula[[x]]),
                                                data = train_df,
                                                family = binomial)),
                       predrisk2 = predRisk(glm(as.formula(formulas$genome_formula[[x]]),
                                                data = train_df,
                                                family = binomial)),
                       cutoff = c(0, 0.075, 1)
      )
      sink(file = NULL)
      
      # variable vs. variable + prset
      sink(paste0(results_path, "/", run_name, "_",
                  formulas$univariate[[x]], "_nri_univariate_vs_prset.txt"))
      reclassification(data = data.frame(train_df), cOutcome = match(outcome_var, names(train_df)),
                       predrisk1 = predRisk(glm(as.formula(formulas$univariate_formula[[x]]),
                                                data = train_df,
                                                family = binomial)),
                       predrisk2 = predRisk(glm(as.formula(formulas$pathway_formula[[x]]),
                                                data = train_df,
                                                family = binomial)),
                       cutoff = c(0, 0.075, 1)
      )
      sink(file = NULL)
    })
    
    # Creating AUC plots ----
    message("Creating AUC plots...\n")
    
    univariate_models <- qread(paste0(results_path, "/", run_name,
                                      "_univariate_bin_results.qs"))
    genome_models <- qread(paste0(results_path, "/", run_name,
                                  "_genome_bin_results.qs"))
    prset_models <- qread(paste0(results_path, "/", run_name,
                                 "_prset_bin_results.qs"))
    
    univariate_auc_cis <- map(seq(univariate_models), ~qread(paste0(results_path, 
                                                                    "/",
                                                                    run_name,
                                                                    "_",
                                                                    formulas$univariate[.x],
                                                                    "_univariate_prs_only_boot.qs")))
    univariate_aucs <- tibble(auc = map(univariate_models, function(x) {
      x$test_auc 
    }) %>%
      flatten_dbl(),
    ll = map(univariate_auc_cis, function(x) {
      x[1] 
    }) %>%
      flatten_dbl(),
    ul = map(univariate_auc_cis, function(x) {
      x[3] 
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "1: Risk Factor + Age, Sex, PCs")
    
    genome_auc_cis <- map(seq(genome_models), ~qread(paste0(results_path, 
                                                            "/",
                                                            run_name,
                                                            "_",
                                                            formulas$univariate[.x],
                                                            "_genome_prs_only_boot.qs")))
    genome_aucs <- tibble(auc = map(genome_models, function(x) {
      x$test_auc
    }) %>%
      flatten_dbl(),
    ll = map(genome_auc_cis, function(x) {
      x[1] 
    }) %>%
      flatten_dbl(),
    ul = map(genome_auc_cis, function(x) {
      x[3] 
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "2: Model 1 + Genome-wide PRS")
    
    prset_auc_cis <- map(seq(prset_models), ~qread(paste0(results_path, 
                                                          "/",
                                                          run_name,
                                                          "_",
                                                          formulas$univariate[.x],
                                                          "_prset_prset_only_boot_lasso.qs")))
    prset_aucs <- tibble(auc = map(prset_models, function(x) {
      x$test_auc
    }) %>%
      flatten_dbl(),
    ll = map(prset_auc_cis, function(x) {
      x[1] 
    }) %>%
      flatten_dbl(),
    ul = map(prset_auc_cis, function(x) {
      x[3] 
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "3: Model 1 + Pathway PRS")
    
    all_aucs <- bind_rows(tibble(auc = prs_only_model$test_auc[[1]],
                                 variable = "prs_std",
                                 model = "Genome-wide PRS + Age, Sex, PCs",
                                 ll = prs_only_boot_ci[1],
                                 ul = prs_only_boot_ci[3]),
                                 # ll = prs_only_boot_ci[4][1] %>% unlist %>% .[2],
                                 # ul = prs_only_boot_ci[4][1] %>% unlist %>% .[3]),
                          tibble(auc = prset_only_model$test_auc[[1]],
                                 variable = "pathway_prs",
                                 model = "Pathway PRS + Age, Sex, PCs",
                                 # ll = prset_only_boot_ci[4][1] %>% unlist %>% .[2],
                                 ll = prset_only_boot_ci[1],
                                 ul = prset_only_boot_ci[3]),
                                 # ul = prset_only_boot_ci[4][1] %>% unlist %>% .[3]),
                          univariate_aucs,
                          genome_aucs,
                          prset_aucs) %>%
      mutate(variable = case_when(
        variable == "prs_std" ~ "PRS",
        variable == "pathway_prs" ~ "Pathway PRS",
        variable == "ascvd" ~ "ASCVD",
        variable == "apoa_std" ~ "ApoA",
        variable == "apob_std" ~ "ApoB",
        variable == "ldl_std" ~ "LDL",
        variable == "hdl_std" ~ "HDL",
        variable == "trigly_std" ~ "Triglycerides",
        variable == "crp_std" ~ "CRP",
        variable == "lpa_std" ~ "Lp(a)"
      )) %>%
      mutate(variable = factor(variable, levels = c("PRS",
                                                    "Pathway PRS",
                                                    "ASCVD",
                                                    "ApoA",
                                                    "ApoB",
                                                    "LDL",
                                                    "HDL",
                                                    "Triglycerides",
                                                    "Lp(a)",
                                                    "CRP")),
             model = factor(model, levels = c("Genome-wide PRS + Age, Sex, PCs",
                                              "Pathway PRS + Age, Sex, PCs",
                                              "1: Risk Factor + Age, Sex, PCs",
                                              "2: Model 1 + Genome-wide PRS",
                                              "3: Model 1 + Pathway PRS")))
    
    # subtype_barplot <- ggplot(subtype_df,
    #                           aes(x = !!sym(outcome_var),
    #                               y = n,
    #                               fill = !!sym(outcome_var))) + 
    #   geom_col() +
    #   geom_text(aes(label = round(prop, 2),
    #                 hjust = ifelse(prop < 0.2, -1, 1)),
    #             size = 3) +
    #   theme_bw() +
    #   theme(axis.title.x = element_blank(),
    #         axis.title.y = element_blank(),
    #         axis.text.x = element_blank(),
    #         panel.grid = element_blank(),
    #         legend.position = "none") +
    #   labs(y = "Count") +
    #   coord_flip()
    
    subtype_pie <- ggplot(subtype_df,
                          aes(x = "",
                              y = prop,
                              fill = !!sym(outcome_var))) + 
      geom_bar(stat = "identity", width = 1, color = "white") +
      # geom_text(aes(label = round(prop, 1),
      #               y = ypos),
      geom_text(aes(label = round(prop, 1)),
                position = position_stack(vjust = 0.5),
                size = 6) +
      # ggrepel::geom_label_repel(aes(label = !!sym(outcome_var),
      #                               y = ypos),
      #                           size = 3,
      #                           nudge_x = 0,
      #                           show.legend = FALSE) +
      coord_polar("y", start = 0) +
      theme_void() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank()) +
      guides(fill = guide_legend(nrow = 2)) +
      scale_fill_manual(
        values = subtype_pie_colors,
        limits = names(subtype_pie_colors)
      )
    subtype_pie
    
    #### For presentation ----
    ggplot(all_aucs, aes(x = forcats::fct_rev(variable), y = auc, fill = forcats::fct_rev(model))) +
      geom_col(position = position_dodge2(width = .9, preserve = "single")) +
      geom_linerange(aes(ymin = ll,
                         ymax = ul),
                     position = position_dodge2(width = .9, preserve = "single")) +
      coord_flip(ylim = c(0.5, 0.75)) + 
      theme_bw() +
      theme(plot.title = element_text(size = 22),
            axis.title.x = element_text(size = 18), 
            axis.title.y = element_blank(),
            panel.grid = element_blank(), 
            # legend.position = "bottom",
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16)) +
      # ggsci::scale_fill_npg() +
      scale_fill_manual(values = c(`Genome-wide PRS + Age, Sex, PCs` = "#70b7f3",
                                   `Pathway PRS + Age, Sex, PCs` = "#146bb4",
                                   `1: Risk Factor + Age, Sex, PCs` = "#f28982",
                                   `2: Model 1 + Genome-wide PRS` = "#ddaaf8",
                                   `3: Model 1 + Pathway PRS` = "#7c19b2")) +
      labs(x = "Variable",
           y = "Test AUC",
           fill = "Model",
           title = plot_title) +
      guides(fill = guide_legend(nrow = 2))
      # annotation_custom(
      #   ggplotGrob(subtype_pie),
      #   xmin = 0.5, xmax = 2, ymin = 0.8, ymax = 1
      # )
    #### ----
    all_auc_plot <- ggplot(all_aucs, aes(x = variable, y = auc, fill = model)) +
      geom_col(position = position_dodge2(width = .9, preserve = "single")) +
      # geom_text(aes(x = variable, y = auc, 
      #               label = glue("{format(round(auc, 2), nsmall = 2)}")),
      #           position = position_dodge2(width = .9, preserve = "single"),
      #           size = 5,
      #           stat = "identity",
      #           vjust = 0) +
      geom_linerange(aes(ymin = ll,
                         ymax = ul),
                     position = position_dodge2(width = .9, preserve = "single")) +
      theme_bw() +
      theme(plot.title = element_text(size = 22),
            axis.title.x = element_blank(), 
            axis.title.y = element_text(size = 18),
            panel.grid = element_blank(), 
            # legend.position = "bottom",
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16)) +
      coord_cartesian(ylim = c(0.5, 1)) +
      # ggsci::scale_fill_npg() +
      scale_fill_manual(values = c(`Genome-wide PRS + Age, Sex, PCs` = "#70b7f3",
                                   `Pathway PRS + Age, Sex, PCs` = "#146bb4",
                                   `1: Risk Factor + Age, Sex, PCs` = "#f28982",
                                   `2: Model 1 + Genome-wide PRS` = "#ddaaf8",
                                   `3: Model 1 + Pathway PRS` = "#7c19b2")) +
      labs(x = "Variable",
           y = "Test AUC",
           fill = "Model",
           title = plot_title) +
      guides(fill = guide_legend(nrow = 2)) +
      annotation_custom(
        ggplotGrob(subtype_pie),
        xmin = 0.5, xmax = 2, ymin = 0.8, ymax = 1
      )
    
    ggsave(paste0(results_path, "/", run_name, "_univariate_test_auc.png"), 
           all_auc_plot,
           width = 10,
           height = 9)
    
    # Creating ROC plots ----
    message("Creating ROC plots...\n")
    map(seq(length(univariate_models)), function(x) {
      switch(outcome_var,
             cad_bin = {
               pretty_outcome_var <- "CAD Case/Control"
             },
             ldl_bin = {
               pretty_outcome_var <- "LDL Subtype"
             }, 
             lpa_bin = {
               pretty_outcome_var <- "Lpa Subtype"
             },
             unstable_cad_bin = {
               pretty_outcome_var <- "Unstable Subtype"
             },
             occlusive_cad_bin = {
               pretty_outcome_var <- "Occlusive Subtype"
             },
             ascvd_cad_bin = {
               pretty_outcome_var <- "ASCVD Subtype"
             }
      )
      
      roc_plot <- ggroc(list(`1: Risk Factor + Age, Sex, PCs` = univariate_models[[x]]$test_roc,
                             `2: Model 1 + Genome-wide PRS` = genome_models[[x]]$test_roc,
                             `3: Model 1 + Pathway PRS` = prset_models[[x]]$test_roc),
                        size = 2,
                        legacy.axes = TRUE) +
        geom_abline(intercept = 0, linetype = 3) +
        ggtitle(glue("ROC Curves for \n {outcome_var} ~ {formulas$univariate[[x]]}")) + 
        theme_bw() +
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 20),
              legend.position = "none") +
        labs(color = "Model",
             x = "False Positive Rate",
             y = "True Positive Rate") +
        guides(color = guide_legend(nrow = 3))
      
      ggsave(paste0(results_path, "/", run_name, "_",
                    formulas$univariate[[x]], "_univariate_test_roc.png"), 
             roc_plot)
    })
    
    # Testing differences in AUC ----
    # genome_test <- map(seq(length(univariate_models)),
    #                    ~ lrtest(univariate_models[[.x]]$model, 
    #                             genome_models[[.x]]$model))
    # 
    # prset_test <- map(seq(length(univariate_models)),
    #                   ~ lrtest(univariate_models[[.x]]$test_roc, 
    #                                    prset_models[[.x]]$test_roc,
    #                                    boot.n = 5000,
    #                                    parallel = TRUE))
    # 
    # genome_delta <- map(seq(length(genome_test)), 
    #                     ~ genome_test[[.x]]$estimate[2] - genome_test[[.x]]$estimate[1]) %>%
    #   unlist()
    # 
    # genome_delta_p <- map(seq(length(genome_test)), 
    #                       ~ genome_test[[.x]]$`Pr(>Chisq)`) %>%
    #   unlist()
    # 
    # prset_delta <- map(seq(length(prset_test)), 
    #                    ~ prset_test[[.x]]$estimate[2] - prset_test[[.x]]$estimate[1]) %>%
    #   unlist()
    test_auc_difference(.results_path = results_path,
                        .run_name = run_name)
    
    # Saving results ----
    qsave(univariate_models, paste0(results_path, "/", run_name,
                                    "_univariate_bin_results.qs"))
    qsave(genome_models, paste0(results_path, "/", run_name,
                                "_genome_bin_results.qs"))
    qsave(prset_models, paste0(results_path, "/", run_name,
                               "_prset_bin_results.qs"))
  } else {
    
    # Univariate models ----
    
    message("Running PRS-only model")
    if(str_detect(run_name, "ascvd_cad")) {
      prs_formula <- glue("{outcome_var} ~ prs_std + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15")
    } else {
      prs_formula <- glue("{outcome_var} ~ prs_std + age + age2 + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + pc11 + pc12 + pc13 + pc14 + pc15")
    }
    
    prs_only_model <- build_and_train_cont(formula = prs_formula,
                                           .train_data = train_df,
                                           .test_data = test_df,
                                           .outcome_var = outcome_var)
    qsave(prs_only_model, paste0(results_path, "/", run_name,
                                 "_prs_only_results.qs"))
    
    message("Running PRSet-only model")
    
    if(str_detect(run_name, "ascvd_cad")) {
      prset_pred_vars <- c(paste0("pc", 1:15), prset_colnames)
    } else {
      prset_pred_vars <- c("age", "age2", "sex",
                           paste0("pc", 1:15), prset_colnames)
    }
    
    prset_only_model <- build_and_lasso_cont(pred_vars = prset_pred_vars,
                                             .train_data = train_df,
                                             .test_data = test_df,
                                             .outcome_var = outcome_var)
    qsave(prset_only_model, paste0(results_path, "/", run_name,
                                   "_prset_only_results.qs"))
    
    message("Running univariate models...\n")
    univariate_models <- map(formulas$univariate_formula, 
                             ~build_and_train_cont(formula = .x,
                                                   .train_data = train_df,
                                                   .test_data = test_df,
                                                   .outcome_var = outcome_var))
    
    # Genome models ----
    message("Running genome models...\n")
    genome_models <- map(formulas$genome_formula,
                         ~build_and_train_cont(formula = .x,
                                               .train_data = train_df,
                                               .test_data = test_df,
                                               .outcome_var = outcome_var))
    
    # Pathway models ----
    message("Running prset models...\n")
    prset_models <- map(seq(nrow(formulas)),
                        ~build_and_lasso_cont(pred_vars = c(prset_pred_vars,
                                                            formulas$univariate %>%
                                                              as.vector() %>%
                                                              .[.x]),
                                              .train_data = train_df,
                                              .test_data = test_df,
                                              .outcome_var = outcome_var))
    
    # Creating R2 plots ----
    message("Creating R2 plots...\n")
    univariate_r2 <- tibble(r2 = map(univariate_models, function(x) {
      x$performance$test_cor
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "1: Risk Factor + Age, Sex, PCs")
    
    genome_r2 <- tibble(r2 = map(genome_models, function(x) {
      x$performance$test_cor
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "2: Model 1 + Genome-wide PRS")
    
    prset_r2 <- tibble(r2 = map(prset_models, function(x) {
      x$performance$test_cor
    }) %>%
      flatten_dbl(),
    variable = formulas$univariate,
    model = "3: Model 1 + Pathway PRS")
    
    all_r2 <- bind_rows(tibble(r2 = prs_only_model$performance$test_cor,
                               variable = "prs_std",
                               model = "Genome-wide PRS + Age, Sex, PCs"),
                        tibble(r2 = prset_only_model$performance$test_cor,
                               variable = "pathway_prs",
                               model = "Pathway PRS + Age, Sex, PCs"),
                        univariate_r2,
                        genome_r2,
                        prset_r2) %>%
      mutate(variable = case_when(
        variable == "prs_std" ~ "PRS",
        variable == "pathway_prs" ~ "Pathway PRS",
        variable == "ascvd" ~ "ASCVD",
        variable == "apoa_std" ~ "ApoA",
        variable == "apob_std" ~ "ApoB",
        variable == "ldl_std" ~ "LDL",
        variable == "hdl_std" ~ "HDL",
        variable == "trigly_std" ~ "Triglycerides",
        variable == "crp_std" ~ "CRP"
      )) %>%
      mutate(variable = factor(variable, levels=c("PRS",
                                                  "Pathway PRS",
                                                  "ASCVD",
                                                  "ApoA",
                                                  "ApoB",
                                                  "LDL",
                                                  "HDL",
                                                  "Triglycerides",
                                                  "CRP")))
    
    all_r2_plot <- ggplot(all_r2, aes(x = variable, y = r2, fill = model)) +
      geom_col(position = position_dodge2(width = .9, preserve = "single")) +
      # geom_text(aes(x = variable, y = r2, 
      #               label = glue("{format(round(r2, 2), nsmall = 2)}")),
      #           position = position_dodge2(width = .9, preserve = "single"),
      #           size = 5,
      #           stat = "identity",
      #           vjust = 0) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            panel.grid = element_blank(), 
            # legend.position ="bottom",
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
            axis.text.y = element_text(size = 16)) +
      # scale_y_continuous(limits = c(0, 0.1), expand = c(0, 0)) +
      ggsci::scale_fill_npg() +
      labs(x = "Variable",
           y = bquote(R^2),
           fill = "Model",
           # title = gsub("cad_subtype/output_data/(.*?)_bin_mod1\\.csv", 
           #              '\\1', pheno_df_path),
           title = plot_title) +
      guides(fill = guide_legend(nrow = 2))
    
    ggsave(paste0(results_path, "/", run_name, "_univariate_test_r2.png"), 
           all_r2_plot,
           width = 10,
           height = 9)
    
    # Saving results ----
    qsave(univariate_models, paste0(results_path, "/", run_name,
                                    "_univariate_results.qs"))
    qsave(genome_models, paste0(results_path, "/", run_name,
                                "_genome_results.qs"))
    qsave(prset_models, paste0(results_path, "/", run_name,
                               "_prset_results.qs"))
  }
  
  return(list(univariate_models = univariate_models,
              genome_models = genome_models,
              prset_models = prset_models))
}

res <- run_analysis(original_ukb_path = argv$original_ukb_path,
                    pheno_df_path = argv$pheno_df_path,
                    results_path = argv$results_path,
                    train_prop = argv$train_prop,
                    run_name = argv$run_name,
                    prsice_run_name = argv$prsice_run_name,
                    prset_run_name = argv$prset_run_name,
                    outcome_var = argv$outcome_var,
                    seed = argv$seed,
                    threshold = argv$threshold)

# # residualized, all pop
# ldl_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
#                            pheno_df_path = "cad_subtype/output_data/all_ldl_bin_mod1.csv",
#                            results_path = "cad_subtype/result_data/ldl_subtype",
#                            train_prop = 0.8,
#                            run_name = "c4d_all_ukb_cad_ldlsub",
#                            prsice_run_name = "c4d_prsice_all_ukb_cad_ldlsub",
#                            prset_run_name = "c4d_prset_all_ukb_cad_ldlsub",
#                            outcome_var = "pheno",
#                            seed = 47,
#                            threshold = 0.5)
# lpa_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
#                            pheno_df_path = "cad_subtype/output_data/all_lpa_bin_mod1.csv",
#                            results_path = "cad_subtype/result_data/lpa_subtype",
#                            train_prop = 0.8,
#                            run_name = "all_ukb_cad_lpasub",
#                            prsice_run_name = "prsice_all_ukb_cad_lpasub",
#                            prset_run_name = "prset_all_ukb_cad_lpasub",
#                            outcome_var = "pheno",
#                            seed = 47,
#                            threshold = 0.5)
# unstable_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
#                                 pheno_df_path = "cad_subtype/output_data/all_unstable_cad_bin_mod1.csv",
#                                 results_path = "cad_subtype/result_data/unstable_subtype",
#                                 train_prop = 0.8,
#                                 run_name = "all_ukb_cad_unstablesub",
#                                 prsice_run_name = "prsice_all_ukb_cad_unstablesub",
#                                 prset_run_name = "prset_all_ukb_cad_unstablesub",
#                                 outcome_var = "pheno",
#                                 seed = 47,
#                                 threshold = 0.5,
#                                 env_sensitivity = TRUE)
# occlusive_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
#                                  pheno_df_path = "cad_subtype/output_data/all_occlusive_cad_bin_mod1.csv",
#                                  results_path = "cad_subtype/result_data/occlusive_subtype",
#                                  train_prop = 0.8,
#                                  run_name = "all_ukb_cad_occlusivesub",
#                                  prsice_run_name = "prsice_all_ukb_cad_occlusivesub",
#                                  prset_run_name = "prset_all_ukb_cad_occlusivesub",
#                                  outcome_var = "pheno",
#                                  seed = 47,
#                                  threshold = 0.5)
# ascvd_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
#                              pheno_df_path = "cad_subtype/output_data/all_ascvd_cad_bin_mod1.csv",
#                              results_path = "cad_subtype/result_data/ascvd_subtype",
#                              train_prop = 0.8,
#                              run_name = "all_ukb_cad_ascvdsub",
#                              prsice_run_name = "prsice_all_ukb_cad_ascvdsub",
#                              prset_run_name = "prset_all_ukb_cad_ascvdsub",
#                              outcome_var = "pheno",
#                              seed = 47,
#                              threshold = 0.5)
#  
# # bin, all pop
# start <- Sys.time()
cad_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                               pheno_df_path = "cad_subtype/output_data/all_cad_bin_mod1.csv",
                               results_path = "cad_subtype/result_data/cad_bin",
                               train_prop = 0.8,
                               run_name = "c4d_all_ukb_cad_bin",
                               prsice_run_name = "c4d_prsice_all_ukb_cad_bin",
                               prset_run_name = "c4d_prset_all_ukb_cad_bin",
                               outcome_var = "cad_bin",
                               seed = 47,
                               threshold = 0.5)
# end <- Sys.time()
# end - start
# 
ldl_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                               pheno_df_path = "cad_subtype/output_data/all_ldl_bin_mod1.csv",
                               results_path = "cad_subtype/result_data/ldl_subtype",
                               train_prop = 0.8,
                               run_name = "c4d_all_ukb_cad_ldlsub_bin",
                               prsice_run_name = "c4d_prsice_all_ukb_cad_ldlsub_bin",
                               prset_run_name = "c4d_prset_all_ukb_cad_ldlsub_bin",
                               outcome_var = "ldl_bin",
                               seed = 47,
                               threshold = 0.5)
lpa_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                               pheno_df_path = "cad_subtype/output_data/all_lpa_bin_mod1.csv",
                               results_path = "cad_subtype/result_data/lpa_subtype",
                               train_prop = 0.8,
                               run_name = "c4d_all_ukb_cad_lpasub_bin",
                               prsice_run_name = "c4d_prsice_all_ukb_cad_lpasub_bin",
                               prset_run_name = "c4d_prset_all_ukb_cad_lpasub_bin",
                               outcome_var = "lpa_bin",
                               seed = 47,
                               threshold = 0.5)

#sensitivity analysis - no lpa gene
lpa_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                               pheno_df_path = "cad_subtype/output_data/all_lpa_bin_mod1.csv",
                               results_path = "cad_subtype/result_data/lpa_subtype",
                               train_prop = 0.8,
                               run_name = "c4d_no_lpa_all_ukb_cad_lpasub_bin",
                               prsice_run_name = "c4d_no_lpa_prsice_all_ukb_cad_lpasub_bin",
                               prset_run_name = "c4d_no_lpa_prset_all_ukb_cad_lpasub_bin",
                               outcome_var = "lpa_bin",
                               seed = 47,
                               threshold = 0.5)

stemi_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                                 pheno_df_path = "cad_subtype/output_data/all_stemi_cad_bin_mod1.csv",
                                 results_path = "cad_subtype/result_data/stemi_subtype",
                                 train_prop = 0.8,
                                 run_name = "c4d_all_ukb_cad_stemisub_bin",
                                 prsice_run_name = "c4d_prsice_all_ukb_cad_stemisub_bin",
                                 prset_run_name = "c4d_prset_all_ukb_cad_stemisub_bin",
                                 outcome_var = "stemi_cad_bin",
                                 seed = 47,
                                 threshold = 0.5)
ascvd_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                                 pheno_df_path = "cad_subtype/output_data/all_ascvd_cad_bin_mod1.csv",
                                 results_path = "cad_subtype/result_data/ascvd_subtype",
                                 train_prop = 0.8,
                                 run_name = "c4d_all_ukb_cad_ascvdsub_bin",
                                 prsice_run_name = "c4d_prsice_all_ukb_cad_ascvdsub_bin",
                                 prset_run_name = "c4d_prset_all_ukb_cad_ascvdsub_bin",
                                 outcome_var = "ascvd_cad_bin",
                                 seed = 47,
                                 threshold = 0.5)
occlusive_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                                     pheno_df_path = "cad_subtype/output_data/all_occlusive_cad_bin_mod1.csv",
                                     results_path = "cad_subtype/result_data/occlusive_subtype",
                                     train_prop = 0.8,
                                     run_name = "c4d_all_ukb_cad_occlusivesub_bin",
                                     prsice_run_name = "c4d_prsice_all_ukb_cad_occlusivesub_bin",
                                     prset_run_name = "c4d_prset_all_ukb_cad_occlusivesub_bin",
                                     outcome_var = "occlusive_cad_bin",
                                     seed = 47,
                                     threshold = 0.5)
unstable_bin_allpop <- run_analysis(original_ukb_path = "cad_subtype/output_data/all_clean_original_ukb.csv",
                                    pheno_df_path = "cad_subtype/output_data/all_unstable_cad_bin_mod1.csv",
                                    results_path = "cad_subtype/result_data/unstable_subtype",
                                    train_prop = 0.8,
                                    run_name = "c4d_all_ukb_cad_unstablesub_bin",
                                    prsice_run_name = "c4d_prsice_all_ukb_cad_unstablesub_bin",
                                    prset_run_name = "c4d_prset_all_ukb_cad_unstablesub_bin",
                                    outcome_var = "unstable_cad_bin",
                                    seed = 47,
                                    threshold = 0.5)

