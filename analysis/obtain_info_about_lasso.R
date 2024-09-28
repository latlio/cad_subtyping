library(tidyverse)
library(qs)
library(glmnet) # required for coef()
library(here)

# find how many snps are in a given model

find_number_of_pathways <- function(prset_models_path) {
  prset_models <- qread(prset_models_path)
  
  map(seq(length(prset_models)), ~prset_models[[.x]]$model %>% 
        coef() %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        setNames(c("variable", "coef")) %>%
        filter(coef != 0,
               !variable %in% c("age", "age2", "sex", paste0("pc", 1:15),
                                "ascvd", "apoa_std", "apob_std", "trigly_std",
                                "ldl_std", "hdl_std", "crp_std", "(Intercept)")) %>%
        nrow()) %>%
    unlist()
}

find_number_of_enriched_pathways <- function(prset_summary_path,
                                             prset_models_path) {
  prset_summary <- fread(here(prset_summary_path)) %>%
    janitor::clean_names() %>%
    filter(!str_detect(set, "MP:")) %>%
    mutate(term = case_when(
      str_detect(set, "GO:") ~ GOfuncR::get_names(set)$go_name,
      # str_detect(set, "MP:") ~ mp_term,
      TRUE ~ set)) %>%
    mutate(term = str_replace_all(term, "_", " ")) %>%
    mutate(term = str_to_title(term)) %>%
    relocate(term, .after = set) %>%
    filter(!is.na(term))
  
  prset_models <- qread(prset_models_path)
  
  map(seq(length(prset_models)), ~prset_models[[.x]]$model %>% 
        coef() %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        setNames(c("variable", "coef")) %>%
        filter(coef != 0,
               !variable %in% c("age", "age2", "sex", paste0("pc", 1:15),
                                "ascvd", "apoa_std", "apob_std", "trigly_std",
                                "ldl_std", "hdl_std", "crp_std", "(Intercept)")) %>%
        mutate(variable = case_when(
          str_detect(variable, "^go") ~ toupper(variable) %>% str_replace(., "_", ":"),
          str_detect(variable, "^reactome") ~ toupper(variable),
          str_detect(variable, "^pid") ~ toupper(variable),
          str_detect(variable, "^kegg") ~ toupper(variable),
          str_detect(variable, "^biocarta") ~ toupper(variable),
          str_detect(variable, "^mp") ~ toupper(variable) %>% str_replace(., "_", ":")
        )) %>%
        left_join(prset_summary,
                  by = c("variable" = "set")) %>%
        filter(competitive_p < 0.05) %>%
        nrow()) %>%
    unlist()
}

extract_identities_of_pathways <- function(prset_models_path) {
  prset_models <- qread(prset_models_path)
  
  map(seq(length(prset_models)), ~prset_models[[.x]]$model %>% 
        coef() %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        setNames(c("variable", "coef")) %>%
        filter(coef != 0,
               !variable %in% c("age", "age2", "sex", paste0("pc", 1:15),
                                "ascvd", "apoa_std", "apob_std", "trigly_std",
                                "ldl_std", "hdl_std", "crp_std", "(Intercept)")))
}

turn_zeros_to_nas_and_remove <- function(x) {
  x[x == 0] <- NA
  
  x_complete <- x[complete.cases(x),]
  return(x_complete)
}

nrow_safely <- safely(nrow)

find_number_of_snps_in_lasso <- function(pathway_snps_path,
                                         prset_summary_path,
                                         prset_models_path) {
  pathway_snps <- qread(pathway_snps_path)
  
  #remove prset_only results (last index)
  pathway_snps_clean <- pathway_snps[-length(pathway_snps)]
  
  prset_summary <- fread(here(prset_summary_path)) %>%
    janitor::clean_names() %>%
    filter(!str_detect(set, "MP:")) %>%
    mutate(term = case_when(
      str_detect(set, "GO:") ~ GOfuncR::get_names(set)$go_name,
      # str_detect(set, "MP:") ~ mp_term,
      TRUE ~ set)) %>%
    mutate(term = str_replace_all(term, "_", " ")) %>%
    mutate(term = str_to_title(term)) %>%
    relocate(term, .after = set) %>%
    filter(!is.na(term))
  
  prset_models <- qread(prset_models_path)
  
  pathways <- map(seq(length(prset_models)), ~prset_models[[.x]]$model %>% 
                    coef() %>%
                    as.matrix() %>%
                    as.data.frame() %>%
                    rownames_to_column() %>%
                    setNames(c("variable", "coef")) %>%
                    filter(coef != 0,
                           !variable %in% c("age", "age2", "sex", paste0("pc", 1:15),
                                            "ascvd", "apoa_std", "apob_std", "trigly_std",
                                            "ldl_std", "hdl_std", "crp_std", "(Intercept)")) %>%
                    mutate(variable = case_when(
                      str_detect(variable, "^go") ~ toupper(variable) %>% str_replace(., "_", ":"),
                      str_detect(variable, "^reactome") ~ toupper(variable),
                      str_detect(variable, "^pid") ~ toupper(variable),
                      str_detect(variable, "^kegg") ~ toupper(variable),
                      str_detect(variable, "^biocarta") ~ toupper(variable),
                      str_detect(variable, "^mp") ~ toupper(variable) %>% str_replace(., "_", ":")
                    )) %>%
                    left_join(prset_summary,
                              by = c("variable" = "set")) %>%
                    filter(competitive_p < 0.05))
  
  filtered_snps <- map(seq(length(pathway_snps_clean)), ~ pathway_snps_clean[[.x]] %>%
                select(all_of(pathways[[.x]]$variable)) %>% 
                filter_all(., any_vars(. != 0)) %>%
                nrow_safely() %>%
                  .$result) %>%
    unlist()
  filtered_snps
}

find_number_of_pathways("cad_subtype/result_data/ascvd_subtype/c4d_all_ukb_cad_ascvdsub_bin_prset_bin_results.qs")
find_number_of_pathways("cad_subtype/result_data/ldl_subtype/c4d_all_ukb_cad_ldlsub_bin_prset_bin_results.qs")
find_number_of_pathways("cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
find_number_of_pathways("cad_subtype/result_data/occlusive_subtype/c4d_all_ukb_cad_occlusivesub_bin_prset_bin_results.qs")
find_number_of_pathways("cad_subtype/result_data/unstable_subtype/c4d_all_ukb_cad_unstablesub_bin_prset_bin_results.qs")
find_number_of_pathways("cad_subtype/result_data/stemi_subtype/c4d_all_ukb_cad_stemisub_bin_prset_bin_results.qs")

find_number_of_enriched_pathways("cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin.summary",
                                 "cad_subtype/result_data/ascvd_subtype/c4d_all_ukb_cad_ascvdsub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.summary",
                                 "cad_subtype/result_data/ldl_subtype/c4d_all_ukb_cad_ldlsub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.summary",
                                 "cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.summary",
                                 "cad_subtype/result_data/occlusive_subtype/c4d_all_ukb_cad_occlusivesub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin.summary",
                                 "cad_subtype/result_data/unstable_subtype/c4d_all_ukb_cad_unstablesub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/stemi_subtype/c4d_prset_all_ukb_cad_stemisub_bin.summary",
                                 "cad_subtype/result_data/stemi_subtype/c4d_all_ukb_cad_stemisub_bin_prset_bin_results.qs")


lpa_pathways <- extract_identities_of_pathways("cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_bin_results.qs")

# *_pathways_snps.qs comes from extract_info.R

find_number_of_snps_in_lasso("cad_subtype/result_data/ascvd_subtype/c4d_ascvd_pathways_snps.qs",
                             "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin.summary",
                             "cad_subtype/result_data/ascvd_subtype/c4d_all_ukb_cad_ascvdsub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/ldl_subtype/c4d_ldl_pathways_snps.qs",
                             "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.summary",
                             "cad_subtype/result_data/ldl_subtype/c4d_all_ukb_cad_ldlsub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/lpa_subtype/c4d_lpa_pathways_snps.qs",
                             "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.summary",
                             "cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/occlusive_subtype/c4d_occlusive_pathways_snps.qs",
                             "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.summary",
                             "cad_subtype/result_data/occlusive_subtype/c4d_all_ukb_cad_occlusivesub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/unstable_subtype/c4d_unstable_pathways_snps.qs",
                             "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin.summary",
                             "cad_subtype/result_data/unstable_subtype/c4d_all_ukb_cad_unstablesub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/stemi_subtype/c4d_stemi_pathways_snps.qs",
                             "cad_subtype/result_data/stemi_subtype/c4d_prset_all_ukb_cad_stemisub_bin.summary",
                             "cad_subtype/result_data/stemi_subtype/c4d_all_ukb_cad_stemisub_bin_prset_bin_results.qs")
stemi
# sensitivity lpa
find_number_of_pathways("cad_subtype/result_data/lpa_subtype/c4d_no_lpa_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
find_number_of_snps_in_lasso("cad_subtype/result_data/lpa_subtype/c4d_sensitivity_lpa_pathways_snps.qs",
                             "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_prset_all_ukb_cad_lpasub_bin.summary",
                             "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
find_number_of_enriched_pathways("cad_subtype/result_data/lpa_subtype/c4d_no_lpa_prset_all_ukb_cad_lpasub_bin.summary",
                                 "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_all_ukb_cad_lpasub_bin_prset_bin_results.qs")
