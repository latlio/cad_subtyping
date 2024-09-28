# Extract pathways, SNPs, and betas for validation in All of Us
# Beatrice is helping me with this

library(tidyverse)
library(qs)
library(data.table)
source("cad_subtype/R/common_prset_fcts.R")

extract_info <- function(prset_only_res_path,
                         prset_res_path,
                         prset_snp_path,
                         save_name) {
  prset_only_res <- qread(prset_only_res_path)
  prset_res <- qread(prset_res_path)
  prset_snp <- fread(prset_snp_path)
  
  if(str_detect(prset_res_path, "ldlsub")) {
    univariate <- c("ascvd", "apoa_std", "trigly_std", "hdl_std", "lpa_std", "crp_std")
  } else if(str_detect(prset_res_path, "lpasub")) {
    univariate <- c("aascvd", "apoa_std", "apob_std", "trigly_std", "ldl_std", "hdl_std",
                    "crp_std")
  } else if(str_detect(prset_res_path, "ascvdsub")) {
    univariate <- c("apoa_std", "apob_std", "trigly_std", "ldl_std", "hdl_std",
                    "crp_std")
  } else {
    univariate <- c("ascvd", "apoa_std", "apob_std", "trigly_std", "ldl_std", "hdl_std", "lpa_std",
                    "crp_std")
  }
  
  #df of selected pathways (pathway id, beta coef) ----
  #select best lambda.1se: index[2]
  
  #prset-only model
  prset_only_pathways <- prset_only_res$model$glmnet.fit$beta[,prset_only_res$model$index[2]] %>%
    as.matrix() %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    setNames(c("pathway_id", "coef")) %>% 
    filter(!is.na(coef)) %>%
    filter(coef != 0) %>%
    arrange(desc(coef)) %>% 
    mutate(pathway_id = case_when(
      str_detect(pathway_id, "^go_") ~ toupper(pathway_id) %>% str_replace(., "_", ":"),
      str_detect(pathway_id, "^mp_") ~ toupper(pathway_id) %>% str_replace(., "_", ":"),
      TRUE ~ toupper(pathway_id)
    ))
  
  #map across each risk factor
  pathways <- map(seq(univariate), ~ prset_res[[.x]]$model$glmnet.fit$beta[,prset_res[[.x]]$model$index[2]] %>%
                    as.matrix() %>%
                    as.data.frame() %>%
                    rownames_to_column() %>%
                    setNames(c("pathway_id", "coef")) %>% 
                    filter(!is.na(coef)) %>%
                    filter(coef != 0) %>%
                    arrange(desc(coef)) %>% 
                    mutate(pathway_id = case_when(
                      str_detect(pathway_id, "^go_") ~ toupper(pathway_id) %>% str_replace(., "_", ":"),
                      str_detect(pathway_id, "^mp_") ~ toupper(pathway_id) %>% str_replace(., "_", ":"),
                      TRUE ~ toupper(pathway_id)
                    ))
  )
  names(pathways) <- univariate
  
  pathways2 <- append(pathways, list("prset-only" = prset_only_pathways))
  
  qsave(pathways2, paste0(save_name, ".qs"))
  
  # df of snps of selected pathways ----
  snps_in_pathways <- map(seq(pathways2), ~ prset_snp %>%
                            select(c(CHR, SNP, BP, pathways2[[.x]] %>%
                                       filter(!pathway_id %in% toupper(univariate),
                                              !pathway_id %in% c("AGE",
                                                                 "SEX",
                                                                 "AGE2",
                                                                 paste0("PC", 1:15))) %>%
                                       pull(pathway_id))))
  names(snps_in_pathways) <- univariate
  
  qsave(snps_in_pathways, paste0(save_name, "_snps.qs"))
}

retrieve_snps_for_pathways <- function(prset_lasso_res,
                                       prset_snp_path,
                                       out_path) {
  lasso_selected_pathways <- save_model_coef(prset_lasso_res$model) %>%
    filter(go_id != "(Intercept)") %>%
    pull(go_id)
  
  prset_snp <- read_tsv(prset_snp_path)
  
  out <- map(lasso_selected_pathways, ~ prset_snp %>%
               select(c(SNP, lasso_selected_pathways)) %>%
               filter(!!sym(.x) == 1) %>%
               select(SNP) %>%
               mutate(path = .x)) %>%
    bind_rows()
  
  write_csv(out, out_path)
}

#have to run cad
extract_info("cad_subtype/result_data/cad_bin/c4d_all_ukb_cad_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/cad_bin/c4d_all_ukb_cad_bin_prset_bin_results.qs",
             "cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin.snp",
             "cad_subtype/result_data/cad_bin/c4d_cad_pathways")

extract_info("cad_subtype/result_data/ascvd_subtype/c4d_all_ukb_cad_ascvdsub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/ascvd_subtype/c4d_all_ukb_cad_ascvdsub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin.snp",
             "cad_subtype/result_data/ascvd_subtype/c4d_ascvd_pathways")

extract_info("cad_subtype/result_data/ldl_subtype/c4d_all_ukb_cad_ldlsub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/ldl_subtype/c4d_all_ukb_cad_ldlsub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.snp",
             "cad_subtype/result_data/ldl_subtype/c4d_ldl_pathways")

extract_info("cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/lpa_subtype/c4d_all_ukb_cad_lpasub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.snp",
             "cad_subtype/result_data/lpa_subtype/c4d_lpa_pathways")

extract_info("cad_subtype/result_data/lpa_subtype/c4d_no_lpa_all_ukb_cad_lpasub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_all_ukb_cad_lpasub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_prset_all_ukb_cad_lpasub_bin.snp",
             "cad_subtype/result_data/lpa_subtype/c4d_sensitivity_lpa_pathways")

extract_info("cad_subtype/result_data/occlusive_subtype/c4d_all_ukb_cad_occlusivesub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/occlusive_subtype/c4d_all_ukb_cad_occlusivesub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.snp",
             "cad_subtype/result_data/occlusive_subtype/c4d_occlusive_pathways")

extract_info("cad_subtype/result_data/unstable_subtype/c4d_all_ukb_cad_unstablesub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/unstable_subtype/c4d_all_ukb_cad_unstablesub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin.snp",
             "cad_subtype/result_data/unstable_subtype/c4d_unstable_pathways")

extract_info("cad_subtype/result_data/stemi_subtype/c4d_all_ukb_cad_stemisub_bin_prset_only_bin_results.qs",
             "cad_subtype/result_data/stemi_subtype/c4d_all_ukb_cad_stemisub_bin_prset_bin_results.qs",
             "cad_subtype/result_data/stemi_subtype/c4d_prset_all_ukb_cad_stemisub_bin.snp",
             "cad_subtype/result_data/stemi_subtype/c4d_stemi_pathways")

# retrieve_snps_for_pathways("cad_subtype/result_data/lpa_subtype/prset_all_ukb_cad_lpasub_bin_results.Rds",
#                            "cad_subtype/result_data/lpa_subtype/all_ukb_cad_lpasub_bin_prset_bin_results.Rds",
#                            "cad_subtype/result_data/lpa_subtype/prset_all_ukb_cad_lpasub_bin.snp",
#                            "lpa_subtype/lpa_pathways")


# find identity of snps-genes
lpa_snps <- qs::qread("/Users/lathanliou/Desktop/Academic/Sinai/Research/OReilly/Code/cad_subtype/result_data/lpa_subtype/c4d_lpa_pathways_snps.qs")
path1_snps <- lpa_snps[[1]] %>% filter(`GO:0001968` != 0) %>% pull(SNP)
library("rsnps")
path1_snps_db <- ncbi_snp_query(path1_snps)

path2_snps <- lpa_snps[[1]] %>% filter(`GO:0034185` != 0) %>% pull(SNP)
path2_snps_db <- ncbi_snp_query(path2_snps)            

path3_snps <- lpa_snps[[1]] %>% filter(`GO:0004252` != 0) %>% pull(SNP)
path3_snps_db <- ncbi_snp_query(path3_snps) 
