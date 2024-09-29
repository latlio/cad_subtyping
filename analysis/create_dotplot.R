# Create dotplots of pathways

library(tidyverse)
library(data.table)
library(GOfuncR)
library(here)
library(janitor)
library(scales)
library(glmnet)
# library(gghighlight)
library(tidytext)

source("cad_subtype/R/common_prset_fcts.R")
# source("common_prset_fcts.R")
create_dotplot <- function(prset_summary_path,
                           lasso_path,
                           filename_top20,
                           filename_lasso) {
  #using read.csv to automatically skip blank lines
  # mp_datadict <- read.csv("cad_subtype/MP_datadict.csv")
  prset_summary <- read_tsv(here(prset_summary_path)) %>%
    janitor::clean_names() %>%
    # left_join(mp_datadict, by = c("set" = "mp_id")) %>%
    filter(!str_detect(set, "MP:")) %>%
    mutate(term = case_when(
      str_detect(set, "GO:") ~ GOfuncR::get_names(set)$go_name,
      # str_detect(set, "MP:") ~ mp_term,
      TRUE ~ set)) %>%
    mutate(term = str_replace_all(term, "_", " ")) %>%
    mutate(term = str_to_title(term)) %>%
    relocate(term, .after = set) %>%
    filter(!is.na(term))
  
  prset_lasso_res <- readRDS(here(lasso_path))
  
  lasso_selected_pathways <- save_model_coef(prset_lasso_res$model) %>%
    filter(go_id != "(Intercept)") %>%
    pull(go_id)
  
  # all pathways
  top20_dotplot <- prset_summary %>%
    filter(competitive_p <= 0.05) %>%
    arrange(desc(prs_r2)) %>%
    slice_head(n = 20) %>%
    ggplot(aes(x = fct_reorder(term, prs_r2), 
               y = prs_r2, size = num_snp, color = competitive_p)) +
    geom_point() +
    scale_x_discrete(labels = label_wrap(20)) +
    scale_size(range = c(5, 20)) +
    coord_flip() +
    theme_bw() +
    labs(y = bquote(R^2),
         size = "Number of \n SNPs",
         color = "Competitive \n p-value") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 16))
  
  ggsave(filename_top20, top20_dotplot,
         width = 10, height = 9)
  
  # lasso-selected
  lasso_dotplot <- prset_summary %>%
    filter(set %in% lasso_selected_pathways,
           competitive_p <= 0.25) %>%
    arrange(competitive_p) %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = fct_reorder(term, prs_r2), 
               y = prs_r2, size = num_snp, color = competitive_p)) +
    geom_point() +
    scale_x_discrete(labels = label_wrap(20)) +
    scale_size(range = c(5, 20)) +
    coord_flip() +
    theme_bw() +
    labs(y = bquote(R^2),
         size = "Number of \n SNPs",
         color = "Competitive \n p-value") +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 16))
  
  ggsave(filename_lasso, lasso_dotplot,
         width = 10, height = 9)
  
  return(list(top20 = top20_dotplot,
              lasso = lasso_dotplot))
}
#cad results is 7.4 GB so I deleted off my computer to save space
all_cad_dotplot <- create_dotplot("cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin.summary",
                                  "cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin_results.Rds",
                                  "cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin_top20_pathways.png",
                                  "cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin_lasso_pathways.png")
all_ascvd_dotplot <- create_dotplot("cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin.summary",
                                    "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin_results.Rds",
                                    "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin_top20_pathways.png",
                                    "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin_lasso_pathways.png")
all_unstable_dotplot <- create_dotplot("cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin.summary",
                                       "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin_results.Rds",
                                       "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin_top20_pathways.png",
                                       "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin_lasso_pathways.png")
all_occlusive_dotplot <- create_dotplot("cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.summary",
                                        "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin_results.Rds",
                                        "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin_top20_pathways.png",
                                        "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin_lasso_pathways.png")
all_ldl_dotplot <- create_dotplot("cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.summary",
                                  "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin_results.Rds",
                                  "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin_top20_pathways.png",
                                  "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin_lasso_pathways.png")
all_lpa_dotplot <- create_dotplot("cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.summary",
                                  "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin_results.Rds",
                                  "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin_top20_pathways.png",
                                  "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin_lasso_pathways.png")

#Customizable figure
prep_data <- function(prset_summary_path,
                      lasso_path) {
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
  
  prset_lasso_res <- readRDS(here(lasso_path))
  
  lasso_selected_pathways <- save_model_coef(prset_lasso_res$model) %>%
    filter(go_id != "(Intercept)") %>%
    pull(go_id)
  
  print(lasso_selected_pathways)
  
  # all pathways
  prset_summary %>%
    filter(competitive_p <= 0.05) %>%
    arrange(desc(prs_r2)) %>%
    slice_head(n = 10) %>%
    mutate(lasso = as.factor(case_when(
      set %in% lasso_selected_pathways ~ "Yes",
      TRUE ~ "No"
    ))) %>%
    select(set, term, prs_r2, full_r2, coefficient,
           standard_error, p, num_snp, competitive_p,
           lasso)
}

cad_pathways <- prep_data("cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin.summary",
                          "cad_subtype/result_data/cad_bin/c4d_prset_all_ukb_cad_bin_results.Rds") %>%
  mutate(subtype = "CAD")
# ascvd_pathways <- prep_data("cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin.summary",
#                             "cad_subtype/result_data/ascvd_subtype/c4d_prset_all_ukb_cad_ascvdsub_bin_results.Rds") %>%
#   mutate(subtype = "ASCVD")
stemi_pathways <- prep_data("cad_subtype/result_data/stemi_subtype/c4d_prset_all_ukb_cad_stemisub_bin.summary",
                            "cad_subtype/result_data/stemi_subtype/c4d_prset_all_ukb_cad_stemisub_bin_results.Rds") %>%
  mutate(subtype = "STEMI")
unstable_pathways <- prep_data("cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin.summary",
                               "cad_subtype/result_data/unstable_subtype/c4d_prset_all_ukb_cad_unstablesub_bin_results.Rds") %>%
  mutate(subtype = "Unstable")
occlusive_pathways <- prep_data("cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.summary",
                                "cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin_results.Rds") %>%
  mutate(subtype = "Occlusive")
ldl_pathways <- prep_data("cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.summary",
                          "cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin_results.Rds") %>%
  mutate(subtype = "LDL")
lpa_pathways <- prep_data("cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.summary",
                          "cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin_results.Rds") %>%
  mutate(subtype = "Lpa")

all_pathways <- ldl_pathways %>% 
  bind_rows(lpa_pathways) %>% 
  bind_rows(occlusive_pathways) %>% 
  bind_rows(unstable_pathways) %>% 
  bind_rows(stemi_pathways) %>% 
  bind_rows(cad_pathways)

all_pathways_plot <- ggplot(all_pathways,
                            aes(x = reorder_within(term, prs_r2, subtype), 
                                y = prs_r2, size = num_snp, color = competitive_p)) +
  geom_point() +
  # geom_label(aes(x = reorder_within(term, prs_r2, subtype), 
  #                y = prs_r2,
  #                label = ifelse(lasso == 1, "*", "")),
  #            inherit.aes = FALSE, nudge_x = 0, 
  #            nudge_y = 0, color = "red") +
  # gghighlight(lasso == 1) +
  # scale_x_discrete(labels = label_wrap(20)) +
  scale_size(range = c(5, 20)) +
  # scale_x_discrete(drop = TRUE) +
  scale_x_reordered() +
  scale_shape_manual(
    values = c("No" = 16, "Yes" = 8)
  ) +
  coord_flip() +
  theme_bw() +
  labs(y = bquote(R^2),
       size = "Number of \n SNPs",
       color = "Competitive \n p-value",
       shape = "Lasso-selected",
       label = NULL) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  facet_wrap(~ subtype, scales = "free",
             nrow = 3,
             ncol = 2) +
  guides(label = "none")
all_pathways_plot

ggsave("cad_subtype/result_data/c4d_pathways2.tiff",
       all_pathways_plot,
       width = 20, height = 9)

# Sensitivity analysis: Lpa gene removed
lpa_sensitivity_pathways <- prep_data("cad_subtype/result_data/lpa_subtype/c4d_no_lpa_prset_all_ukb_cad_lpasub_bin.summary",
                                      "cad_subtype/result_data/lpa_subtype/c4d_no_lpa_prset_all_ukb_cad_lpasub_bin_results.Rds") %>%
  mutate(subtype = "Lpa")

lpa_sens_plot <- ggplot(lpa_sensitivity_pathways,
       aes(x = reorder_within(term, prs_r2, subtype), 
           y = prs_r2, size = num_snp, color = competitive_p)) +
  geom_point() +
  scale_size(range = c(5, 20)) +
  scale_x_reordered() +
  scale_shape_manual(
    values = c("No" = 16, "Yes" = 8)
  ) +
  coord_flip() +
  theme_bw() +
  labs(y = bquote(R^2),
       size = "Number of \n SNPs",
       color = "Competitive \n p-value",
       shape = "Lasso-selected",
       label = NULL) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  guides(label = "none")
ggsave("cad_subtype/result_data/lpa_subtype/lpa_sensitivity_pathways.tiff",
       plot = lpa_sens_plot, 
       height = 9,
       width = 10)
