# Drug Target Connector

library(tidyverse)
library(biomaRt)
library(conflicted)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
library(data.table)
library(qs)
library(here)
library(networkD3)

# Load in data ----
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

drug_connector_raw <- fread("cad_subtype/wholedatabase_for_targetor.txt")

genes <- getBM(attributes= c("hgnc_symbol",
                             "ensembl_gene_id"),
               values = drug_connector_raw$gene,
               mart = ensembl)

drug_connector <- drug_connector_raw %>%
  left_join(genes, 
            by = c("gene" = "hgnc_symbol")) %>%
  separate_wider_delim(atc, 
                       names = c("atc_id",
                                 "drug",
                                 "cid"),
                       delim = "|") %>%
  mutate(drug = str_replace(drug, "NAME:", ""))

relevant_cad_drugs <- c("ISOSORBIDE_MONONITRATE",
                        "RIOCIGUAT",
                        "ISOSORBIDE_DINITRATE",
                        "APROTININ",
                        "AMBRISENTAN",
                        "MACITENTAN",
                        "ACETYLDIGITOXIN",
                        "EZETIMIBE",
                        "BELINOSTAT",
                        "OXPRENOLOL",
                        "BUPRANOLOL",
                        "THEOPHYLLINE",
                        "DESLANOSIDE",
                        "CARBIMAZOLE",
                        "PARACETAMOL",
                        "FENOFIBRATE",
                        "OXYTETRACYCLINE",
                        "FOLIC_ACID",
                        "VINCRISTINE",
                        "SITAXENTAN",
                        "DEMECOLCINE",
                        "SELENIUM_COMPOUNDS",
                        "PROCAINAMIDE",
                        "MITOBRONITOL",
                        "CHOLINE_ALFOSCERATE",
                        "PAZOPANIB",
                        "PHENOSULFONPHTHALEIN",
                        "DAUNORUBICIN",
                        "ALBENDAZOLE",
                        "DACTINOMYCIN",
                        "ALCURONIUM",
                        "BENAZEPRIL",
                        "PRAVASTATIN",
                        "METHYLROSANILINE",
                        "GREPAFLOXACIN",
                        "BUDIPINE",
                        "FLUORESCEIN",
                        "ROSUVASTATIN",
                        "BOSENTAN",
                        "LISINOPRIL",
                        "MOEXIPRIL",
                        "SPIRAPRIL",
                        "DELAPRIL",
                        "ZOFENOPRIL",
                        "QUINAPRIL",
                        "CILAZAPRIL",
                        "FOSINOPRIL",
                        "ATORVASTATIN",
                        "FLUVASTATIN", "LOVASTATIN",
                        "PITAVASTATIN", "SIMVASTATIN")

magma <- fread("cad_subtype/result_data/CAD.sumstat.genes.out") %>%
  janitor::clean_names()

# Lpa analysis ----
lpa_prset_summary <- fread(here("cad_subtype/result_data/lpa_subtype/c4d_prset_all_ukb_cad_lpasub_bin.summary")) %>%
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

lpa_prset_snps <- qread("cad_subtype/result_data/lpa_subtype/c4d_lpa_pathways_snps.qs")
lpa_prsice_snps <- fread("cad_subtype/result_data/lpa_subtype/c4d_prsice_all_ukb_cad_lpasub_bin.snp")

lpa_lasso_snps <- tibble(rsid = lpa_prset_snps[[1]] %>%
                           select(SNP,
                                  "GO:0001968",
                                  "GO:0034185",
                                  "GO:0004252") %>%
                           filter(if_any(c("GO:0001968",
                                           "GO:0034185",
                                           "GO:0004252"), ~.x != 0)) %>%
                           pull(SNP)) %>%
  mutate(dbsnp = "dbsnp") %>%
  relocate(dbsnp, .before = "rsid")

lpa_prsice_snps %>%
  filter(SNP %in% lpa_lasso_snps$rsid) %>%
  nrow()

write_delim(lpa_lasso_snps, file = "cad_subtype/result_data/lpa_lasso_snps.txt",
            delim = "\t", col_names = FALSE)

## after running it in snp nexus
lpa_lasso_genes <- fread("cad_subtype/result_data/lpa_subtype/near_gens_lpa lasso.txt") %>%
  janitor::clean_names()

lpa_genes_drugs <- lpa_lasso_genes %>%
  left_join(drug_connector, 
            by = c("overlapped_gene" = "gene")) %>%
  left_join(magma %>%
              select(gene, zstat, p),
            by = c("ensembl_gene_id" = "gene")) %>%
  filter(p < 0.05)
  # filter(drug %in% relevant_cad_drugs)

top_20_lpa_drugs <- lpa_genes_drugs %>%
  group_by(drug) %>%
  count() %>% 
  arrange(desc(n)) %>%
  ungroup() %>%
  slice_head(n = 20) %>%
  pull(drug)

snp_pathways_df <- lpa_prset_snps[[1]] %>%
  select(SNP,
         "GO:0001968",
         "GO:0034185",
         "GO:0004252") %>%
  filter(if_any(c("GO:0001968",
                  "GO:0034185",
                  "GO:0004252"), ~.x != 0)) %>% 
  pivot_longer(-SNP,
               names_to = "pathway",
               values_to = "contains_snp") %>%
  filter(contains_snp == 1) %>%
  group_by(SNP) %>%
  summarise(lasso_selected_pathway = paste(pathway, collapse = ", "))
  
lpa_genes_drugs_clean <- lpa_genes_drugs %>%
  select(variation_id, chromosome, position, overlapped_gene, ensembl_gene_id,
         type, drug, activity_type, zstat, p) %>%
  left_join(snp_pathways_df,
            by = c("variation_id" = "SNP")) %>%
  write_csv(., file = "cad_subtype/result_data/lpa_subtype/lpa_snps_genes_drugs.csv")

pathway1 <- lpa_prset_snps[[1]] %>%
  select(SNP,
         "GO:0001968",
         "GO:0034185",
         "GO:0004252") %>%
  filter(if_any(c("GO:0001968",
                  "GO:0034185",
                  "GO:0004252"), ~.x != 0)) %>%
  mutate(pathway1 = case_when(`GO:0001968` == 1 ~ "GO:0001968",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(pathway1)) %>%
  left_join(lpa_lasso_genes,
            by = c("SNP" = "variation_id")) %>%
  select(overlapped_gene, pathway1) %>%
  rename(source = pathway1,
         target = overlapped_gene) %>%
  distinct()

pathway2 <- lpa_prset_snps[[1]] %>%
  select(SNP,
         "GO:0001968",
         "GO:0034185",
         "GO:0004252") %>%
  filter(if_any(c("GO:0001968",
                  "GO:0034185",
                  "GO:0004252"), ~.x != 0)) %>%
  mutate(pathway2 = case_when(`GO:0034185` == 1 ~ "GO:0034185",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(pathway2)) %>%
  left_join(lpa_lasso_genes,
            by = c("SNP" = "variation_id")) %>%
  select(overlapped_gene, pathway2) %>%
  rename(source = pathway2,
         target = overlapped_gene) %>%
  distinct()

pathway3 <- lpa_prset_snps[[1]] %>%
  select(SNP,
         "GO:0001968",
         "GO:0034185",
         "GO:0004252") %>%
  filter(if_any(c("GO:0001968",
                  "GO:0034185",
                  "GO:0004252"), ~.x != 0)) %>%
  mutate(pathway3 = case_when(`GO:0004252` == 1 ~ "GO:0004252",
                              TRUE ~ NA_character_)) %>%
  filter(!is.na(pathway3)) %>%
  left_join(lpa_lasso_genes,
            by = c("SNP" = "variation_id")) %>%
  select(overlapped_gene, pathway3) %>%
  rename(source = pathway3,
         target = overlapped_gene) %>%
  distinct()

#need 1: source = pathway, target = gene
#need 2: source = gene, target = drug

# lpa_sankey_df <- pathway1 %>%
#   bind_rows(pathway2) %>%
#   bind_rows(pathway3) %>%
#   bind_rows(lpa_genes_drugs %>%
#               select(overlapped_gene, drug) %>%
#               filter(drug %in% top_20_drugs) %>%
#               rename(source = overlapped_gene,
#                      target = drug)) %>%
#   mutate(value = 1)

lpa_sankey_df <- lpa_genes_drugs %>%
  select(overlapped_gene, drug, zstat) %>%
  filter(drug %in% top_20_lpa_drugs) %>%
  rename(source = overlapped_gene,
         target = drug)

nodes <- data.frame(name = c(as.character(lpa_sankey_df$source), 
                             as.character(lpa_sankey_df$target)) %>% 
                      unique())

# remove nodes that don't have drugs associated or pathways associated
# basically nodes that don't link a gene to a drug
# anti_nodes <- lpa_sankey_df %>%
#   filter(id_target < 80,
#          id_source <= 3) %>%
#   pull(id_target)

lpa_sankey_df$id_source = match(lpa_sankey_df$source, nodes$name)-1 
lpa_sankey_df$id_target = match(lpa_sankey_df$target, nodes$name)-1

sankeyNetwork(Links = lpa_sankey_df, Nodes = nodes,
              Source = "id_source", Target = "id_target",
              Value = "zstat", NodeID = "name")

# LDL drug analysis ----
ldl_prset_summary <- fread(here("cad_subtype/result_data/ldl_subtype/c4d_prset_all_ukb_cad_ldlsub_bin.summary")) %>%
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

ldl_prset_snps <- qread("cad_subtype/result_data/ldl_subtype/c4d_ldl_pathways_snps.qs")
  
ldl_lasso_snps <- tibble(rsid = ldl_prset_snps[[1]] %>%
                           select(SNP,
                                  "GO:0042632",
                                  "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE",
                                  "GO:0034361",
                                  "GO:0034364",
                                  "REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE",
                                  "GO:0030301",
                                  "GO:0050750",
                                  "GO:0034185") %>%
                           filter(if_any(c("GO:0042632",
                                           "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE",
                                           "GO:0034361",
                                           "GO:0034364",
                                           "REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE",
                                           "GO:0030301",
                                           "GO:0050750",
                                           "GO:0034185"), ~.x != 0)) %>%
                           pull(SNP)) %>%
  mutate(dbsnp = "dbsnp") %>%
  relocate(dbsnp, .before = "rsid")

write_delim(ldl_lasso_snps, file = "cad_subtype/result_data/ldl_subtype/ldl_lasso_snps.txt",
            delim = "\t", col_names = FALSE)

## after running it in snp nexus
ldl_lasso_genes <- fread("cad_subtype/result_data/ldl_subtype/near_gens_ldl lasso.txt") %>%
  janitor::clean_names()

ldl_genes_drugs <- ldl_lasso_genes %>%
  left_join(drug_connector, 
            by = c("overlapped_gene" = "gene")) %>%
  left_join(magma %>%
              select(gene, zstat, p),
            by = c("ensembl_gene_id" = "gene")) %>%
  filter(p < 0.05)
# filter(drug %in% relevant_cad_drugs)

top_20_ldl_drugs <- ldl_genes_drugs %>%
  group_by(drug) %>%
  count() %>% 
  arrange(desc(n)) %>%
  ungroup() %>%
  # filter(!drug %in% c("EMETINE",
  #                     "TROGLITAZONE"))
  slice_head(n = 20) %>%
  pull(drug)

ldl_snp_pathways_df <- ldl_prset_snps[[1]] %>%
  select(SNP,
         "GO:0042632",
         "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE",
         "GO:0034361",
         "GO:0034364",
         "REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE",
         "GO:0030301",
         "GO:0050750",
         "GO:0034185") %>%
  filter(if_any(c("GO:0042632",
                  "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE",
                  "GO:0034361",
                  "GO:0034364",
                  "REACTOME_PLASMA_LIPOPROTEIN_CLEARANCE",
                  "GO:0030301",
                  "GO:0050750",
                  "GO:0034185"), ~.x != 0)) %>% 
  pivot_longer(-SNP,
               names_to = "pathway",
               values_to = "contains_snp") %>%
  filter(contains_snp == 1) %>%
  group_by(SNP) %>%
  summarise(lasso_selected_pathway = paste(pathway, collapse = ", "))

ldl_genes_drugs_clean <- ldl_genes_drugs %>%
  select(variation_id, chromosome, position, overlapped_gene, ensembl_gene_id,
         type, drug, activity_type, zstat, p) %>%
  left_join(ldl_snp_pathways_df,
            by = c("variation_id" = "SNP")) %>%
  write_csv(., file = "cad_subtype/result_data/ldl_subtype/ldl_snps_genes_drugs.csv")

ldl_sankey_df <- ldl_genes_drugs %>%
  select(overlapped_gene, drug, zstat) %>%
  filter(drug %in% top_20_ldl_drugs) %>%
  rename(source = overlapped_gene,
         target = drug)

ldl_nodes <- data.frame(name = c(as.character(ldl_sankey_df$source), 
                             as.character(ldl_sankey_df$target)) %>% 
                      unique())

ldl_sankey_df$id_source = match(ldl_sankey_df$source, ldl_nodes$name)-1 
ldl_sankey_df$id_target = match(ldl_sankey_df$target, ldl_nodes$name)-1

ldl_sankey <- sankeyNetwork(Links = ldl_sankey_df, Nodes = ldl_nodes,
              Source = "id_source", Target = "id_target",
              Value = "zstat", NodeID = "name")
saveNetwork(ldl_sankey, "cad_subtype/result_data/ldl_subtype/ldl_sankey.html")

# occlusive

occlusive_prset_summary <- fread(here("cad_subtype/result_data/occlusive_subtype/c4d_prset_all_ukb_cad_occlusivesub_bin.summary")) %>%
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

occlusive_prset_snps <- qread("cad_subtype/result_data/occlusive_subtype/c4d_occlusive_pathways_snps.qs")

occlusive_lasso_snps <- tibble(rsid = occlusive_prset_snps[[1]] %>%
                           select(SNP,
                                  "GO:0034185") %>%
                           filter(if_any(c("GO:0034185"), ~.x != 0)) %>%
                           pull(SNP)) %>%
  mutate(dbsnp = "dbsnp") %>%
  relocate(dbsnp, .before = "rsid")

write_delim(occlusive_lasso_snps, file = "cad_subtype/result_data/occlusive_subtype/occlusive_lasso_snps.txt",
            delim = "\t", col_names = FALSE)

## after running it in snp nexus
occlusive_lasso_genes <- fread("cad_subtype/result_data/occlusive_subtype/near_gens_occlusive lasso.txt") %>%
  janitor::clean_names()

occlusive_genes_drugs <- occlusive_lasso_genes %>%
  left_join(drug_connector, 
            by = c("overlapped_gene" = "gene")) %>%
  left_join(magma %>%
              select(gene, zstat, p),
            by = c("ensembl_gene_id" = "gene")) %>%
  filter(p < 0.05)
# filter(drug %in% relevant_cad_drugs)

top_20_occlusive_drugs <- occlusive_genes_drugs %>%
  group_by(drug) %>%
  count() %>% 
  arrange(desc(n)) %>%
  ungroup() %>%
  slice_head(n = 20) %>%
  pull(drug)

occlusive_snp_pathways_df <- occlusive_prset_snps[[1]] %>%
  select(SNP,
         "GO:0034185") %>%
  filter(if_any(c("GO:0034185"), ~.x != 0)) %>% 
  pivot_longer(-SNP,
               names_to = "pathway",
               values_to = "contains_snp") %>%
  filter(contains_snp == 1) %>%
  group_by(SNP) %>%
  summarise(lasso_selected_pathway = paste(pathway, collapse = ", "))

occlusive_genes_drugs_clean <- occlusive_genes_drugs %>%
  select(variation_id, chromosome, position, overlapped_gene, ensembl_gene_id,
         type, drug, activity_type, zstat, p) %>%
  left_join(occlusive_snp_pathways_df,
            by = c("variation_id" = "SNP")) %>%
  write_csv(., file = "cad_subtype/result_data/occlusive_subtype/occlusive_snps_genes_drugs.csv")

occlusive_sankey_df <- occlusive_genes_drugs %>%
  select(overlapped_gene, drug, zstat) %>%
  filter(drug %in% top_20_occlusive_drugs) %>%
  rename(source = overlapped_gene,
         target = drug)

occlusive_nodes <- data.frame(name = c(as.character(occlusive_sankey_df$source), 
                                 as.character(occlusive_sankey_df$target)) %>% 
                          unique())

occlusive_sankey_df$id_source = match(occlusive_sankey_df$source, occlusive_nodes$name)-1 
occlusive_sankey_df$id_target = match(occlusive_sankey_df$target, occlusive_nodes$name)-1

occlusive_sankey <- sankeyNetwork(Links = occlusive_sankey_df, Nodes = occlusive_nodes,
                            Source = "id_source", Target = "id_target",
                            Value = "zstat", NodeID = "name")
saveNetwork(occlusive_sankey, "cad_subtype/result_data/occlusive_subtype/occlusive_sankey.html")
