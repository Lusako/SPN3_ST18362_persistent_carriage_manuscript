#Load_packages
#BiocManager::install("ggtree")
pacman::p_load(char = c("ggpubr",'lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", 
                        "boot", "magrittr", "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs","sjstats",
                        "sjPlot", "oddsratio", "cowplot","lattice","ggcorrplot", "effects","glmmTMB","lme4",
                        "stringr", "broom", "ggsignif", "ggcorrplot", "rstatix", "stats", "gt", "webshot", "irr",
                        "writexl", "gt", "adegenet", "gplots", "ggnewscale", "ape","ggtree", "treeio",
                        "ggplotify", "tidytree", "svglite"))

#Import_data_for_plate_sweeps 
ENA_data <- import(here("data", "ENA_accession_numbers.tsv"))
Manifest_data_sw <- import(here("data", "new_manifest.xlsx")) %>%
  filter(`SAMPLE DESCRIPTION` == "PLATE SWEEP") %>%
  rename("SAMPLE_TYPE" = "SAMPLE TYPE") %>%
  mutate(SAMPLE_TYPE = case_when(SAMPLE_TYPE == "PNS" ~ "PNS",
                                 SAMPLE_TYPE == "RPNS" ~ "PNS"))

Manifest_metadata_sw <- import(here("data", "manifest_metadata.xlsx")) %>%
  filter(`SAMPLE DESCRIPTION` == "PLATE SWEEP", pid == "PD069O")  %>%
  rename("SAMPLE_TYPE" = "SAMPLE TYPE") %>%
  mutate(SAMPLE_TYPE = case_when(SAMPLE_TYPE == "PNS" ~ "PNS",
                                 SAMPLE_TYPE == "RPNS" ~ "PNS"))

ENA_Manifest_data_sw <- left_join(Manifest_data_sw,ENA_data)
sw_sequence_metadata <- left_join(Manifest_metadata_sw,ENA_Manifest_data_sw)

#Import_data_for_single_isolates
ENA_data <- import(here("data", "ENA_accession_numbers.tsv"))
Manifest_data <- import(here("data", "new_manifest.xlsx")) %>%
  filter(`SAMPLE DESCRIPTION` == "ISOLATE") %>%
  rename("SAMPLE_TYPE" = "SAMPLE TYPE") %>%
  mutate(SAMPLE_TYPE = case_when(SAMPLE_TYPE == "PNS" ~ "PNS",
                                 SAMPLE_TYPE == "RPNS" ~ "PNS"))

Manifest_metadata <- import(here("data", "manifest_metadata.xlsx")) %>%
  filter(`SAMPLE DESCRIPTION` == "ISOLATE", pid == "SPN3_PID") %>%
  select(-`CONC. (ng/ul)`, -`SAMPLE DESCRIPTION`) %>%
  rename("SAMPLE_TYPE" = "SAMPLE TYPE") %>%
  mutate(SAMPLE_TYPE = case_when(SAMPLE_TYPE == "PNS" ~ "PNS",
                                 SAMPLE_TYPE == "RPNS" ~ "PNS"))


ENA_Manifest_data <- left_join(Manifest_data,ENA_data)
sc_sequence_metadata <- left_join(Manifest_metadata,ENA_Manifest_data)

#Import_data_for_pneumodude_followup
pneumov_12 <- import(here("data", "pneumov_12.xlsx")) %>%
  filter(!is.na(date))

#Import_data_carriage_duration
Duration_12 <- import(here("data", "pneumov_12_duration.xlsx")) %>%
  filter(serogroup == "VT", !is.na(dur)) %>%
  select(serotype, dur)

#Import_data_serocalls
Serocall_6995 <- 
  import(here("../","Serocall_results.xlsx")) %>%
  rename("Lane name" = "Lane_id")

Serocall_PD069O_meta_data <- left_join(sw_sequence_metadata, Serocall_6995)

#Import_data_Virulence_genes
Virulence_genes_PD069O <- import(here("../Virulence_genes.xlsx")) %>% 
  select(-Lane_id)

Virulence_genes_PD069O_2 <- 
  Virulence_genes_PD069O %>%
  pivot_longer(!Day, names_to = "Virulence genes", values_to = "count")  %>% 
  mutate(Day = as.character(Day))
