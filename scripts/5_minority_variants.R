#minority variants min cov 20
#exclude #45907_2#174 and "45897_2#160" due to multiple carriage and "45907_2#141" due to failed QC
snps_snpeff_lofreq <- 
  import(here("../../lofreq_GPSC10_ref/", "lofreq_final.csv")) %>%
  rename("gene_name" = "ANN[*].GENE",
         "effect" = "ANN[*].EFFECT") %>%
  mutate(Lane_id = str_extract(`45897_1#211 CHROM`, "^[^\\s]+"),
         CHROM = str_extract(`45897_1#211 CHROM`, "\\s.+")) %>% 
  select(-`45897_1#211 CHROM`) %>% 
  filter(AF <= 0.5 & AF >= 0.15, Lane_id != "45907_2#174") #, !grepl("GL",gene_name)

#write_xlsx(snps_snpeff_lofreq,"/Users/lusako/Downloads/lofreq_final.xlsx")
coverage <- 
  import(here("../../lofreq_QC_and_scripts/", "coverage_report.xls")) %>%
  select(Lane_id, genome_coverage, minAF)

Kraken <- 
  import(here("../../lofreq_QC_and_scripts/", "qc_summary.csv")) %>%
  rename("Lane_id" = "Species") %>%
  select(-Total) %>%
  pivot_longer(!Lane_id, names_to = "specie", values_to = "percentage") %>%
  filter(percentage >= 0.1, Lane_id != "45907_2#174")

Kraken_spn <- 
  Kraken %>%
  filter(specie == "Streptococcus pneumoniae")

meta_data <- 
  sw_sequence_metadata %>%
  rename("Lane_id" = "Lane name") %>%
  select(pid,Lane_id,day,dens) %>%
  filter(Lane_id != "45907_2#174", Lane_id != "45907_2#141", Lane_id != "45897_2#160")

list_df = list(meta_data,snps_snpeff_lofreq,coverage,Kraken_spn) 
snps_snpeff_lofreq_meta_data <- list_df %>% reduce(left_join, by='Lane_id')

######
Kraken_meta_data <- 
  left_join(meta_data, Kraken)

A <- 
  Kraken_meta_data %>%
  group_by(pid,specie,day,percentage) %>%
  tally() %>%
  mutate(day = as.factor(day)) %>%
  
  ggplot(aes(x = day, y = percentage, fill = specie)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"),width=0.9) +
  ylim(0,100) +
  scale_x_discrete(limits=c("0","2","28","87","119","149","178","205","253","289","315","337")) +
  labs(title="", x = "Day", y = "Percentage of reads from whole plate sweeps sequencing") +
  scale_fill_discrete(limits = c("Streptococcus pneumoniae","Staphylococcus aureus", 
                                 "Enterococcus casseliflavus","Enterococcus casseliflavus",
                                 "Unclassified","Staphylococcus epidermidis","Streptococcus mitis",
                                 "Streptococcus pseudopneumoniae","Haemophilus influenzae",
                                 "Streptococcus oralis", "Streptococcus sp. I-G2", "Homo sapiens")) +
  theme_classic(base_size = 15)

ggsave(here("output", "Sup_Fig3_species.svg"),
       plot = (A), 
       width = 15, height = 10, unit="in", dpi = 300)

#Genic
in_gene <-
  snps_snpeff_lofreq_meta_data %>%
  filter(effect != "downstream_gene_variant", effect != "upstream_gene_variant") %>% 
  unite("REF_ALT", REF:ALT, remove = FALSE) %>% 
  mutate(effect_2 = str_extract_all(effect, "[a-z]+", simplify = TRUE)[,1])

#unique variants from day 2 with day 0 as reference
unique_snp_day <- 
  in_gene %>%
  filter(day == 0 | day == 2 | day == 28 | day == 87 | day == 119 | day == 149 | day == 178 | 
           day == 205| day == 253 | day == 289 | day == 315 | day == 337) %>%
  select(day, POS, REF_ALT, effect, gene_name) %>%
  arrange(day) %>%
  distinct(POS, REF_ALT, .keep_all = TRUE) %>%
  filter(day != 0)

#iSNVs
unique_snp_day_meta_data <-
  left_join(unique_snp_day, in_gene) %>%
  filter(!is.na(gene_name), gene_name != "GL183_00005", day != 2, day != 149, day != 337, 
         effect != "conservative_inframe_deletion",effect != "disruptive_inframe_deletion", 
         effect != "frameshift_variant", effect != "frameshift_variant&stop_gained", 
         effect != "conservative_inframe_insertion", effect != "disruptive_inframe_insertion")
#GL183_00005 is CP046355 whole genome

B <- 
  unique_snp_day_meta_data %>%
  group_by(gene_name) %>%
  tally() %>%
  
  ggplot(aes(x = reorder(gene_name,-n), y = n)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of intra-host Single Nucleotide Variants", x = "GENE") +
  theme_classic(base_size = 25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(here("output", "Fig 6A_isnvs.svg"),
       plot = (B), 
       width = 15, height = 10, unit="in", dpi = 300)

C <- 
  unique_snp_day_meta_data %>%
  group_by(gene_name, effect) %>%
  tally() %>%
  
  ggplot(aes(x = reorder(gene_name,-n), y = n, fill = effect)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of intra-host Single Nucleotide Variants", x = "GENE") +
  theme_classic(base_size = 25) +
  guides(fill=guide_legend("EFFECT")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.7, 0.878))

ggsave(here("output", "Fig 6B_isnvs.svg"),
       plot = (C), 
       width = 15, height = 10, unit="in", dpi = 300)

#Substitution Rate (per 1000 nucleotide)
num_days <- 
  unique_snp_day_meta_data %>%
  group_by(day, gene_name) %>%
  tally() %>%
  ungroup() %>%
  group_by(gene_name) %>%
  tally() %>%
  ungroup() %>%
  rename("num_days" = "n") %>%
  mutate(num_days = as.numeric(num_days))

num_nuc_sub <- 
  unique_snp_day_meta_data %>%
  group_by(POS, REF_ALT, gene_name) %>%
  tally() %>%
  ungroup() %>%
  group_by(gene_name)%>%
  tally() %>%
  ungroup() %>%
  rename("num_nuc_sub" = "n") %>%
  mutate(num_nuc_sub = as.numeric(num_nuc_sub))

gene_size_product <- 
  import(here("../../", "CP046355_gene_sizes_and_cds_products_with_notes.txt")) %>%
  rename("gene_name"= "Gene/LocusTag",
         "gene_size" = "Size", #number of nucleoitides
         "CDS_product" = "CDS Product") %>% 
  filter(gene_name != "GL183_00005", gene_name != "tnpB") #GL183_00005 is CP046355 whole genome
#summarise tnpB

gene_size_product_in_gene  <- 
  rbind(gene_size_product, list('tnpB', 338, 
                                "IS66 family insertion sequence element accessory protein TnpB",
                                "Derived by automated computational analysis using gene prediction method: Protein Homology."))

list_df = list(num_days, num_nuc_sub, gene_size_product_in_gene) 
gene_size_num_days <- list_df %>% reduce(left_join, by='gene_name') %>%
  mutate(sub_day_1000 = (num_nuc_sub/gene_size)*1000/num_days)

#write_xlsx(gene_size_num_days,"/Users/lusako/Downloads/gene_size_num_days.xlsx")

color_palette <- c(
  "1" = "#377EB8",  # Blue
  "2" = "#F781BF",  # Pink 
  "3" = "#4DAF4A",  # Green
  "4" = "#FF7F00",  # Orange
  "5" = "#FFFF33",  # Yellow 
  "6" = "#A65628",  # Brown
  "7" = "#E41A1C"   # Red
)

D <- 
  gene_size_num_days %>%
  mutate(num_days = as.character(num_days)) %>%
  
  ggplot(aes(x = reorder(gene_name,-sub_day_1000), y = sub_day_1000, fill = num_days)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) +  # Bar plot with dodged bars for different samples
  scale_fill_manual(values = color_palette) +  # Apply the color palette
  labs(#title = "Substitution Rate Per Gene and Sample",
    x = "Gene",
    y = "Substitution Rate (per 1000 nucleotide)") +
  #facet_wrap(~num_days) +
  theme_minimal() +  # Minimal theme for a clean look
  theme_classic(base_size = 25) +
  guides(fill=guide_legend("Sampled times")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.8, 0.7))

ggsave(here("output", "Fig 6C_iSNVs.svg"),
       plot = (D), 
       width = 15, height = 10, unit="in", dpi = 300)

#Intergenic
intergenic <-
  snps_snpeff_lofreq_meta_data %>%
  filter(effect == "downstream_gene_variant" | effect == "upstream_gene_variant") %>% 
  unite("REF_ALT", REF:ALT, remove = FALSE) %>% 
  mutate(effect_2 = str_extract_all(effect, "[a-z]+", simplify = TRUE)[,1])

#unique variants from day 2 with day 0 as reference
unique_intergenic_snp_day <- 
  intergenic %>%
  filter(day == 0 | day == 2 | day == 28 | day == 87 | day == 119 | day == 149 | day == 178 | 
           day == 205| day == 253 | day == 289 | day == 315 | day == 337) %>%
  select(day, POS, REF_ALT, effect, gene_name) %>%
  arrange(day) %>%
  distinct(POS, REF_ALT, .keep_all = TRUE) %>%
  filter(day != 0)

unique_intergenic_snp_day_meta_data <-
  left_join(unique_intergenic_snp_day, intergenic) %>%
  filter(!is.na(gene_name), gene_name != "GL183_00005", day != 2, day != 149, day != 337,
         REF_ALT %in% c("C_T", "G_A", "T_C", "A_G", "G_T", "T_G", "T_A", "A_C", "A_T",
                        "C_A", "C_G", "G_C"))
#GL183_00005 is CP046355 whole genome

E <- 
  unique_intergenic_snp_day_meta_data %>%
  group_by(gene_name) %>%
  tally() %>%
  
  ggplot(aes(x = reorder(gene_name, -n), y = n)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of intra-host Single Nucleotide Variants", x = "GENE") +
  theme_classic(base_size = 25) +
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(here("output", "Fig 6D_intergenic_iSNVs.svg"),
       plot = (E), 
       width = 15, height = 10, unit="in", dpi = 300)

F <- 
  unique_intergenic_snp_day_meta_data %>%
  group_by(gene_name, effect) %>%
  tally() %>%
  
  ggplot(aes(x = gene_name, y = n, fill = effect)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of intra-host Single Nucleotide Variants", x = "GENE") +
  theme_classic(base_size = 25) +
  guides(fill=guide_legend("EFFECT")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.7, 0.878))

ggsave(here("output", "Fig 6E_intergenic_iSNVs.svg"),
       plot = (F), 
       width = 15, height = 10, unit="in", dpi = 300)

#Substitution Rate (per 1000 nucleotide)
num_days_intergenic <- 
  unique_intergenic_snp_day_meta_data %>%
  group_by(day, gene_name) %>%
  tally() %>%
  ungroup() %>%
  group_by(gene_name) %>%
  tally() %>%
  ungroup() %>%
  rename("num_days" = "n") %>%
  mutate(num_days = as.numeric(num_days))

num_nuc_sub_intergenic <- 
  unique_intergenic_snp_day_meta_data %>%
  group_by(POS, REF_ALT, gene_name) %>%
  tally() %>%
  ungroup() %>%
  group_by(gene_name)%>%
  tally() %>%
  ungroup() %>%
  rename("num_nuc_sub" = "n") %>%
  mutate(num_nuc_sub = as.numeric(num_nuc_sub))

gene_size <- 
  import(here("../../snpEff/", "CP046355_gene_sizes_and_cds_products_with_notes.txt")) %>%
  rename("gene_name"= "Gene/LocusTag",
         "gene_size" = "Size", #number of nucleoitides
         "CDS_product" = "CDS Product") %>% 
  filter(gene_name != "GL183_00005", gene_name != "pstA", gene_name != "galE", gene_name != "pstC") #GL183_00005 is CP046355 whole genome
#summarise pstA, galE and pstC

#creating a dataframe

new_genes <- data.frame(
  gene_name = c('pstA', 'galE', 'pstC'),
  gene_size = c(851, 861, 891),
  CDS_product = c("phosphate ABC transporter permease PstA","UDP-glucose 4-epimerase GalE", 
                  "phosphate ABC transporter permease subunit PstC"),
  Note = c("Derived by automated computational analysis using gene prediction method: Protein Homology."))

df2 <- as.data.frame(gene_size)

gene_size_intergenic <- rbind(df2, new_genes)

list_df = list(num_days_intergenic, num_nuc_sub_intergenic, gene_size_intergenic) 
gene_size_num_days_intergenic <- list_df %>% reduce(left_join, by='gene_name') %>%
  mutate(sub_day_1000 = (num_nuc_sub/gene_size)*1000/num_days)

color_palette_2 <- c(
  "1" = "#377EB8",  # Blue
  "2" = "#F781BF",  # Pink 
  "3" = "#4DAF4A",  # Green
  "4" = "#FF7F00",  # Orange
  "5" = "#FFFF33",  # Yellow 
  "6" = "#A65628",  # Brown
  "7" = "#E41A1C",   # Red
  "8" = "#000000"   # Black
)

G <- 
  gene_size_num_days_intergenic %>%
  mutate(num_days = as.character(num_days)) %>%
  
  ggplot(aes(x = reorder(gene_name,-sub_day_1000), y = sub_day_1000, fill = num_days)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), width = 0.7) +  # Bar plot with dodged bars for different samples
  scale_fill_manual(values = color_palette) +  # Apply the color palette
  labs(#title = "Substitution Rate Per Gene and Sample",
    x = "Gene",
    y = "Substitution Rate (per 1000 nucleotide)") +
  #facet_wrap(~num_days) +
  theme_minimal() +  # Minimal theme for a clean look
  theme_classic(base_size = 25) +
  guides(fill=guide_legend("Sampled times")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.8, 0.7))

ggsave(here("output", "Fig 6F_intergenic_iSNVs.svg"),
       plot = (G), 
       width = 15, height = 10, unit="in", dpi = 300)

#Supplementary table 2
gene_G <- 
  gene_size_num_days_intergenic %>%
  mutate(gene_region = "intergenic")

gene_G_2 <- gene_size_num_days %>%
  mutate(gene_region = "genic")

intra_host_single_nucleotide_variants <- rbind(gene_G_2, gene_G)

#Supplementary table 2
write_xlsx(intra_host_single_nucleotide_variants,"/Users/lusako/Downloads/Supplementary_table_2.xlsx")
