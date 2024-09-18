#Mutations from day 14 with day 2 as reference
unique_snp_day <- 
  snps_snpeff_meta_data %>%
  filter(day == 2 | day == 14 | day == 21 | day == 87 | day == 119 | day == 149 | 
           day == 205| day == 253 | day == 315 | day == 337) %>%
  select(day, POS, REF_ALT, TYPE, FTYPE, EFFECT, EFFECT_2, GENE, PRODUCT) %>%
  arrange(day) %>%
  distinct(POS, REF_ALT, .keep_all = TRUE) %>%
  filter(day != 2)

A <- 
  unique_snp_day %>%
  group_by(GENE) %>%
  tally() %>%
  filter(!is.na(GENE)) %>%
  
  ggplot(aes(x = GENE, y = n)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of nucleotide substitutions", x = "GENE") +
  theme_classic(base_size = 25) +
  #guides(fill=guide_legend("EFFECT")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(here("output", "Fig3B_snps.svg"),
       plot = (A), 
       width = 15, height = 10, unit="in", dpi = 300)

B <- 
  unique_snp_day %>%
  group_by(GENE, EFFECT_2) %>%
  tally() %>%
  filter(!is.na(GENE)) %>%
  
  ggplot(aes(x = GENE, y = n, fill = EFFECT_2)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  labs(title="", y="number of nucleotide substitutions", x = "GENE") +
  theme_classic(base_size = 25) +
  guides(fill=guide_legend("EFFECT")) + 
  #scale_y_break(c(8, 30), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.9, 0.878))

ggsave(here("output", "Fig4C_snps.svg"),
       plot = ((B| plot_layout(width = c(1)))), 
       width = 15, height = 10, unit="in", dpi = 300)