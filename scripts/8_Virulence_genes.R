Sup_Fig_2 <- 
  Virulence_genes_PD069O_2 %>% 
  
  ggplot(aes(Day, `Virulence genes`, fill = count)) + 
  geom_tile() + 
  #theme_minimal() + 
  scale_fill_gradient(low="#0069ec", high="#ff2722") + 
  labs(title = "",x ="Time(days)", y ="Virulence genes") +
  scale_x_discrete(limits=c("2","14","21","87","119", "149", "205", "225", "253", "315", "337")) +
  theme_classic(base_size = 25)

ggsave(here("output", "Fig_Virulence_genes.svg"),
       plot = (Sup_Fig_2), 
       width = 15, height = 10, unit="in", dpi = 300)
