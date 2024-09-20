#Carriage duration
A <- 
  Duration_12 %>%
  
  ggplot(aes(x= serotype, y=dur)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size = 3, alpha=0.9) +
  theme_classic(base_size = 25) + 
  theme(plot.title = element_text(face = "bold", size = 25),
        axis.title.y = element_text(face = "bold", size = 25)) +
  scale_x_discrete(limits=c("3","4","6A","6B","9V", "14", "18C","19A", "19F", "23F"))  +
  labs(title="", x = "Serotypes", y = "Pneumococcal carriage duration (days)") +
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300))

ggsave(here("output", "Fig1A_carriage_prevalence.svg"),
       plot = (A), 
       width = 15, height = 10, unit="in", dpi = 300)

#median duration among serotype 3
Duration_12 %>%
  filter(serotype == 3) %>%
  summarise(Median = median(dur), SD = sd(dur),
            CI_L = Median - (SD * 1.96)/sqrt(50),
            CI_U = Median + (SD * 1.96)/sqrt(50))

#FigB line graph of pneumococcal carriage density overtime (log10CFU/ml)
f <- 
  pneumov_12 %>% 
  filter(!is.na(dens), pid == "PD069O") %>%
  mutate(dens = log10(dens),
         month = month(date, label = TRUE)) %>%
  group_by(day,month) %>%
  summarise(Median = median(dens), SD = sd(dens),
            CI_L = Median - (SD * 1.96)/sqrt(50),
            CI_U = Median + (SD * 1.96)/sqrt(50))

df <- bind_rows(
  f %>% mutate(new = "ST18362"))

B <- 
  ggplot() +
  geom_point(data = df, aes(x=factor(day), y = Median, color = new)) +
  geom_line(data = df, aes(x= factor(day), y = Median, group = new, color = new)) +
  ylim(0,8) +
  geom_errorbar(data = df, aes(x=factor(day), y = Median, color = new, ymin = CI_L, ymax = CI_U), width = 0.3, linewidth = 0.8) +
  theme_bw(base_size = 25, base_family = "Lato") +
  theme(legend.position = c(0.85, 0.878), legend.title = element_blank(),
        legend.box = "vertical", plot.title=element_text(size=25, face = "bold")) +
  scale_color_manual(values=c("ST18362" = "red")) + #, "all Serotype 3" = '#42d4f4')
  labs(title="",x ="Time (day)", y = "Pneumococcal carriage density \n (log10CFU/ml)")

ggsave(here("output", "Fig1B_carriage_density_dur.svg"),
       plot = (B), 
       width = 15, height = 10, unit="in", dpi = 300)
