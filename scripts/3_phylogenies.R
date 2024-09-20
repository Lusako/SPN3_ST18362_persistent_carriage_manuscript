#ST3_WG_tree Gubbins tree
#import data
ST3_WG_metadata <- import(here("../ST3_WG_metadata.csv")) %>%
  select(Lane_id, GPSC, ST, Country, Continent)

#import Newick files
filename_new4 <- "../ST3_WG_tree.nwk"
ST3_WG_tree <- ape::read.tree(filename_new4)
phytools::read.newick(filename_new4)

plot(ST3_WG_tree)

ST3_WG_tree = drop.tip(ST3_WG_tree, "Reference")

#Check the data and tree
head(ST3_WG_tree$tip.label)
head(ST3_WG_metadata$Lane_id)
ST3_WG_metadata$Lane_id %in% ST3_WG_tree$tip.label
ST3_WG_tree$tip.label %in% ST3_WG_metadata$Lane_id
ST3_WG_metadata$Lane_id[!ST3_WG_tree$tip.label %in% ST3_WG_metadata$Lane_id]

ST3_WG_tree <- groupClade(ST3_WG_tree, c(800, 1001, 891, 656, 819))

A <- 
  ST3_WG_tree %>%
  
  ggtree(layout='rectangular', aes(color=group)) %<+% 
  ST3_aln_tree_metadata +
  #geom_tiplab() +
  geom_tippoint(aes(color=group),
                size = 8) +
  ylim(1,650) +
  scale_color_manual(name = "Sequence type",
                     values = c('#42d4f4', '#3cb44b', '#ffe119', '#f58231','#4363d8','#e6194B'),
                     breaks = c("1","2","3","4","0","5"),
                     labels = c("Clade with ST700","Clade with ST180", "Clade with ST260",
                                "Clade with ST458", "Clade with Other ST", "ST18362"),
                     na.value = "white") +
  geom_hilight(node = 800, fill = "steelblue", extend = 300, alpha = .3)

GPSC <- data.frame("GPSC" = ST3_WG_metadata[,c("GPSC")])
rownames(GPSC) <- ST3_WG_metadata$Lane_id

Country <- data.frame("Country" = ST3_WG_metadata[,c("Country")])
rownames(Country) <- ST3_WG_metadata$Lane_id

Continent <- data.frame("Continent" = ST3_WG_metadata[,c("Continent")])
rownames(Continent) <- ST3_WG_metadata$Lane_id

h_E <-
  gheatmap(A, GPSC, offset= 2, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 20,
           legend_title="GPSC", font.size = 15) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 30)) +
  scale_fill_manual(name = "GPSC",
                    values = c('#f032e6', '#3cb44b', '#ffe119', '#f58231', '#4363d8', '#911eb4', '#42d4f4',
                               '#e6194B', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', "#e75480",
                               '#800000', '#aaffc3', '#808000', '#ffd8b1', "steelblue","black", "darkblue" ,
                               "brown","grey","#fffac8", "darkgreen", "darkgrey"),
                    breaks = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    labels = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    na.value = "white")

h_E <- h_E + new_scale_fill()

h_E2 <-
  gheatmap(h_E, Continent, offset= 1200, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 30,
           legend_title="Continent", font.size = 15) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 30)) + 
  theme(legend.title = element_text(size = 35), legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm")) +
  guides(color=guide_legend(override.aes = list(size=15)))

ggsave(here("output", "Fig2A_ST3_global.svg"),
       plot = ((h_E2| plot_layout(width = c(1)))), 
       width = 30, height = 30, unit="in", dpi = 300)

#GPSC10 new
#import data
GPSC10_Global_new <- import(here("../GPSC10_metadata_insilico_serotype.xlsx"))
ST3_aln_tree_metadata <- import(here("../ST3_aln_tree_metadata.csv")) %>%
  select(Lane_id, GPSC)

#import Newick files
filename_new <- "../GPSC10_Global_new.nwk"
GPSC10_tree_new <- ape::read.tree(filename_new)
phytools::read.newick(filename_new)

plot(GPSC10_tree_new)

#Check the data and tree
head(GPSC10_tree_new$tip.label)
head(GPSC10_Global_new$Lane_id)
GPSC10_Global_new$Lane_id %in% GPSC10_tree_new$tip.label
GPSC10_tree_new$tip.label %in% GPSC10_Global_new$Lane_id
GPSC10_Global_new$Lane_id[!GPSC10_tree_new$tip.label %in% GPSC10_Global_new$Lane_id]


tree_new <- tree_subset(GPSC10_tree_new, "45897_1_60", levels_back = 4)
tree_new[["tip.label"]]

ST3_tree_1 <- groupOTU(GPSC10_tree_new, .node=c("45897_1_58","45897_2_302","45897_1_62","45907_1_184","45897_1_50",
                                                "45897_1_30","45897_2_278","45897_1_54","45897_1_52","45897_2_286", "45897_1_60"))

B <- 
  ST3_tree_1 %>%
  
  ggtree(layout='circular', aes(color=group)) %<+% 
  GPSC10_Global_new +
  geom_tippoint(aes(color=group),
                size = 0.05) +
  xlim(-1000, NA) +
  ylim(NA, 1000) +
  scale_color_manual(name = "MLST",
                     values = c('#e6194B', 'darkgreen'),
                     breaks = c("1","0"),
                     labels = c("ST18362","Other ST"),
                     na.value = "white") #+
#geom_hilight(node = 1216, fill = "steelblue", extend = 300, alpha = .5) 
#geom_cladelabel(node=1280, label="Lokiarchaeota", color="green")

Serotype <- data.frame("serotype" = GPSC10_Global_new[,c("serotype")])
rownames(Serotype) <- GPSC10_Global_new$Lane_id

h_a <-
  gheatmap(B, Serotype, offset= 2, width=0.08, color = NA,
           colnames = FALSE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 30,
           legend_title="Serotype", font.size = 15) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 35)) +
  scale_fill_manual(name = "Serotype",
                    values = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#f032e6',
                               '#42d4f4', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
                               '#800000', '#aaffc3', '#808000', '#ffd8b1'),
                    breaks = c("3","6A","6C","7B","10A","11A","13","14","15B/15C",
                               "17F","19A","19F","22F","23A","23B","23F","24F"),
                    labels = c("3","6A","6C","7B","10A","11A","13","14","15B/15C",
                               "17F","19A","19F","22F","23A","23B","23F","24F"),
                    na.value = "white")

ggsave(here("output", "Fig2B_GPSC10_global_tree.svg"),
       plot = ((h_a| plot_layout(width = c(1)))), 
       width = 12, height = 8, unit="in", dpi = 300)

#ST3_aln_tree
#import data
ST3_aln_tree_metadata <- import(here("../ST3_aln_tree_metadata.csv")) %>%
  select(Lane_id, GPSC)

#import Newick files
filename_new2 <- "../ST3_aln_tree.nwk"
ST3_aln_tree <- ape::read.tree(filename_new2)
phytools::read.newick(filename_new2)

plot(ST3_aln_tree)

ST3_aln_tree = drop.tip(ST3_aln_tree, "ST3_cps_locus")

#Check the data and tree
head(ST3_aln_tree$tip.label)
head(ST3_aln_tree_metadata$Lane_id)
ST3_aln_tree_metadata$Lane_id %in% ST3_aln_tree$tip.label
ST3_aln_tree$tip.label %in% ST3_aln_tree_metadata$Lane_id
ST3_aln_tree_metadata$Lane_id[!ST3_aln_tree$tip.label %in% ST3_aln_tree_metadata$Lane_id]

tree_new <- tree_subset(ST3_aln_tree, "45897_1_58", levels_back = 4)
tree_new[["tip.label"]]

ST3_tree_2 <- groupOTU(ST3_aln_tree, .node=c("45897_1_58","45897_2_302","45897_1_62","45907_1_184","45897_1_50",
                                             "45897_1_30","45897_2_278","45897_1_54","45897_1_52","45897_2_286", "45897_1_60"))

ggtree(ST3_tree_2, aes(color=group)) + geom_tiplab()

#tree_10 <- groupClade(GPSC10_tree_new, c(1214))

C <- 
  ST3_tree_2 %>%
  
  ggtree(layout='rectangular', aes(color=group)) %<+% 
  ST3_aln_tree_metadata +
  #geom_tiplab() +
  geom_tippoint(aes(color=group),
                size = 10) +
  scale_color_manual(name = "MLST",
                     values = c('#e6194B', 'darkgreen'),
                     breaks = c("1","0"),
                     labels = c("ST18362","Other ST"),
                     na.value = "white")

GPSC <- data.frame("GPSC" = ST3_aln_tree_metadata[,c("GPSC")])
rownames(GPSC) <- ST3_aln_tree_metadata$Lane_id

h_C <-
  gheatmap(C, GPSC, offset= 2, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 30,
           legend_title="GPSC", font.size = 15) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 35)) +
  scale_fill_manual(name = "GPSC",
                    values = c('#f032e6', '#3cb44b', '#ffe119', '#f58231', '#4363d8', '#911eb4', '#42d4f4',
                               '#e6194B', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
                               '#800000', '#aaffc3', '#808000', '#ffd8b1', "steelblue","black", "darkblue" ,
                               "brown","grey","#e75480", "darkgreen", "darkgrey"),
                    breaks = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    labels = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    na.value = "white") +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30),
        legend.key.size = unit(2, "cm"))


ggsave(here("output", "Fig4A_cps_locus_ST3.svg"),
       plot = ((h_C| plot_layout(width = c(1)))), 
       width = 30, height = 30, unit="in", dpi = 300)

#PD069O_tree_AMR
#import data
PD069O_tree_metadata <- import(here("../PD069O_tree_metadata.csv")) %>%
  mutate(Day = as.character(Day), Benzylpenicillin = as.character(Benzylpenicillin),
         Ceftriaxone = as.character(Ceftriaxone))

#import Newick files
filename_new3 <- "../GPSC10_clean_full_final_tree.nwk"
PD069O_tree <- ape::read.tree(filename_new3)
phytools::read.newick(filename_new3)

plot(PD069O_tree)

PD069O_tree <- drop.tip(
  ladderize(
    root.phylo(PD069O_tree, outgroup = "Reference", resolve.root = TRUE),
    right = FALSE), tip = c("Reference"))

plot(PD069O_tree)

#Check the data and tree
head(PD069O_tree$tip.label)
head(PD069O_tree_metadata$Lane_id)
PD069O_tree_metadata$Lane_id %in% PD069O_tree$tip.label
PD069O_tree$tip.label %in% PD069O_tree_metadata$Lane_id
PD069O_tree_metadata$Lane_id[!PD069O_tree$tip.label %in% PD069O_tree_metadata$Lane_id]


D <- 
  PD069O_tree %>%
  
  ggtree(layout='rectangular') %<+% 
  PD069O_tree_metadata +
  #geom_tiplab() +
  geom_tippoint(aes(color=Day),
                size = 10) +
  scale_color_manual(name = "Days",
                     values = c('#f032e6', '#3cb44b', '#ffe119', '#f58231', '#4363d8', '#911eb4', '#42d4f4',
                                '#e6194B', '#bfef45', '#fabed4', '#469990'),
                     breaks = c("2","14","21","87","119","149","205","225","253",
                                "315","337"),
                     labels = c("2","14","21","87","119","149","205","225","253",
                                "315","337"),
                     na.value = "white")

Tetracycline <- data.frame("Tetracycline" = PD069O_tree_metadata[,c("Tetracycline")])
rownames(Tetracycline) <- PD069O_tree_metadata$Lane_id

`Co-trimoxazole` <- data.frame("Co-trimoxazole" = PD069O_tree_metadata[,c("Co-trimoxazole")])
rownames(`Co-trimoxazole`) <- PD069O_tree_metadata$Lane_id

Oxacillin <- data.frame("Oxacillin" = PD069O_tree_metadata[,c("Oxacillin")])
rownames(Oxacillin) <- PD069O_tree_metadata$Lane_id

Erythromycin <- data.frame("Erythromycin" = PD069O_tree_metadata[,c("Erythromycin")])
rownames(Erythromycin) <- PD069O_tree_metadata$Lane_id

Benzylpenicillin <- data.frame("Benzylpenicillin" = PD069O_tree_metadata[,c("Benzylpenicillin")])
rownames(Benzylpenicillin) <- PD069O_tree_metadata$Lane_id

Ceftriaxone <- data.frame("Ceftriaxone" = PD069O_tree_metadata[,c("Ceftriaxone")])
rownames(Ceftriaxone) <- PD069O_tree_metadata$Lane_id

h_D <-
  gheatmap(D, Tetracycline, offset= 0.1, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.5,
           legend_title="Tetracycline", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Tetracycline",
                    values = c('#ff2722'),
                    breaks = c("R"),
                    labels = c("Resistant"),
                    na.value = "white")

h_D <- h_D + new_scale_fill()

h_D2 <-
  gheatmap(h_D, `Co-trimoxazole`, offset= 0.65, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.75,
           legend_title="Co-trimoxazole", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Co-trimoxazole",
                    values = c('#ff2722'),
                    breaks = c("R"),
                    labels = c("Resistant"),
                    na.value = "white")

h_D2 <- h_D2 + new_scale_fill()

h_D3 <-
  gheatmap(h_D2, Oxacillin, offset= 1.2, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.2,
           legend_title="Oxacillin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Oxacillin",
                    values = c('#ff2722'),
                    breaks = c("R"),
                    labels = c("Resistant"),
                    na.value = "white")

h_D3 <- h_D3 + new_scale_fill()

h_D4 <-
  gheatmap(h_D3, Erythromycin, offset= 1.75, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.65,
           legend_title="Erythromycin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Erythromycin",
                    values = c('#0069ec'),
                    breaks = c("S"),
                    labels = c("Susceptible"),
                    na.value = "white")

h_D4 <- h_D4 + new_scale_fill()

h_D5 <-
  gheatmap(h_D4, Benzylpenicillin, offset= 2.3, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.8,
           legend_title="Benzylpenicillin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Benzylpenicillin",
                    values = c('#ff7470ff', '#ff7470ff'),
                    breaks = c("0.5","0.75"),
                    labels = c("0.5", "0.5"),
                    na.value = "white")

h_D5 <- h_D5 + new_scale_fill()

h_D6 <-
  gheatmap(h_D5, Ceftriaxone, offset= 2.85, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 0.5,
           legend_title="Ceftriaxone", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 2)) +
  scale_fill_manual(name = "Ceftriaxone",
                    values = c('#0069ec'),
                    breaks = c("0.19"),
                    labels = c("0.19"),
                    na.value = "white") + 
  theme(legend.title = element_text(size = 25), legend.text = element_text(size = 25),
        legend.key.size = unit(1, "cm"))


ggsave(here("output", "Sup_Fig1A_PD069O_tree.svg"),
       plot = ((h_D6| plot_layout(width = c(1)))), 
       width = 20, height = 15, unit="in", dpi = 300)

#GPS serotype 3
#import data
ST3_Whole_genome_metadata <- import(here("../ST3_whole_genome_Metadata.csv")) %>%
  mutate(penMic = as.character(penMic))

#import Newick files
filename_new2 <- "../ST3_whole_genome.nwk"
ST3_Whole_genome_tree <- ape::read.tree(filename_new2)
phytools::read.newick(filename_new2)

plot(ST3_Whole_genome_tree)

ST3_Whole_genome_tree = drop.tip(ST3_Whole_genome_tree, "GCA_000210955_1_ASM21095v1_genomic")

#Check the data and tree
head(ST3_Whole_genome_tree$tip.label)
head(ST3_Whole_genome_metadata$Lane_id)
ST3_Whole_genome_metadata$Lane_id %in% ST3_Whole_genome_tree$tip.label
ST3_Whole_genome_tree$tip.label %in% ST3_Whole_genome_metadata$Lane_id
ST3_Whole_genome_metadata$Lane_id[!ST3_Whole_genome_tree$tip.label %in% ST3_Whole_genome_metadata$Lane_id]

tree_new <- tree_subset(ST3_Whole_genome_tree, "45897_1_58", levels_back = 4)
tree_new[["tip.label"]]

ST3_tree_3 <- groupOTU(ST3_Whole_genome_tree, .node=c("45897_1_58","45897_2_302","45897_1_62","45907_1_184","45897_1_50",
                                                      "45897_1_30","45897_2_278","45897_1_54","45897_1_52","45897_2_286", "45897_1_60"))

ggtree(ST3_tree_3, aes(color=group)) + geom_tiplab()

#tree_10 <- groupClade(GPSC10_tree_new, c(1214))

E <- 
  ST3_tree_3 %>%
  
  ggtree(layout='rectangular', aes(color=group)) %<+% 
  ST3_Whole_genome_metadata +
  #geom_tiplab() +
  geom_tippoint(aes(color=group),
                size = 5) +
  scale_color_manual(name = "MLST",
                     values = c('#e6194B', 'darkgreen'),
                     breaks = c("1","0"),
                     labels = c("ST18362","Other ST"),
                     na.value = "white")

GPSC <- data.frame("GPSC" = ST3_Whole_genome_metadata[,c("GPSC")])
rownames(GPSC) <- ST3_Whole_genome_metadata$Lane_id

Tetracycline <- data.frame("Tetracycline" = ST3_Whole_genome_metadata[,c("Tetracycline")])
rownames(Tetracycline) <- ST3_Whole_genome_metadata$Lane_id

`Co-trimoxazole` <- data.frame("Co_Trimoxazole" = ST3_Whole_genome_metadata[,c("Co_Trimoxazole")])
rownames(`Co-trimoxazole`) <- ST3_Whole_genome_metadata$Lane_id

Clindamycin <- data.frame("Clindamycin" = ST3_Whole_genome_metadata[,c("Clindamycin")])
rownames(Clindamycin) <- ST3_Whole_genome_metadata$Lane_id

Erythromycin <- data.frame("Erythromycin" = ST3_Whole_genome_metadata[,c("Erythromycin")])
rownames(Erythromycin) <- ST3_Whole_genome_metadata$Lane_id

Penicillin <- data.frame("penMic" = ST3_Whole_genome_metadata[,c("penMic")])
rownames(Penicillin) <- ST3_Whole_genome_metadata$Lane_id

Chloramphenicol <- data.frame("Chloramphenicol" = ST3_Whole_genome_metadata[,c("Chloramphenicol")])
rownames(Chloramphenicol) <- ST3_Whole_genome_metadata$Lane_id

Fluoroquinolones <- data.frame("Fluoroquinolones" = ST3_Whole_genome_metadata[,c("Fluoroquinolones")])
rownames(Fluoroquinolones) <- ST3_Whole_genome_metadata$Lane_id

h_E <-
  gheatmap(E, GPSC, offset= 0.1, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 50,
           legend_title="GPSC", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "GPSC",
                    values = c('#f032e6', '#3cb44b', '#ffe119', '#f58231', '#4363d8', '#911eb4', '#42d4f4',
                               '#e6194B', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
                               '#800000', '#aaffc3', '#808000', '#ffd8b1', "steelblue","black", "darkblue" ,
                               "brown","grey","#e75480", "darkgreen", "darkgrey"),
                    breaks = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    labels = c("1","3","5","9","10","11","12","43","51",
                               "56","83","92","103","234","253","309", "334","363","371","646","666",
                               "790","791", "Not assigned"),
                    na.value = "white") 

h_E <- h_E + new_scale_fill()

h_E2 <-
  gheatmap(h_E, Tetracycline, offset= 1200, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Tetracycline", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Tetracycline",
                    values = c('#ff2722', "#0069ec"),
                    breaks = c("R", "S"),
                    labels = c("Resistant", "Susceptible"),
                    na.value = "white")

h_E2 <- h_E2 + new_scale_fill()

h_E3 <-
  gheatmap(h_E2, `Co-trimoxazole`, offset= 2400, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Co-trimoxazole", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Co-trimoxazole",
                    values = c('#ff2722', "#0069ec", "yellow"),
                    breaks = c("R", "S", "I"),
                    labels = c("Resistant", "Susceptible", "Intermediate"),
                    na.value = "white")

h_E3 <- h_E3 + new_scale_fill()

h_E4 <-
  gheatmap(h_E3, Clindamycin, offset= 3600, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Clindamycin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Clindamycin",
                    values = c('#ff2722', "#0069ec"),
                    breaks = c("R", "S"),
                    labels = c("Resistant", "Susceptible"),
                    na.value = "white")

h_E4 <- h_E4 + new_scale_fill()

h_E5 <-
  gheatmap(h_E4, Erythromycin, offset= 4800, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Erythromycin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Erythromycin",
                    values = c('#ff2722', "#0069ec"),
                    breaks = c("R", "S"),
                    labels = c("Resistant", "Susceptible"),
                    na.value = "white")

h_E5 <- h_E5 + new_scale_fill()

h_E6 <-
  gheatmap(h_E5, Penicillin, offset= 6000, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Penicillin", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Penicillin",
                    values = c('#0069ec', '#0069ec', '#0069ec', '#ff7470ff','#ff7470ff',
                               '#ff7470ff', '#ff2722'),
                    breaks = c("0","0.03","0.06","0.12","0.25","0.5", "2"),
                    labels = c("0","0.03","0.06","0.12","0.25","0.5", "2"),
                    na.value = "white")

h_E6 <- h_E6 + new_scale_fill()

h_E7 <-
  gheatmap(h_E6, Chloramphenicol, offset= 7200, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 90,
           legend_title="Chloramphenicol", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Chloramphenicol",
                    values = c('#ff2722', "#0069ec"),
                    breaks = c("R", "S"),
                    labels = c("Resistant", "Susceptible"),
                    na.value = "white")

h_E7 <- h_E7 + new_scale_fill()

h_E8 <-
  gheatmap(h_E7, Fluoroquinolones, offset= 8400, width=0.08, color = NA,
           colnames = TRUE, colnames_position = "top", 
           colnames_angle = 90, colnames_offset_y = 80,
           legend_title="Fluoroquinolones", font.size = 8) +
  scale_x_ggtree() +
  scale_y_continuous(expand=c(0, 70)) +
  scale_fill_manual(name = "Fluoroquinolones",
                    values = c('#ff2722', "#0069ec"),
                    breaks = c("R", "S"),
                    labels = c("Resistant", "Susceptible"),
                    na.value = "white") + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"))


ggsave(here("output", "Sup_Fig1B_PD069O_tree.svg"),
       plot = ((h_E8| plot_layout(width = c(1)))), 
       width = 20, height = 15, unit="in", dpi = 300)

