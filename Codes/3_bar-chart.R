
# Load the .rds object 

CIMDEC_run_rarF <- readRDS("R_objects/CIMDEC_run_rarF.rds")

library("phyloseq")

# fonctions pour la selection des couleurs pour les graphiques 
get_cols <- function (n){
  col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
           "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
           "#ccebc5", "#ffed6f")
  
  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
            "#a6cee3", "#006d2c", "#4d4d4d", "#d73027",
            "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
            "#ffff99", "#b15928")
  colorRampPalette(col2)(n)
}


######### Summary stats ########

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) # devtools::install_github("microsud/microbiomeutilities")

summarize_phyloseq(CIMDEC_run_raw)
summarize_phyloseq(CIMDEC_run_rarF)

# Histogramme pour montrer la distribution des lectures par echantillons
x<-plot_bar(CIMDEC_run_rarF)

read_depth <- x$data %>% 
  ggplot(aes(x=Abundance, y=Ech))+
  geom_col()+
  facet_wrap(vars(Type), nrow = 1, scales = "free_y");read_depth

# Comme le nb de lecture dans le temoin est tres bas (negligeable) le temoin a ete enleve du graphique final
CIMDEC_run_rarF_c <- subset_samples(CIMDEC_run_rarF, Type != "3_Témoin")

################################################
CIMDEC_run_rarF_dat <- psmelt(CIMDEC_run_rarF_c)

str(CIMDEC_run_rarF_dat)

predefined_species = c("Botrytis cinerea complex",
                       "Elsinoe ampelina",
                       "Phyllosticta ampelicida",
                       "Diaporthe ampelina complex",
                       "Sclerotinia sclerotiorum complex",
                       "Alternaria dauci",
                       "Colletotrichum acutatum",
                       "Plasmopara pusilla",
                       "Neocercospora carotae complex",
                       "Greeneria uvicola",
                       "Coniella diplodiella")

#################################### 
library("ggtext")
library("dplyr")

CIMDEC_run_rarF_dat <- CIMDEC_run_rarF_dat %>%
  mutate(Species = case_when(
    Species %in% 
      predefined_species ~ Species,  # Keep valid species as is
    TRUE ~ "Other"  # Replace other species with "other"
  ));CIMDEC_run_rarF_dat

#################################### 

library("tidyr")
otu_rel_abund <- CIMDEC_run_rarF_dat %>%
  group_by(Ech, Type) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon") %>% 
  mutate(Ech = factor(Ech, 
                              levels=c("A.dauci", 
                                       "C.carotae", 
                                       "S.sclerotiorum", 
                                       "B.cinerea", 
                                       "B.pseudocinerea", 
                                       "C.acutatum",
                                       "C.diplodiella",
                                       "D.ampelina",
                                       "E.ampelina",
                                       "G.uvicola",
                                       "G.bidwellii",
                                       "P.viticola",
                                       "Témoin")))

###################################

mock_abund<-otu_rel_abund %>%
  filter(level=="Species") %>%
  group_by(Ech, Type, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Ech, Type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = factor(taxon, 
                        levels=c("Botrytis cinerea complex",
                                 "Elsinoe ampelina",
                                 "Phyllosticta ampelicida",
                                 "Diaporthe ampelina complex",
                                 "Sclerotinia sclerotiorum complex",
                                 "Alternaria dauci",
                                 "Colletotrichum acutatum",
                                 "Plasmopara pusilla",
                                 "Neocercospora carotae complex",
                                 "Greeneria uvicola",
                                 "Coniella diplodiella",
                                 "Other")))

summary <- mock_abund %>% 
  group_by(Ech, Type, taxon) %>%
  summarize(mean_rel_abund) %>% 
  pivot_wider(names_from = Type, values_from = mean_rel_abund)


write.csv(summary, "Data/summary.csv")


#write.csv(summary, "Data/mocks/summary_mock.csv", row.names=FALSE)
library("ggplot2")
# plot the stacked bar chart 
CIMDEC_run_stacked<-ggplot(data=mock_abund, 
                           aes(x=Ech, 
                               y=mean_rel_abund, 
                               fill=taxon)) +
  geom_col(colour = "black", width=0.8, linewidth=0.1) +
  facet_wrap(vars(Type), nrow = 1, scales = "free")+
  #theme(strip.text = element_blank())+
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=get_cols(11))+
  theme(legend.position="bottom")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=5))+
  #theme(legend.position = c(0.9, -0.05),
  #      legend.justification = c(0.9, -0.05))+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=0.4));CIMDEC_run_stacked


################################################

library("cowplot")
depth<-plot_grid(read_depth,
                 CIMDEC_run_stacked, 
                 align="vh",
                 labels = c("A", "B"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2,
                 rel_heights = c(1, 2));depth


library("ggpubr")
ggsave(file="Figures/Barchart_run_2.pdf", 
       width=9, height=6, units="in", dpi=300)

#########################################################################################
#### pour run 1 - verifier la contamination du temoin, de G. bidwelli et P. viticola ####
#########################################################################################

# Temoin
Temoin <- subset_samples(CIMDEC_run_rarF, Ech =="Témoin")

summarize_phyloseq(Temoin)

A1<-plot_tree(tax_glom(Temoin, 
                       taxrank="Genus"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              label.tips="Genus", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              size="abundance",
              plot.margin =0.9, sizebase=1)+
  scale_size_continuous(range = c(0.0000001, 4))+
  facet_wrap(~Ech, scales="free_x")+ 
  guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(1, 1, 1, 1),
    legend.key.size = unit(0.2, 'cm'),
    legend.background = element_rect(fill = "transparent"));A1

###### G.bidwellii

G.bidwellii <- subset_samples(CIMDEC_run_rarF, Ech =="G.bidwellii")

B1<-plot_tree(tax_glom(G.bidwellii, 
                       taxrank="Genus"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              label.tips="Genus", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              size="abundance",
              plot.margin =0.9,
              sizebase=1)+
  scale_size_continuous(range = c(0.0000001, 4))+
  facet_wrap(~Ech, scales="free_x")+ 
  guides(color = guide_legend(nrow = 2)) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(1, 1, 1, 1),
    legend.key.size = unit(0.2, 'cm'),
    legend.background = element_rect(fill = "transparent")
  );B1


######

P.viticola <- subset_samples(CIMDEC_run_rarF, Ech =="P.viticola")

C1<-plot_tree(tax_glom(P.viticola, 
                       taxrank="Genus"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              label.tips="Genus", 
              text.size=3, 
              base.spacing=0.1,
              justify="jagged",
              size="abundance",
              plot.margin =0.9,
              sizebase=1)+
scale_size_continuous(range = c(0.0000001, 4))+
  facet_wrap(~Ech, scales="free_x")+
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(1, 1, 1, 1),
    legend.key.size = unit(0.2, 'cm'),
    legend.background = element_rect(fill = "transparent")
  );C1

######

library("cowplot")

Tree_1<-plot_grid(A1,
                  B1, 
                  C1, 
                  align="vh",
                  labels = c("A", "B", "C"),
                  hjust = -1,
                  vjust= 2,
                  nrow = 1);Tree_1

ggsave(file="Figures/arbre_run1.pdf", 
       width=7, height=9, units="in", dpi=300)

