###############################################################
########## preprocessing for the mock experiments #############
###############################################################

# Import the raw phyloseq mock_raw object 

CIMDEC_run_raw <- readRDS("R_objects/CIMDEC_run_raw.rds")

###############################################################################################
######### Summary stats ########

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) # devtools::install_github("microsud/microbiomeutilities")

summarize_phyloseq(CIMDEC_run_raw)

Dep1<-plot_read_distribution(CIMDEC_run_raw, groups = "Type", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10',limits=c(1000, 1000000))+
  scale_fill_manual(values=c("#111111"))+
  geom_vline(xintercept = 50000, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "", y = "Count");Dep1

CIMDEC_run_raw <- prune_samples(sample_sums(CIMDEC_run_raw) >= 0, CIMDEC_run_raw)

summarize_phyloseq(CIMDEC_run_raw)

Dep2<-plot_read_distribution(CIMDEC_run_raw, groups = "Type", 
                             plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits=c(1000, 1000000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 50000, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "Reads per samples", y = "Count");Dep2


library("cowplot")
depth<-plot_grid(Dep1+theme(legend.position="none"),
                 Dep2+theme(legend.position="none"), 
                 align="vh",
                 labels = c("A", "B"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2)

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final

library("ggpubr")
ggsave(file="Figures/FigS2_depth_final.pdf", 
       width=8, height=5, units="in", dpi=300)


##############################################################################################

library("phyloseq")
library("stringr")
taxM <- data.frame(phyloseq::tax_table(CIMDEC_run_raw))
tax.cleanM <- data.frame(row.names = row.names(taxM),
                         Kingdom = str_replace(taxM[,1], "k__",""),
                         Phylum = str_replace(taxM[,2], "p__",""),
                         Class = str_replace(taxM[,3], "c__",""),
                         Order = str_replace(taxM[,4], "o__",""),
                         Family = str_replace(taxM[,5], "f__",""),
                         Genus = str_replace(taxM[,6], "g__",""),
                         Species = str_replace(taxM[,7], "s__",""),
                         stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanM[,i] <- as.character(tax.cleanM[,i])}

tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "_", " "))

########################################################################

#Sclerotinia sclerotiorum complex
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Sclerotinia sclerotiorum", "Sclerotinia sclerotiorum complex"))%>%
  mutate(Species=str_replace(Species, "Sclerotinia trifoliorum", "Sclerotinia sclerotiorum complex"))



phyloseq::tax_table(CIMDEC_run_raw) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(CIMDEC_run_raw)

# For the mock experiments, the data were first rarefied
CIMDEC_run_rar = rarefy_even_depth(CIMDEC_run_raw, 
                                   rngseed=1024,
                                   replace=FALSE)

CIMDEC_run_rarF = prune_taxa(taxa_sums(CIMDEC_run_rar) > 9, CIMDEC_run_rar); CIMDEC_run_rarF




saveRDS(CIMDEC_run_rarF, "R_objects/CIMDEC_run_rarF.rds")
