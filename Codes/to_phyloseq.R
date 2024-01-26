#######################################################################
############## build the phyloseq object for the runs #################
#######################################################################

# BiocManager::install("phyloseq")

## load the biom_table, taxonomy and metadata ####
CIMDEC_biom<- read.csv("Data/CIMDEC_otu_table_CIMDEC_plaque_2.csv", 
                       header=TRUE, 
                       sep=",")
head(CIMDEC_biom)

CIMDEC_taxo<- read.csv("Data/CIMDEC_phyloseq_taxonomy_CIMDEC_plaque_2.csv", 
                       header=TRUE, 
                       sep=",")
head(CIMDEC_taxo)

CIMDEC_meta <- read.csv("Data/CIMDEC_Meta_data_plaque_2.csv", 
                        header=TRUE, 
                        sep=",") #\t
head(CIMDEC_meta)

library("dplyr")
CIMDEC_meta$year<-as.factor(CIMDEC_meta$year)

CIMDEC_meta <- CIMDEC_meta %>% 
  tibble::column_to_rownames("Sample_ID")

library("phyloseq")
samplesS = sample_data(CIMDEC_meta)

# define the row names from the otu column ####

CIMDEC_biom <- CIMDEC_biom %>%
  tibble::column_to_rownames("OTU") 

CIMDEC_taxo <- CIMDEC_taxo %>%
  tibble::column_to_rownames("OTU") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

CIMDEC_biom <- as.matrix(CIMDEC_biom)
CIMDEC_taxo <- as.matrix(CIMDEC_taxo)

class(CIMDEC_biom)
class(CIMDEC_taxo)

# convert to phyloseq objects ####

CIMDEC_OTU = otu_table(CIMDEC_biom, taxa_are_rows = TRUE)
CIMDEC_TAX = phyloseq::tax_table(CIMDEC_taxo)


CIMDEC_run_raw <- phyloseq(CIMDEC_OTU, 
                           CIMDEC_TAX,
                           samplesS)


CIMDEC_run_raw

sample_variables(CIMDEC_run_raw)

#create and merge 
library("ape")
random_tree = rtree(ntaxa(CIMDEC_run_raw), 
                    rooted=TRUE, 
                    tip.label=taxa_names(CIMDEC_run_raw))

CIMDEC_run_raw <- phyloseq(CIMDEC_OTU, 
                           CIMDEC_TAX, 
                           samplesS, 
                           random_tree)
CIMDEC_run_raw

saveRDS(CIMDEC_run_raw, file= "R_objects/CIMDEC_run_raw.rds")

###############################################################
########## preprocessing for the mock experiments #############
###############################################################

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

#Botrytis cinerea complex
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Botrytis cinerea", "Botrytis cinerea complex"))%>%
  mutate(Species=str_replace(Species, "Botrytis pseudocinerea", "Botrytis cinerea complex"))

#Neocercospora carotae complex
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Cercospora carotae", "Neocercospora carotae complex"))%>%
  mutate(Species=str_replace(Species, "Neocercospora carotae", "Neocercospora carotae complex"))

#Diaporthe ampelina complex
#souches a valider, identification pheno seulement - confirmer par sequencage si ampelina ou ambigua
#D. ambigua peut aussi causer des symptomes de phomopsis sur les bois et tissus verts 
#complex temporaire jusqu'a confirmation souches
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Diaporthe ampelina", "Diaporthe ampelina complex"))%>%
  mutate(Species=str_replace(Species, "Diaporthe ambigua", "Diaporthe ampelina complex"))


phyloseq::tax_table(CIMDEC_run_raw) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(CIMDEC_run_raw)

# Tous les ASVs qui ne sont pas supportes par un minimum de 10 lectures sont enleves du jeu de donnees

CIMDEC_run_rarF = prune_taxa(taxa_sums(CIMDEC_run_raw) > 9, CIMDEC_run_raw); CIMDEC_run_rarF

saveRDS(CIMDEC_run_rarF, "R_objects/CIMDEC_run_rarF.rds")


