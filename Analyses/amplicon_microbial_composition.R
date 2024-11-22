# Microbial Community Composition Analysis - 30 bald Amplicon sequencing
# Purpose: Reading in processed amplicon sequencing data and quanitifying microbial community diversity and composition 
# author: Joshua Fowler
# Date: Feb 28, 2024
library(tidyverse)
library(GUniFrac)
library(iNEXT)
library(vegan)
library(ape)
library(readxl)
# library(mia)


path <- c("~/Dropbox/UofMiami/Archbold 30 Bald Experiment - Processed Amplicon Sequencing Data")
taxonomy <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession") 


# loading in a key that connects exttraction tube id to our soil source sample id
extraction_key <- readxl::read_xlsx(path = paste0(path,"/30-Bald-libraryprep_samplepooling.xlsx"), sheet = "bald_id")

# loading in assigned taxonomy based on unite classifier for fungi and silva classifier for bacteria
taxonomy_ITS.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_ITS.tsv" )) %>% 
  mutate(clustering = "97", amplicon = "ITS") %>% rename(feature = `Feature ID`)
taxonomy_ITS.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_ITS.tsv" ))%>% 
  mutate(clustering = "99", amplicon = "ITS") %>% rename(feature = `Feature ID`)
taxonomy_ITS.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_ITS.tsv" ))%>% 
  mutate(clustering = "ESV", amplicon = "ITS") %>% rename(feature = `Feature ID`)
taxonomy_16S.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_16S.tsv" ))%>% 
  mutate(clustering = "97", amplicon = "16S") %>% rename(feature = `Feature ID`)
taxonomy_16S.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_16S.tsv" ))%>% 
  mutate(clustering = "99", amplicon = "16S") %>% rename(feature = `Feature ID`)
taxonomy_16S.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_16S.tsv" ))%>% 
  mutate(clustering = "ESV", amplicon = "16S") %>% rename(feature = `Feature ID`)





# loading in the feature tables showing OTU id's for each sample
featuretable_ITS.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-97-feature-table-ITS.tsv" ), skip = 1)
featuretable_ITS.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-99-feature-table-ITS.tsv" ), skip = 1)
featuretable_ITS.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/feature-table-ESV-ITS.tsv" ))


featuretable_16S.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-97-feature-table-16S.tsv" ), skip = 1)
featuretable_16S.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-99-feature-table-16S.tsv" ), skip = 1)
featuretable_16S.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/feature-table-ESV-16S.tsv" ))

# defining some functions to make calculating rarefaction curves and plotting easing for all the samples
community.matrix <- function(x){
  id <- x[[1]]
  s <- t(as.matrix(x[,sapply(x, is.numeric)]))
  colnames(s) <- id
  return(s)
  }
 
rarefaction_plot <-function(x){
  x %>% 
    group_by(Site) %>% 
    mutate(label = if_else(Sample == max(Sample), as.character(Label),NA_character_ )) %>% 
    ggplot()+
    geom_line(aes(x = Sample, y = Species, group = Sequencing_ID, color = as.factor(extraction_round)))+
    ggrepel::geom_text_repel(aes(x = Sample, y = Species, label = label, color = as.factor(extraction_round)), min.segment.length = unit(0, 'lines'), na.rm = TRUE)+
    guides(color = "none")+
    labs(x = "Sequencing Depth", y = "Number of OTUs") + theme_light()
}


# calculating rarefaction curves
rarecurve_16S.ESV_raw <- rarecurve(community.matrix(featuretable_16S.ESV), tidy = TRUE) 
rarecurve_16S.ESV <- rarecurve_16S.ESV_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))

rarecurve_16S.99_raw <- rarecurve(community.matrix(featuretable_16S.99), tidy = TRUE) 
rarecurve_16S.99 <- rarecurve_16S.99_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))

rarecurve_16S.97_raw <- rarecurve(community.matrix(featuretable_16S.97), tidy = TRUE) 
rarecurve_16S.97 <- rarecurve_16S.97_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))


rarecurve_ITS.ESV_raw <- rarecurve(community.matrix(featuretable_ITS.ESV), tidy = TRUE)
rarecurve_ITS.ESV <- rarecurve_ITS.ESV_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))

rarecurve_ITS.99_raw <- rarecurve(community.matrix(featuretable_ITS.99), tidy = TRUE) 
rarecurve_ITS.99 <- rarecurve_ITS.99_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))

rarecurve_ITS.97_raw <- rarecurve(community.matrix(featuretable_ITS.97), tidy = TRUE) 
rarecurve_ITS.97 <- rarecurve_ITS.97_raw %>%
  left_join(extraction_key, by = join_by(Site == Sequencing_ID)) %>% mutate(Sequencing_ID = Site, Label = case_when(!is.na(bald) ~ bald, is.na(bald) ~ Notes))


# plotting for each OTU clustering level
plot_rare_16S.ESV <- rarecurve_16S.ESV %>%  rarefaction_plot
plot_rare_16S.99 <- rarecurve_16S.99 %>%  rarefaction_plot
plot_rare_16S.97 <- rarecurve_16S.97 %>%  rarefaction_plot

# plot_rare_16S.ESV <- rarecurve_16S.ESV %>% filter(extraction_round == 2) %>%  rarefaction_plot
# plot_rare_16S.ESV <- rarecurve_16S.ESV %>% filter(extraction_round == 3) %>%  rarefaction_plot


plot_rare_ITS.ESV <- rarecurve_ITS.ESV %>%  rarefaction_plot
plot_rare_ITS.99 <- rarecurve_ITS.99 %>%  rarefaction_plot
plot_rare_ITS.97 <- rarecurve_ITS.97 %>%  rarefaction_plot

plot_rare_16S.ESV 
plot_rare_16S.99 
plot_rare_16S.97 

plot_rare_ITS.ESV
plot_rare_ITS.99 
plot_rare_ITS.97 


# calculating rarefied datasets to 2000 sequencing

rare_16S.ESV <- Rarefy(community.matrix(featuretable_16S.ESV), depth = 2500) 
rare_16S.99 <- Rarefy(community.matrix(featuretable_16S.99), depth = 2500) 
rare_16S.97 <- Rarefy(community.matrix(featuretable_16S.97), depth = 2500) 

rare_ITS.ESV <- Rarefy(community.matrix(featuretable_ITS.ESV), depth = 1000) 
rare_ITS.99 <- Rarefy(community.matrix(featuretable_ITS.99), depth = 1000) 
rare_ITS.97 <- Rarefy(community.matrix(featuretable_ITS.97), depth = 1000) 

# filtering out species that are found in the control samples
sample_16S_list <- extraction_key %>% filter(taxa == "16S", !(Sequencing_ID %in% c("16S-Control","16S-53-control", "16S-16-1", "16S-Cont-Apr2", "16S-53-PCRcontrol", "16S-Control2","16S-53", "16S-53A", "16S-53C-Mar29", "16S-OGC"))) %>% distinct(Sequencing_ID) %>% deframe() 



sample_ITS_list <- extraction_key %>% filter(taxa == "ITS", !(Sequencing_ID %in% c("ITS-Control","ITS-53-control", "ITS-16-1", "ITS-Cont-Apr2", "ITS-53-PCRcontrol", "ITS-Control2","ITS-53", "ITS-53A", "ITS-53C-Mar29", "ITS-OGC", "ITS2-ITS2Mar29"))) %>% distinct(Sequencing_ID) %>% deframe() 


cleaned_16S.ESV <- subset(rare_16S.ESV$otu.tab.rff, rownames(rare_16S.ESV$otu.tab.rff) %in% sample_16S_list)
cleaned_16S.99 <- subset(rare_16S.99$otu.tab.rff, rownames(rare_16S.99$otu.tab.rff) %in% sample_16S_list)
cleaned_16S.97 <- subset(rare_16S.97$otu.tab.rff, rownames(rare_16S.97$otu.tab.rff) %in% sample_16S_list)


cleaned_ITS.ESV <- subset(rare_ITS.ESV$otu.tab.rff, rownames(rare_ITS.ESV$otu.tab.rff) %in% sample_ITS_list)
cleaned_ITS.99 <- subset(rare_ITS.99$otu.tab.rff, rownames(rare_ITS.99$otu.tab.rff) %in% sample_ITS_list)
cleaned_ITS.97 <- subset(rare_ITS.97$otu.tab.rff, rownames(rare_ITS.97$otu.tab.rff) %in% sample_ITS_list)



# cleaned_16S.ESV <- rare_16S.ESV$otu.tab.rff[,rare_16S.ESV$otu.tab.rff["16S-Control",]!=0]
# cleaned_16S.ESV <- cleaned_16S.ESV[sample_16S_list,]
# 
# cleaned_16S.99<- rare_16S.99$otu.tab.rff[,rare_16S.99$otu.tab.rff["16S-Control",]!=0]
# cleaned_16S.99 <- cleaned_16S.99[sample_16S_list,]
# 
# cleaned_16S.97 <- rare_16S.97$otu.tab.rff[,rare_16S.97$otu.tab.rff["16S-Control",]==0]
# cleaned_16S.97 <- cleaned_16S.97[sample_16S_list,]
# 

# cleaned_ITS.ESV <- rare_ITS.ESV$otu.tab.rff[,rare_ITS.ESV$otu.tab.rff["ITS-Control",]==0]
# cleaned_ITS.ESV <- cleaned_ITS.ESV[sample_ITS_list,]
# 
# cleaned_ITS.99<- rare_ITS.99$otu.tab.rff[,rare_ITS.99$otu.tab.rff["ITS-Control",]==0]
# cleaned_ITS.99 <- cleaned_ITS.99[sample_ITS_list,]
# 
# cleaned_ITS.97 <- rare_ITS.97$otu.tab.rff[,rare_ITS.97$otu.tab.rff["ITS-Control",]==0]
# cleaned_ITS.97 <- cleaned_ITS.97[sample_ITS_list,]

dist_16S.ESV <- vegdist(cleaned_16S.ESV, method = "jaccard", binary=TRUE)
dist_16S.99 <- vegdist(cleaned_16S.99, method = "jaccard", binary=TRUE)
dist_16S.97 <- vegdist(cleaned_16S.97, method = "jaccard", binary=TRUE)


dist_ITS.ESV <- vegdist(cleaned_ITS.ESV, method = "jaccard", binary=TRUE)
dist_ITS.99 <- vegdist(cleaned_ITS.99, method = "jaccard", binary=TRUE)
dist_ITS.97 <- vegdist(cleaned_ITS.97, method = "jaccard", binary=TRUE)



# performing a PCoA to characteris the composition of each soil sources' microbial community
PCOA_16S.ESV <- wcmdscale(d = dist_16S.ESV, eig = T)
PCOA_16S.99 <- wcmdscale(d = dist_16S.99, eig = T)
PCOA_16S.97 <- wcmdscale(d = dist_16S.97, eig = T)

PCOA_ITS.ESV <- wcmdscale(d = dist_ITS.ESV, eig = T)
PCOA_ITS.99 <- wcmdscale(d = dist_ITS.99, eig = T)
PCOA_ITS.97 <- wcmdscale(d = dist_ITS.97, eig = T)

print(PCOA_16S.ESV)
plot(PCOA_16S.ESV)
plot(PCOA_16S.99)
plot(PCOA_16S.97)

plot(PCOA_ITS.ESV)
plot(PCOA_ITS.99)
plot(PCOA_ITS.97)



# saving first two pcoa axes for each plot

axes_16S <- as_tibble(PCOA_16S.ESV$points, rownames = NA) %>% 
  rownames_to_column(var = "sample_id") %>% mutate(amplicon = substr(sample_id, 1, 3), 
                                                   soil_tag = substr(sample_id, 5, 100)) %>% 
  mutate(soil_source = case_when(soil_tag == "94-Mar29"~"94",
                                 soil_tag == "16-0"~"16",
                                 TRUE ~ soil_tag)) %>% 
  dplyr::select(sample_id, amplicon, soil_source, Dim1, Dim2)

axes_ITS <- as_tibble(PCOA_ITS.ESV$points, rownames = NA) %>% 
  rownames_to_column(var = "sample_id") %>% mutate(amplicon = substr(sample_id, 1, 3), 
                                                   soil_tag = substr(sample_id, 5, 100)) %>% 
  mutate(soil_source = case_when(soil_tag == "24N-Apr2"~"24N",
                                 soil_tag == "6E-ITS46EMar29"~"6E",
                                 soil_tag == "60-ITS160Apr2"~"160",
                                 soil_tag == "ITS46E-ITS46EMar29"~"46E",
                                 soil_tag == "ITS97-ITS97Mar29"~"97",
                                 soil_tag == "ITS5E-ITS5EApr2"~"5E", 
                                 soil_tag == "3C-ITS53CMar29" ~ "53",
                                 soil_tag == "E-ITS5EApr2" ~ "5E",
                                 soil_tag == "7-ITS97Mar29" ~ "97",
                                 soil_tag == "65-Mar29-8ul" ~ "65",
                                 TRUE ~ soil_tag)) %>% 
  dplyr::select(sample_id, amplicon,  soil_source, Dim1, Dim2)
  



# Now I'm gonna calculate shannon diversity for each sample

get_diversity <- function(x, clustering){
  data.frame(sample_id = rownames(x),
             shannon_diversity = diversity(x),
             richness = apply(x>0, 1, sum),
             clustering = clustering) %>% 
    mutate(amplicon = substr(sample_id, 1, 3), 
           soil_tag = substr(sample_id, 5, 100)) 
}

shannon_16S.ESV <- get_diversity(cleaned_16S.ESV, clustering = "ESV")
shannon_16S.99 <- get_diversity(cleaned_16S.99, clustering = "99")
shannon_16S.97 <- get_diversity(cleaned_16S.97, clustering = "97")

shannon_16S <- bind_rows(shannon_16S.ESV, shannon_16S.99, shannon_16S.97) %>% 
  mutate(soil_source = case_when(soil_tag == "94-Mar29"~"94",
                                 soil_tag == "16-0"~"16",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)

shannon_ITS.ESV <- get_diversity(cleaned_ITS.ESV, clustering = "ESV")
shannon_ITS.99 <- get_diversity(cleaned_ITS.99, clustering = "99")
shannon_ITS.97 <- get_diversity(cleaned_ITS.97, clustering = "97") 

shannon_ITS <- bind_rows(shannon_ITS.ESV, shannon_ITS.99, shannon_ITS.97) %>% 
  mutate(soil_source = case_when(soil_tag == "24N-Apr2"~"24N",
                                 soil_tag == "6E-ITS46EMar29"~"6E",
                                 soil_tag == "60-ITS160Apr2"~"160",
                                 soil_tag == "ITS46E-ITS46EMar29"~"46E",
                                 soil_tag == "ITS97-ITS97Mar29"~"97",
                                 soil_tag == "ITS5E-ITS5EApr2"~"5E", 
                                 soil_tag == "3C-ITS53CMar29" ~ "53",
                                 soil_tag == "E-ITS5EApr2" ~ "5E",
                                 soil_tag == "7-ITS97Mar29" ~ "97",
                                 soil_tag == "65-Mar29-8ul" ~ "65",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)
  




composition_16S <- axes_16S %>% 
  left_join(shannon_16S)


composition_ITS <- axes_ITS %>% 
  left_join(shannon_ITS)

composition_df <- bind_rows(composition_16S, composition_ITS) %>% 
  select(-sample_id) %>% 
  pivot_wider(names_from = amplicon, values_from = c(Dim1, Dim2, richness, shannon_diversity))

write_csv(composition_df, "30bald_composition_metrics.csv")

# a few plots to look at differences across sites and clustering
composition_long <-  bind_rows(composition_16S, composition_ITS) %>% 
  select(-sample_id) %>% 
  pivot_longer(cols = c(Dim1, Dim2, shannon_diversity, richness), names_to = "metric", values_to = "value")

ggplot(composition_long %>% filter(metric == "richness")) +
  geom_tile(aes(x = soil_source, y = clustering, fill = value)) +
  facet_grid(amplicon~metric, scales = "free")

ggplot(composition_long %>% filter(metric == "shannon_diversity")) +
  geom_tile(aes(x = soil_source, y = clustering, fill = value)) +
  facet_grid(amplicon~metric, scales = "free")

ggplot(composition_long %>% filter(metric == "shannon_diversity")) +
  geom_tile(aes(x = soil_source, y = clustering, fill = value)) +
  facet_grid(amplicon~metric, scales = "free")

ggplot(composition_long %>% filter(metric == "Dim1")) +
  geom_tile(aes(x = soil_source, y = clustering, fill = value)) +
  facet_grid(amplicon~metric, scales = "free")

ggplot(composition_long %>% filter(metric == "Dim2")) +
  geom_tile(aes(x = soil_source, y = clustering, fill = value)) +
  facet_grid(amplicon~metric, scales = "free")



# vizualizing the taxonomy across the samples
join_taxonomy <- function(table, taxa, clustering){
  as_tibble(table, rownames = "sample_id") %>% 
    mutate(soil_tag = substr(sample_id, 5, 100)) %>% 
    pivot_longer(cols = c(-sample_id, -soil_tag), names_to = "feature", values_to = "reads") %>% 
    left_join(taxa, relationship = "many-to-many") 
}

cleaned_taxa_16S.ESV <- join_taxonomy(cleaned_16S.ESV, taxonomy_16S.ESV) %>% 
  mutate(soil_source = case_when(soil_tag == "94-Mar29"~"94",
                                 soil_tag == "16-0"~"16",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)
cleaned_taxa_16S.99 <- join_taxonomy(cleaned_16S.99, taxonomy_16S.99)%>% 
  mutate(soil_source = case_when(soil_tag == "94-Mar29"~"94",
                                 soil_tag == "16-0"~"16",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)
cleaned_taxa_16S.97 <- join_taxonomy(cleaned_16S.97, taxonomy_16S.97) %>% 
  mutate(soil_source = case_when(soil_tag == "94-Mar29"~"94",
                                 soil_tag == "16-0"~"16",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)

cleaned_taxa_ITS.ESV <- join_taxonomy(cleaned_ITS.ESV, taxonomy_ITS.ESV)%>% 
  mutate(soil_source = case_when(soil_tag == "24N-Apr2"~"24N",
                                 soil_tag == "6E-ITS46EMar29"~"6E",
                                 soil_tag == "60-ITS160Apr2"~"160",
                                 soil_tag == "ITS46E-ITS46EMar29"~"46E",
                                 soil_tag == "ITS97-ITS97Mar29"~"97",
                                 soil_tag == "ITS5E-ITS5EApr2"~"5E", 
                                 soil_tag == "3C-ITS53CMar29" ~ "53",
                                 soil_tag == "E-ITS5EApr2" ~ "5E",
                                 soil_tag == "7-ITS97Mar29" ~ "97",
                                 soil_tag == "65-Mar29-8ul" ~ "65",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)
cleaned_taxa_ITS.99 <- join_taxonomy(cleaned_ITS.99, taxonomy_ITS.99)%>% 
  mutate(soil_source = case_when(soil_tag == "24N-Apr2"~"24N",
                                 soil_tag == "6E-ITS46EMar29"~"6E",
                                 soil_tag == "60-ITS160Apr2"~"160",
                                 soil_tag == "ITS46E-ITS46EMar29"~"46E",
                                 soil_tag == "ITS97-ITS97Mar29"~"97",
                                 soil_tag == "ITS5E-ITS5EApr2"~"5E", 
                                 soil_tag == "3C-ITS53CMar29" ~ "53",
                                 soil_tag == "E-ITS5EApr2" ~ "5E",
                                 soil_tag == "7-ITS97Mar29" ~ "97",
                                 soil_tag == "65-Mar29-8ul" ~ "65",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)
cleaned_taxa_ITS.97 <- join_taxonomy(cleaned_ITS.97, taxonomy_ITS.97)%>% 
  mutate(soil_source = case_when(soil_tag == "24N-Apr2"~"24N",
                                 soil_tag == "6E-ITS46EMar29"~"6E",
                                 soil_tag == "60-ITS160Apr2"~"160",
                                 soil_tag == "ITS46E-ITS46EMar29"~"46E",
                                 soil_tag == "ITS97-ITS97Mar29"~"97",
                                 soil_tag == "ITS5E-ITS5EApr2"~"5E", 
                                 soil_tag == "3C-ITS53CMar29" ~ "53",
                                 soil_tag == "E-ITS5EApr2" ~ "5E",
                                 soil_tag == "7-ITS97Mar29" ~ "97",
                                 soil_tag == "65-Mar29-8ul" ~ "65",
                                 TRUE ~ soil_tag)) %>% select(-soil_tag)






taxonomy_df <- bind_rows(cleaned_taxa_ITS.97, cleaned_taxa_ITS.99, cleaned_taxa_ITS.ESV,
                         cleaned_taxa_16S.97, cleaned_taxa_16S.99, cleaned_taxa_16S.ESV) %>% 
  separate(Taxon, into = c(taxonomy), sep = ";") %>% 
  mutate(across(all_of(taxonomy), ~ word(.x, 2, sep = fixed("__"))))

abund_16S.ESV <- taxonomy_df %>% filter(amplicon == "16S" & clustering == "ESV") %>% 
  filter(reads>0)

abund_ITS.ESV <- taxonomy_df %>% filter(amplicon == "ITS" & clustering == "ESV") %>% 
  filter(reads>0)



ggplot(abund_16S.ESV)+
  geom_bar(aes(x = soil_source, y = reads, fill = Phylum), stat = "identity")+
  theme_minimal() +theme(axis.text.x = element_text(hjust = 1, angle = 45))


ggplot(abund_ITS.ESV)+
  geom_bar(aes(x = soil_source, y = reads, fill = Phylum), stat = "identity")+
  theme_minimal() + theme(axis.text.x = element_text(hjust = 1, angle = 45))



# plotting just the top 5 most abundant
top5_16S.ESV <- abund_16S.ESV %>% 
  group_by(soil_source, sample_id) %>% 
  slice_max(order_by = reads, n = 5)



top5_ITS.ESV <- abund_ITS.ESV %>% 
  group_by(soil_source, sample_id) %>% 
  slice_max(order_by = reads, n = 5)

ggplot(top5_16S.ESV)+
  geom_bar(aes(x = soil_source, y = reads, fill = Family), stat = "identity")+
  theme_minimal() +theme(axis.text.x = element_text(hjust = 1, angle = 45))


ggplot(top5_ITS.ESV)+
  geom_bar(aes(x = soil_source, y = reads, fill = Family), stat = "identity")+
  theme_minimal() + theme(axis.text.x = element_text(hjust = 1, angle = 45))












# old stuff





# visualizing pca loadings

propvar_16S.ESV <- summary(PCA_16S.ESV)$cont$importance[2,]
cumvar_16S.ESV <- summary(PCA_16S.ESV)$cont$importance[3,]

barplot(cumvar_16S.ESV, col = "red");abline(a = .8, b = 0)
barplot(propvar_16S.ESV, type = "p")



propvar_16S.97 <- summary(PCA_16S.97)$cont$importance[2,]
cumvar_16S.97 <- summary(PCA_16S.97)$cont$importance[3,]

barplot(cumvar_16S.97, col = "red");abline(a = .8, b = 0)
barplot(propvar_16S.97, type = "p")



propvar_ITS.ESV <- summary(PCA_ITS.ESV)$cont$importance[2,]
cumvar_ITS.ESV <- summary(PCA_ITS.ESV)$cont$importance[3,]

barplot(cumvar_ITS.ESV, col = "red");abline(a = .8, b = 0)
barplot(propvar_ITS.ESV, type = "p")


propvar_ITS.97 <- summary(PCA_ITS.97)$cont$importance[2,]
cumvar_ITS.97 <- summary(PCA_ITS.97)$cont$importance[3,]

barplot(cumvar_ITS.97, col = "red");abline(a = .8, b = 0)
barplot(propvar_ITS.97, type = "p")



barplot(PCA_16S.97$CA$eig)
ordiplot(PCA_16S.97, display = "sites", choices = c(1,2))
orditorp(PCA_ITS.97, display = "sites", choices = c(1,2))



summary(PCA_ITS.97)
plot(PCA_ITS.97$CA$eig, type = "l")
ordiplot(PCA_ITS.97, choices = c(4,3), display = "sites")
orditorp(PCA_ITS.97, choices = c(4,3), display = "sites")


barplot(PCA_16S.ESV$CA$eig)
evplot()



# Calculating diversity
div_16S.97 <- tibble(sample_id = names(diversity(cleaned_16S.97, index = "shannon")), diversity = diversity(cleaned_16S.97, index = "shannon"))

ggplot(div_16S.97)+
  geom_point(aes(x = sample_id, y = diversity))


div_ITS.97 <- tibble(sample_id = names(diversity(cleaned_ITS.97, index = "shannon")), diversity = diversity(cleaned_ITS.97, index = "shannon"))

ggplot(div_ITS.97)+
  geom_point(aes(x = sample_id, y = diversity))



# calculating proportion of different families
taxa_cleaned_16S.97 <- as.data.frame(t(cleaned_16S.97)) %>% rownames_to_column(., "OTU_id") %>% 
  mutate(clustering = "97", amplicon = "16S")

taxa_cleaned_ITS.97 <- as.data.frame(t(cleaned_ITS.97)) %>% rownames_to_column(., "OTU_id") %>% 
  mutate(clustering = "97", amplicon = "ITS")


  
taxa_cleaned <- bind_rows(taxa_cleaned_16S.97, taxa_cleaned_ITS.97) %>% 
  pivot_longer(cols = -c(OTU_id, clustering, amplicon), names_to = "Sample_id", values_to = "count") %>% 
  left_join(taxonomy_df, by = c("OTU_id" = "Feature ID", "clustering" = "clustering", "amplicon" = "amplicon")) %>% 
  mutate(Sample_id_fixed = word(Sample_id, 2, sep = "\\-"))



ggplot(taxa_cleaned)+
  geom_bar(aes(x = Sample_id_fixed, y = count, fill = Phylum), position = "fill", stat = "sum")+
  facet_wrap(~amplicon)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# looking at the contamination in the controls
dirty_16S.ESV <- rare_16S.ESV$otu.tab.rff[!rownames(rare_16S.ESV$otu.tab.rff) %in% sample_16S_list,]


taxa_dirty_16S.ESV <- as.data.frame(t(dirty_16S.ESV)) %>% rownames_to_column(., "OTU_id") %>% 
  mutate(clustering = "ESV", amplicon = "16S") %>% 
  pivot_longer(cols = -c(OTU_id, clustering, amplicon), names_to = "Sample_id", values_to = "count") %>% 
  left_join(taxonomy_df, by = c("OTU_id" = "Feature ID", "clustering" = "clustering", "amplicon" = "amplicon")) %>% 
  mutate(Sample_id_fixed = word(Sample_id, 2, sep = "\\-"))


ggplot(taxa_dirty_16S.ESV)+
  geom_bar(aes(x = Sample_id_fixed, y = count, fill = Phylum), position = "fill", stat = "sum")+
  facet_wrap(~amplicon)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
