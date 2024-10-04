# Microbial Community Composition Analysis - 30 bald Amplicon sequencing
# Purpose: Reading in processed amplicon sequencing data and quanitifying microbial community diversity and composition 
# author: Joshua Fowler
# Date: Feb 28, 2024
library(tidyverse)
library(GUniFrac)
library(iNEXT)
library(vegan)
# library(mia)


path <- c("~/Dropbox/UofMiami/Archbold 30 Bald Experiment - Processed Amplicon Sequencing Data")
taxonomy <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession") 

# loading in assigned taxonomy based on unite classifier for fungi and silva classifier for bacteria
taxonomy_ITS.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_ITS.tsv" )) %>% 
  mutate(clustering = "97", amplicon = "ITS")
taxonomy_ITS.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_ITS.tsv" ))%>% 
  mutate(clustering = "99", amplicon = "ITS")
taxonomy_ITS.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_ITS.tsv" ))%>% 
  mutate(clustering = "ESV", amplicon = "ITS")
taxonomy_16S.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_16S.tsv" ))%>% 
  mutate(clustering = "97", amplicon = "16S")
taxonomy_16S.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_16S.tsv" ))%>% 
  mutate(clustering = "99", amplicon = "16S")
taxonomy_16S.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_16S.tsv" ))%>% 
  mutate(clustering = "ESV", amplicon = "16S")



taxonomy_df <- bind_rows(taxonomy_ITS.97, taxonomy_ITS.99, taxonomy_ITS.ESV,
                        taxonomy_16S.97, taxonomy_16S.99, taxonomy_16S.ESV) %>% 
  separate(Taxon, into = c(taxonomy), sep = ";") %>% mutate(across(taxonomy, ~ word(.x, 2, sep = fixed("__"))))



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
    mutate(label = if_else(Sample == max(Sample), as.character(Site),NA_character_ )) %>% 
    ggplot()+
    geom_line(aes(x = Sample, y = Species, group = Site, color = Site))+
    ggrepel::geom_text_repel(aes(x = Sample, y = Species, label = label), min.segment.length = unit(0, 'lines'), na.rm = TRUE)+
    guides(color = "none")+
    labs(x = "Sequencing Depth", y = "Number of OTUs") + theme_light()
}

# calculating rarefaction curves
rarecurve_16S.ESV <- rarecurve(community.matrix(featuretable_16S.ESV), tidy = TRUE) 
rarecurve_16S.99 <- rarecurve(community.matrix(featuretable_16S.99), tidy = TRUE) 
rarecurve_16S.97 <- rarecurve(community.matrix(featuretable_16S.97), tidy = TRUE) 

rarecurve_ITS.ESV <- rarecurve(community.matrix(featuretable_ITS.ESV), tidy = TRUE) 
rarecurve_ITS.99 <- rarecurve(community.matrix(featuretable_ITS.99), tidy = TRUE) 
rarecurve_ITS.97 <- rarecurve(community.matrix(featuretable_ITS.97), tidy = TRUE) 

# plotting for each OTU clustering level
plot_rare_16S.ESV <- rarecurve_16S.ESV %>%  rarefaction_plot
plot_rare_16S.99 <- rarecurve_16S.99 %>%  rarefaction_plot
plot_rare_16S.97 <- rarecurve_16S.97 %>%  rarefaction_plot

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

rare_16S.ESV <- Rarefy(community.matrix(featuretable_16S.ESV), depth = 2000) 
rare_16S.99 <- Rarefy(community.matrix(featuretable_16S.99), depth = 2000) 
rare_16S.97 <- Rarefy(community.matrix(featuretable_16S.97), depth = 2000) 

rare_ITS.ESV <- Rarefy(community.matrix(featuretable_ITS.ESV), depth = 1000) 
rare_ITS.99 <- Rarefy(community.matrix(featuretable_ITS.99), depth = 1000) 
rare_ITS.97 <- Rarefy(community.matrix(featuretable_ITS.97), depth = 1000) 

# filtering out species that are found in the control samples
sample_16S_list <- c( "16S-16-0","16S-16-1","16S-1N","16S-1S","16S-2","16S-21N","16S-21S","16S-24N",
                  "16S-24S","16S-26","16S-28","16S-35","16S-36","16S-37","16S-38","16S-39",
                  "16S-41N","16S-45S","16S-46E","16S-49",
                  "16S-53C-Mar29","16S-5E","16S-62","16S-65","16S-7N","16S-84","16S-94-Mar29","16S-95",
                  "16S-97","16S-99")

sample_ITS_list<-c("ITS-16-1","ITS-1N","ITS-21N","ITS-21S","ITS-24N-Apr2","ITS-24S","ITS-26","ITS-28",
                   "ITS-36","ITS-37","ITS-38","ITS-39","ITS-41N","ITS-45S","ITS-49",
                   "ITS-53A","ITS-53B","ITS-62","ITS-65-Mar29-8ul","ITS-7N","ITS-84","ITS-94","ITS-95",
                   "ITS160-ITS160Apr2","ITS2-ITS2Mar29","ITS46E-ITS46EMar29",
                   "ITS53C-ITS53CMar29","ITS5E-ITS5EApr2","ITS97-ITS97Mar29")

cleaned_16S.ESV <- rare_16S.ESV$otu.tab.rff[,rare_16S.ESV$otu.tab.rff["16S-Control",]==0]
cleaned_16S.ESV <- cleaned_16S.ESV[sample_16S_list,]

cleaned_16S.99<- rare_16S.99$otu.tab.rff[,rare_16S.99$otu.tab.rff["16S-Control",]==0]
cleaned_16S.99 <- cleaned_16S.99[sample_16S_list,]

cleaned_16S.97 <- rare_16S.97$otu.tab.rff[,rare_16S.97$otu.tab.rff["16S-Control",]==0]
cleaned_16S.97 <- cleaned_16S.97[sample_16S_list,]


cleaned_ITS.ESV <- rare_ITS.ESV$otu.tab.rff[,rare_ITS.ESV$otu.tab.rff["ITS-Control",]==0]
cleaned_ITS.ESV <- cleaned_ITS.ESV[sample_ITS_list,]

cleaned_ITS.99<- rare_ITS.99$otu.tab.rff[,rare_ITS.99$otu.tab.rff["ITS-Control",]==0]
cleaned_ITS.99 <- cleaned_ITS.99[sample_ITS_list,]

cleaned_ITS.97 <- rare_ITS.97$otu.tab.rff[,rare_ITS.97$otu.tab.rff["ITS-Control",]==0]
cleaned_ITS.97 <- cleaned_ITS.97[sample_ITS_list,]


PCA_16S.ESV <- rda(cleaned_16S.ESV)
PCA_16S.99 <- rda(cleaned_16S.99)
PCA_16S.97 <- rda(cleaned_16S.97)


PCA_ITS.ESV <- rda(rare_ITS.ESV$otu.tab.rff)
PCA_ITS.99 <- rda(rare_ITS.99$otu.tab.rff)
PCA_ITS.97 <- rda(rare_ITS.97$otu.tab.rff)



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
