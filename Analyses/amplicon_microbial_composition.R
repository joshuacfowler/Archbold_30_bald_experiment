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

# loading in assigned taxonomy based on unite classifier for fungi and silva classifier for bacteria
taxonomy_ITS.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_ITS.tsv" ))
taxonomy_ITS.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_ITS.tsv" ))
taxonomy_ITS.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_ITS.tsv" ))
taxonomy_16S.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_97_16S.tsv" ))
taxonomy_16S.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_99_16S.tsv" ))
taxonomy_16S.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/taxonomy_ESV_16S.tsv" ))

# loading in the feature tables showing OTU id's for each sample
featuretable_ITS.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-97-feature-table-ITS.tsv" ), skip = 1)
featuretable_ITS.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-99-feature-table-ITS.tsv" ), skip = 1)
featuretable_ITS.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/feature-table-ESV-ITS.tsv" ))


featuretable_16S.97 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-97-feature-table-16S.tsv" ), skip = 1)
featuretable_16S.99 <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/denovo-99-feature-table-16S.tsv" ), skip = 1)
featuretable_16S.ESV <- read_tsv(file = paste0(path,"/Processed_amplicon_microbiome_data", "/feature-table-ESV-16S.tsv" ))

featuretable_16S.ESV <- na.omit(as.matrix(featuretable_16S.ESV[,sapply(featuretable_16S.ESV, is.numeric) ]))


# defining some functions to make calculating rarefaction curves and plotting easing for all the samples
community.matrix <- function(x){t(as.matrix(x[,sapply(x, is.numeric)]))}
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
rare_16S.ESV <- rarecurve(community.matrix(featuretable_16S.ESV), se = TRUE, tidy = TRUE) 
rare_16S.99 <- rarecurve(community.matrix(featuretable_16S.99), tidy = TRUE) 
rare_16S.97 <- rarecurve(community.matrix(featuretable_16S.97), tidy = TRUE) 

rare_ITS.ESV <- rarecurve(community.matrix(featuretable_ITS.ESV), tidy = TRUE) 
rare_ITS.99 <- rarecurve(community.matrix(featuretable_ITS.99), tidy = TRUE) 
rare_ITS.97 <- rarecurve(community.matrix(featuretable_ITS.97), tidy = TRUE) 

# plotting for each OTU clustering level
plot_rare_16S.ESV <- rare_16S.ESV %>%  rarefaction_plot
plot_rare_16S.99 <- rare_16S.99 %>%  rarefaction_plot
plot_rare_16S.97 <- rare_16S.97 %>%  rarefaction_plot

plot_rare_ITS.ESV <- rare_ITS.ESV %>%  rarefaction_plot
plot_rare_ITS.99 <- rare_ITS.99 %>%  rarefaction_plot
plot_rare_ITS.97 <- rare_ITS.97 %>%  rarefaction_plot

plot_rare_16S.ESV
plot_rare_16S.99 
plot_rare_16S.97 

plot_rare_ITS.ESV
plot_rare_ITS.99 
plot_rare_ITS.97 


