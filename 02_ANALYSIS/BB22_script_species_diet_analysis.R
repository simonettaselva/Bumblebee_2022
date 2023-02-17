################################################
# GOAL 0: Species Diet Analysis
# by Simonetta Selva
#
# Created: December 14, 2022
# Project: Bumblebee 2022
################################################

# AIM: Characterize the diet compositional and structural (taxonomic, functional, and phylogenetic diversity) 
# and chemical properties of the two bumblebee species in both urban and rural landscapes.

# information: every subsection works in itself


# SPECIES ABBUNDANCES IN POLLEN AND PHYLOGENETIC TREE ----

## preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pals)
library(V.PhyloMaker)
library(ape)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set with 

# initialize region as column
BB22.full$region <- c() 

# add site and region as columns
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <-paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep="")
  BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep="")
}

## 1. create a data frame with taxa (species, genus, family) ----

  phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

## 2. phylogentic tree ----

  # phylogenetic hypotheses under three scenarios based on a backbone phylogeny 
  tree.result <- phylo.maker(phylo, scenarios=c("S1","S2","S3"))
  
  # plot the phylogenies with node ages displayed with sceanario 3
  tree <- plot.phylo(tree.result$scenario.3, cex = 0.5, main = "Phylogenetic tree of species in pollen"); tree
  
  # get order of species in tree
  phylo.order <- data.frame(sps=tree.result$scenario.3$tip.label)
  phylo.order$order <- seq(1, length(phylo.order$sps))
  colnames(phylo.order) <- c("plant.species", " order")
  phylo.order$plant.species <- sub("_", " ", phylo.order$plant.species)
  
## 3. find how many taxa were visited per region ----
  
  # create new data frame with needed variables
  BB22.full.bubble <- BB22.full%>%
    dplyr::group_by(site, plant.species, bbspecies)%>%
    dplyr::summarise(Abundance = sum(Abundance),
                     region = region,
                     landscape = landscape)
  
  # find how many taxa were visited per region
  BB22.full.taxa <- BB22.full.bubble %>%
    dplyr::group_by(region, bbspecies) %>%
    dplyr::summarise(plant.species = plant.species)%>%
    distinct()
  
  BB22.full.taxa <- BB22.full.taxa%>%
    dplyr::summarise(Nr.taxa = n())
  
  write_csv(BB22.full.taxa, "BB22_NrTaxa_region.csv")

## 4. plotting ----
  
  # merge two data frames and order along plant species in phylo-tree
  BB22.full.bubble_ordered <- merge(x = BB22.full.bubble, y = phylo.order, by.x = "plant.species")%>%
    mutate(plant.species = as_factor(plant.species))
  BB22.full.bubble_ordered <- BB22.full.bubble_ordered[order(BB22.full.bubble_ordered$` order`),]
  BB22.full.bubble_ordered$plant.species <- factor(BB22.full.bubble_ordered$plant.species, 
                                                   levels = unique(BB22.full.bubble_ordered$plant.species[order(BB22.full.bubble_ordered$` order`)]))

  # plot bubble plot with relative abundances along order of species in tree
  setwd(output)
  
  ### FIGURE 1 ####
  # region level
  palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
  ggplot(BB22.full.bubble_ordered, aes(x = region, y =BB22.full.bubble_ordered$plant.species)) + 
    geom_point(aes(size = Abundance, color = landscape), alpha=0.5) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species")+    
    theme_classic(base_size = 20) + 
    guides(alpha = "none") +
    scale_color_manual(values = palette.landscape, labels = c("rural", "urban"), name = "Landscape") +
    guides(color = guide_legend(override.aes=list(alpha = 1)))
  ggsave(paste("./01_Goal 0/FIGURE_1.png", sep = ""), width = 16, height = 16)
  
  # on site level (not used)
  palette.site <- kelly(18)[3:18] #create color palette for sites
  ggplot(BB22.full.bubble_ordered, aes(x = site, y =BB22.full.bubble_ordered$plant.species, color = site)) + 
    geom_point(aes(size = Abundance, fill = site, alpha=0.5)) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species") +
    theme(axis.text.x = element_text(angle = 90))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_color_manual(values = palette.site, guide = "none")+ #no legend
    scale_fill_manual(values = palette.site, guide = "none") + #no legend
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # landscape level (not used)
  palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
  ggplot(BB22.full.bubble_ordered, aes(x = landscape, y =BB22.full.bubble_ordered$plant.species, color = landscape)) + 
    geom_point(aes(size = Abundance, fill = landscape, alpha=0.5)) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species")+ 
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_color_manual(values = palette.landscape, guide = "none") +
    scale_fill_manual(values = palette.landscape, guide = "none")
  setwd(input)
