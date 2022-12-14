################################################
# GOAL 0
# by Simonetta Selva
#
# Created: December 14, 2022
# Project: Bumblebee 2022
################################################

# AIM: Characterize the diet compositional and structural (taxonomic, functional, and phylogenetic diversity) 
# and chemical properties of the two bumblebee species in both urban and rural landscapes.

#reset environment
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

############################
#### POLLEN COMPOSITION ####


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

#### plot plant families per site, region and species

# 1. find the 30 most abundant families
  families.overview <- BB22.full%>%
    dplyr::group_by(family)%>%
    dplyr::summarise(cum.abund = sum(Abundance))
  tot <- sum(families.overview$cum.abund)
  families.overview <- families.overview %>%
    dplyr::summarise(family = family,
                     cum.abund = cum.abund,
                     rel.abund = cum.abund/tot)
  # plot family abundance
  ggplot(families.overview, aes(y=reorder(family, -cum.abund), x=cum.abund)) + 
    geom_bar(stat="identity")+ theme_classic()
  
  # use cumulative abundance to select 30 most abundant species
  rare.families <- families.overview %>% top_n(nrow(families.overview)-30, -cum.abund)
  BB22.full$family.agg <- BB22.full$family
  # group all rare families into "other families"
  for(h in rare.families$family){
    for(i in 1:nrow(BB22.full)){
      if(BB22.full$family.agg[i] == h){
        BB22.full$family.agg[i] <- "Other families"
      }
    }
  }

# 2. produce data frame with families' abundances per site
    families.site <- BB22.full %>%
      dplyr::group_by(family.agg, site, bbspecies) %>%
      dplyr::summarise(cum.abund = sum(Abundance))
    families.site$family.agg <- as.factor(families.site$family.agg)
    
# 3. plot families per site
    # color palette for families
    palette.fams=c("#375E97", "#80BD9E", "#FA812F", "#F34A4A", "#07575B", "#66A5AD", "#C4DFE6", "#FAAF08", 
                   "#336B87", "#5D535E", "#DFE166", "#1995AD", "#258039", "#73605B", "#4897D8", "#DDBC95",
                   "#A3A599", "#A1D6E2", "#88A550", "#75B1A9", "#D9B44A", "#4F6457", "#ACD0C0", "#0F1B07", 
                   "#F7EFE2", "#D09683", "#A1BE95", "#F62A00", "#20948B", "#9B4F0F", "#CB0000")
    # filled bar plot per site with abundance of families
    setwd(output)
    ggplot(families.site, aes(fill=family.agg, y=cum.abund, x=site)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Plant families per Site") +
      labs(fill='Plant families') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.fams, limits = unique(BB22.full$family.agg))
    ggsave(paste("./proportion in pollen/PlantFamilies_per_Site.png", sep = ""), width = 16, height = 8, device = "png")
    
# 4. produce data frame with families' abundances per region
    families.region <- BB22.full %>%
      dplyr::group_by(family.agg, region, bbspecies) %>%
      dplyr::summarise(cum.abund = sum(Abundance))
    families.region$family.agg <- as.factor(families.region$family.agg)
    
# 5. plot families per region    
    ggplot(families.region, aes(fill=family.agg, y=cum.abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Plant families per region") +
      labs(fill='Plant families') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.fams, limits = unique(BB22.full$family.agg))
    ggsave(paste("./proportion in pollen/PlantFamilies_per_Region.png", sep = ""), width = 16, height = 8, device = "png", )
    setwd(input) 

#### plot plant growth form per site, region and species
# 1. produce data frame with exotic/native per site
    ex.nat.site <- BB22.full %>%
      dplyr::group_by(native_exotic, site, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    ex.nat.site$native_exotic <- as.factor(ex.nat.site$native_exotic)
    
# 2. plot native/exotic per site
    # create color palette
    palette.ex =c("#518A45", "#BDA509")
    
    setwd(output)
    ggplot(ex.nat.site, aes(fill=native_exotic, y=abund, x=site)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Origin status per site") + ylab("Proportion of species in the pollen")+
      labs(fill='Origin status') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'))
    ggsave(paste("./proportion in pollen/OriginStatus_per_Site.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with exotic/native per region
    ex.nat.region <- BB22.full %>%
      dplyr::group_by(native_exotic, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    ex.nat.region$native_exotic <- as.factor(ex.nat.region$native_exotic)
    
# 4. plot native/exotic per region    
    ggplot(ex.nat.region, aes(fill=native_exotic, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Origin status per region") + ylab("Proportion of species in the pollen")+
      labs(fill='Origin status') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'))
    ggsave(paste("./proportion in pollen/OriginStatus_per_Region.png", sep = ""), width = 16, height = 8)
    setwd(input)
    

#### growth type per site and species
# 1. produce data frame with growth form per site
    growth.site <- BB22.full %>%
      dplyr::group_by(growth_form_category, site, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    growth.site$growth_form_category <- as.factor(growth.site$growth_form_category)
    
# 2. plot growth form per site
    # create color palette
    palette.growth =c("#375E97", "#518A45", "#FA812F", "#F34A4A", "#5A99AD", "#A3A3A3")
    
    setwd(output)
    ggplot(growth.site, aes(fill=growth_form_category, y=abund, x=site)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Growth Form") + ylab("Proportion of species in the pollen")+ 
      labs(fill='Growth form') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.growth)
    ggsave(paste("./proportion in pollen/GrowthForm_per_Site.png", sep = ""), width = 16, height = 8)
    
# 3. produce data frame with growth form per region
    growth.region <- BB22.full %>%
      dplyr::group_by(growth_form_category, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    growth.region$growth_form_category <- as.factor(growth.region$growth_form_category)
    
# 4. plot growth form per region    
    ggplot(growth.region, aes(fill=growth_form_category, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Growth Form") + ylab("Proportion of species in the pollen")+ 
      labs(fill='Growth form') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.growth)
    ggsave(paste("./proportion in pollen/GrowthForm_per_Region.png", sep = ""), width = 16, height = 8)
    setwd(input)
    
#### blossom class per site and species
# 1. produce data frame with blossom class per site
    blossom.site <- BB22.full %>%
      dplyr::group_by(structural_blossom_class, site, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    blossom.site$structural_blossom_class <- as.factor(blossom.site$structural_blossom_class)

# 2. plot blossom class per site
    # create color palette
    palette.bloss=c("#375E97", "#80BD9E", "#FA812F", "#758A30", "#07575B", "#D95F56", "#C4DFE6")
    
    setwd(output)
    ggplot(blossom.site, aes(fill=structural_blossom_class, y=abund, x=site)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Blossom Class per site") + ylab("Proportion of species in the pollen")+
      labs(fill='Blossom Class') +
      theme(axis.text.x = element_text(angle = 90))+
      scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"))
    ggsave(paste("./proportion in pollen/./proportion in pollen/BlossomClass_per_Site.png", sep = ""), width = 16, height = 8)
    
# 3. produce data frame with blossom class per region
    blossom.region <- BB22.full %>%
      dplyr::group_by(structural_blossom_class, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    blossom.region$structural_blossom_class <- as.factor(blossom.region$structural_blossom_class)
    
# 4. plot blossom class per region
    ggplot(blossom.region, aes(fill=structural_blossom_class, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Blossom Class per site") + ylab("Proportion of species in the pollen")+
      labs(fill='Blossom Class') +
      theme(axis.text.x = element_text(angle = 90))+
      scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"))
    ggsave(paste("./proportion in pollen/./proportion in pollen/./proportion in pollen/BlossomClass_per_Region.png", sep = ""), width = 16, height = 8)
    
    setwd(input)

    
############################
#### PHYLOGENETIC TREE #####    

# load required library
library(V.PhyloMaker)
library(ape)

# make a list of all occurring species
# species <- BB22.full$plant.species

# create a data frame with taxa (species, genus, family)
phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

# run the phylo-function (load data since it takes a long time)
tree.result <- phylo.maker(phylo, scenarios=c("S1","S2","S3"))

# plot the phylogenies with node ages displayed with sceanario 3
tree <- plot.phylo(tree.result$scenario.3, cex = 0.5, main = "Phylogenetic tree of species in pollen")
tree
# write.tree(result$scenario.3, "tree.tre")

# bubble plot relative abundances

# add: on region level
# add:on landscape level
BB22.full.site <- BB22.full%>%
  dplyr::group_by(site, plant.species, bbspecies)%>%
  dplyr::summarise(Abundance = sum(Abundance))
p2 <- ggplot(BB22.full.site, aes(x = site, y =BB22.full.site$plant.species, color = site)) + 
 geom_point(aes(size = Abundance, fill = site, alpha=0.5)) + facet_wrap(~bbspecies) +theme_classic()
