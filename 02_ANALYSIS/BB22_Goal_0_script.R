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
BB22.full <- read_csv("BB22_full.csv") # import data set 

# initialize region as column
BB22.full$region <- c() 

# add site and region as columns
  for (i in 1:nrow(BB22.full)) {
    BB22.full$site[i] <-paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep="")
    BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep="")
  }

# divide data set into leg and body pollen
BB22.full.leg <- BB22.full[BB22.full$bborgan=="L",]
BB22.full.body <- BB22.full[BB22.full$bborgan=="B",]

#### LEG POLLEN ####
#### plot plant families per site, region and species

# 1. find the 30 most abundant families
  families.overview <- BB22.full.leg%>%
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
  BB22.full.leg$family.agg <- BB22.full.leg$family
  # group all rare families into "other families"
  for(h in rare.families$family){
    for(i in 1:nrow(BB22.full.leg)){
      if(BB22.full.leg$family.agg[i] == h){
        BB22.full.leg$family.agg[i] <- "Other families"
      }
    }
  }

# 2. produce data frame with families' abundances per site with leg pollen
    families.site <- BB22.full.leg %>%
      dplyr::group_by(family.agg, site, bbspecies) %>%
      dplyr::summarise(cum.abund = sum(Abundance))
    families.site$family.agg <- as.factor(families.site$family.agg)
    
# 3. plot families per site
    # color palette for families
    palette.fams=c("#07575B", "#5D535E", "#C4DFE6", "#FAAF08", "#336B87" , "#DFE166", "#1995AD", 
                  "#4897D8", "#80BD9E", "#FA812F", "#66A5AD", "#F34A4A", "#258039", "#375E97", "#73605B", "#DDBC95",
                   "#A3A599", "#A1D6E2", "#88A550", "#75B1A9", "#D9B44A", "#4F6457", "#ACD0C0", "#0F1B07", 
                   "#F7EFE2", "#D09683", "#A1BE95", "#F62A00", "#20948B", "#9B4F0F", "#CB0000")
    # filled bar plot per site with abundance of families
    setwd(output)
    ggplot(families.site, aes(fill=family.agg, y=cum.abund, x=site)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Plant families per site LEG") +
      labs(fill='Plant families') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.fams, name = "Plant families")
    ggsave(paste("./proportion in pollen/PlantFamilies_per_Site_leg.png", sep = ""), width = 16, height = 8, device = "png")
    
# 4. produce data frame with families' abundances per region
    families.region <- BB22.full.leg %>%
      dplyr::group_by(family.agg, region, bbspecies) %>%
      dplyr::summarise(cum.abund = sum(Abundance))
    families.region$family.agg <- as.factor(families.region$family.agg)
    
# 5. plot families per region
    ggplot(families.region, aes(fill=family.agg, y=cum.abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Plant families per region LEG") +
      labs(fill='Plant families') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.fams, name = "Plant families")
    ggsave(paste("./proportion in pollen/PlantFamilies_per_Region_leg.png", sep = ""), width = 16, height = 8, device = "png", )
    setwd(input) 

#### plot plant growth form per site, region and species
# 1. produce data frame with exotic/native per site
    ex.nat.site <- BB22.full.leg %>%
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
      ggtitle("Origin status per site LEG") + ylab("Proportion of species in the pollen")+
      labs(fill='Origin status') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name = "")
    ggsave(paste("./proportion in pollen/OriginStatus_per_Site_leg.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with exotic/native per region
    ex.nat.region <- BB22.full.leg %>%
      dplyr::group_by(native_exotic, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    ex.nat.region$native_exotic <- as.factor(ex.nat.region$native_exotic)
    
# 4. plot native/exotic per region    
    ggplot(ex.nat.region, aes(fill=native_exotic, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Origin status per region LEG") + ylab("Proportion of species in the pollen")+
      labs(fill='Origin status') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name="")
    ggsave(paste("./proportion in pollen/OriginStatus_per_Region_leg.png", sep = ""), width = 16, height = 8)
    setwd(input)
    

#### growth type per site and species
# 1. produce data frame with growth form per site
    growth.site <- BB22.full.leg %>%
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
      ggtitle("Growth Form LEG") + ylab("Proportion of species in the pollen")+ 
      labs(fill='Growth form') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.growth, name="")
    ggsave(paste("./proportion in pollen/GrowthForm_per_Site_leg.png", sep = ""), width = 16, height = 8)
    
# 3. produce data frame with growth form per region
    growth.region <- BB22.full.leg %>%
      dplyr::group_by(growth_form_category, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    growth.region$growth_form_category <- as.factor(growth.region$growth_form_category)
    
# 4. plot growth form per region    
    ggplot(growth.region, aes(fill=growth_form_category, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Growth Form LEG") + ylab("Proportion of species in the pollen")+ 
      labs(fill='Growth form') +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values=palette.growth, name="")
    ggsave(paste("./proportion in pollen/GrowthForm_per_Region_leg.png", sep = ""), width = 16, height = 8)
    setwd(input)
    
#### blossom class per site and species
# 1. produce data frame with blossom class per site
    blossom.site <- BB22.full.leg %>%
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
      ggtitle("Blossom Class per site LEG") + ylab("Proportion of species in the pollen")+
      labs(fill='Blossom Class') +
      theme(axis.text.x = element_text(angle = 90))+
      scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"),
                        name = "")
    ggsave(paste("./proportion in pollen/BlossomClass_per_Site_leg.png", sep = ""), width = 16, height = 8)
    
# 3. produce data frame with blossom class per region
    blossom.region <- BB22.full.leg %>%
      dplyr::group_by(structural_blossom_class, region, bbspecies) %>%
      dplyr::summarise(abund = sum(binom.abund))
    blossom.region$structural_blossom_class <- as.factor(blossom.region$structural_blossom_class)
    
# 4. plot blossom class per region
    ggplot(blossom.region, aes(fill=structural_blossom_class, y=abund, x=region)) + 
      geom_bar(position="fill", stat="identity")+ 
      theme_classic(base_size = 20) +
      facet_wrap(~bbspecies)+ 
      ggtitle("Blossom Class per site LEG") + ylab("Proportion of species in the pollen")+
      labs(fill='Blossom Class') +
      theme(axis.text.x = element_text(angle = 90))+
      scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"),
                        name = "")
    ggsave(paste("./proportion in pollen/BlossomClass_per_Region_leg.png", sep = ""), width = 16, height = 8)
    setwd(input)

    
#### BODY POLLEN ####
#### plot plant families per site, region and species

# 1. find the 30 most abundant families
families.overview <- BB22.full.body%>%
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
BB22.full.body$family.agg <- BB22.full.body$family
# group all rare families into "other families"
for(h in rare.families$family){
  for(i in 1:nrow(BB22.full.body)){
    if(BB22.full.body$family.agg[i] == h){
      BB22.full.body$family.agg[i] <- "Other families"
    }
  }
}

# 2. produce data frame with families' abundances per site with body pollen
families.site <- BB22.full.body %>%
  dplyr::group_by(family.agg, site, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.site$family.agg <- as.factor(families.site$family.agg)

# 3. plot families per site
# color palette for families
palette.fams=c("#07575B", "#5D535E", "#C4DFE6", "#DFE166","#336B87","#FAAF08", "#1995AD", 
                        "#4897D8", "#80BD9E", "#FA812F", "#66A5AD", "#F34A4A", "#258039", "#375E97", "#73605B", "#DDBC95",
                        "#A3A599", "#A1D6E2", "#88A550", "#75B1A9", "#D9B44A", "#4F6457", "#ACD0C0", "#0F1B07", 
                        "#F7EFE2", "#D09683", "#A1BE95", "#F62A00", "#20948B", "#9B4F0F", "#CB0000")
# filled bar plot per site with abundance of families
setwd(output)
ggplot(families.site, aes(fill=family.agg, y=cum.abund, x=site)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Plant families per site BODY") +
  labs(fill='Plant families') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.fams, name = "Plant families")
ggsave(paste("./proportion in pollen/PlantFamilies_per_Site_body.png", sep = ""), width = 16, height = 8, device = "png")

# 4. produce data frame with families' abundances per region
families.region <- BB22.full.body %>%
  dplyr::group_by(family.agg, region, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.region$family.agg <- as.factor(families.region$family.agg)

# 5. plot families per region
ggplot(families.region, aes(fill=family.agg, y=cum.abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Plant families per region BODY") +
  labs(fill='Plant families') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.fams, name = "Plant families")
ggsave(paste("./proportion in pollen/PlantFamilies_per_Region_body.png", sep = ""), width = 16, height = 8, device = "png", )
setwd(input) 

#### plot plant growth form per site, region and species
# 1. produce data frame with exotic/native per site
ex.nat.site <- BB22.full.body %>%
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
  ggtitle("Origin status per site BODY") + ylab("Proportion of species in the pollen")+
  labs(fill='Origin status') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name = "")
ggsave(paste("./proportion in pollen/OriginStatus_per_Site_body.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with exotic/native per region
ex.nat.region <- BB22.full.body %>%
  dplyr::group_by(native_exotic, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
ex.nat.region$native_exotic <- as.factor(ex.nat.region$native_exotic)

# 4. plot native/exotic per region    
ggplot(ex.nat.region, aes(fill=native_exotic, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Origin status per region BODY") + ylab("Proportion of species in the pollen")+
  labs(fill='Origin status') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.ex, labels=c('Exotic', 'Native'), name="")
ggsave(paste("./proportion in pollen/OriginStatus_per_Region_body.png", sep = ""), width = 16, height = 8)
setwd(input)


#### growth type per site and species
# 1. produce data frame with growth form per site
growth.site <- BB22.full.body %>%
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
  ggtitle("Growth Form BODY") + ylab("Proportion of species in the pollen")+ 
  labs(fill='Growth form') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.growth, name="")
ggsave(paste("./proportion in pollen/GrowthForm_per_Site_body.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with growth form per region
growth.region <- BB22.full.body %>%
  dplyr::group_by(growth_form_category, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
growth.region$growth_form_category <- as.factor(growth.region$growth_form_category)

# 4. plot growth form per region    
ggplot(growth.region, aes(fill=growth_form_category, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Growth Form BODY") + ylab("Proportion of species in the pollen")+ 
  labs(fill='Growth form') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=palette.growth, name="")
ggsave(paste("./proportion in pollen/GrowthForm_per_Region_body.png", sep = ""), width = 16, height = 8)
setwd(input)

#### blossom class per site and species
# 1. produce data frame with blossom class per site
blossom.site <- BB22.full.body %>%
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
  ggtitle("Blossom Class per site BODY") + ylab("Proportion of species in the pollen")+
  labs(fill='Blossom Class') +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"),
                    name = "")
ggsave(paste("./proportion in pollen/BlossomClass_per_Site_body.png", sep = ""), width = 16, height = 8)

# 3. produce data frame with blossom class per region
blossom.region <- BB22.full.body %>%
  dplyr::group_by(structural_blossom_class, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
blossom.region$structural_blossom_class <- as.factor(blossom.region$structural_blossom_class)

# 4. plot blossom class per region
ggplot(blossom.region, aes(fill=structural_blossom_class, y=abund, x=region)) + 
  geom_bar(position="fill", stat="identity")+ 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies)+ 
  ggtitle("Blossom Class per site BODY") + ylab("Proportion of species in the pollen")+
  labs(fill='Blossom Class') +
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values=palette.bloss, labels=c('Bell Trumpet', 'Brush', "Dish Bowl", "Flag", "Gullet", "Stalk Disk", "Tube"),
                    name = "")
ggsave(paste("./proportion in pollen/BlossomClass_per_Region_body.png", sep = ""), width = 16, height = 8)
setwd(input)
    
    
    
############################
#### PHYLOGENETIC TREE #####    

#reset environment
rm(list=ls())

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load required library
library(V.PhyloMaker)
library(ape)

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

# make a list of all occurring species
# species <- BB22.full.leg$plant.species

# only using leg pollen as it reflects better the bees choice of plant species
BB22.full.leg <- BB22.full[BB22.full$bborgan=="L",]

# create a data frame with taxa (species, genus, family)
phylo <- data.frame(species = BB22.full.leg$plant.species, genus = BB22.full.leg$genus, family = BB22.full.leg$family)

# run the phylo-function
tree.result <- phylo.maker(phylo, scenarios=c("S1","S2","S3"))

# plot the phylogenies with node ages displayed with sceanario 3
tree <- plot.phylo(tree.result$scenario.3, cex = 0.5, main = "Phylogenetic tree of species in pollen"); tree
 
# get order of species in tree
phylo.order <- data.frame(sps=tree.result$scenario.3$tip.label)
phylo.order$order <- seq(1, length(phylo.order$sps))
colnames(phylo.order) <- c("plant.species", " order")
phylo.order$plant.species <- sub("_", " ", phylo.order$plant.species)

# create new data frame with needed variables
BB22.full.leg.bubble <- BB22.full.leg%>%
  dplyr::group_by(site, plant.species, bbspecies)%>%
  dplyr::summarise(Abundance = sum(Abundance),
                   region = region,
                   landscape = landscape)

# merge two data frames and order along plant species in phylo-tree
BB22.full.leg.bubble_ordered <- merge(x = BB22.full.leg.bubble, y = phylo.order, by.x = "plant.species")%>%
  mutate(plant.species = as_factor(plant.species))
BB22.full.leg.bubble_ordered <- BB22.full.leg.bubble_ordered[order(BB22.full.leg.bubble_ordered$` order`),]
BB22.full.leg.bubble_ordered$plant.species <- factor(BB22.full.leg.bubble_ordered$plant.species, 
                                                     levels = unique(BB22.full.leg.bubble_ordered$plant.species[order(BB22.full.leg.bubble_ordered$` order`)]))

# plot bubble plot with relative abundances along order of species in tree
setwd(output)
  # site level
  library(pals)
  palette.site <- kelly(18)[3:18] #create color palette for sites
  ggplot(BB22.full.leg.bubble_ordered, aes(x = site, y =BB22.full.leg.bubble_ordered$plant.species, color = site)) + 
   geom_point(aes(size = Abundance, fill = site, alpha=0.5)) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species") +
    theme(axis.text.x = element_text(angle = 90))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_color_manual(values = palette.site, guide = "none")+ #no legend
    scale_fill_manual(values = palette.site, guide = "none") #no legend
  ggsave(paste("./proportion in pollen/Phylo_Bubble_Site_leg.png", sep = ""), width = 16, height = 16)
  
  # region level
  palette.region <- kelly(18)[10:16] #create color palette for region
  ggplot(BB22.full.leg.bubble_ordered, aes(x = region, y =BB22.full.leg.bubble_ordered$plant.species, color = region)) + 
    geom_point(aes(size = Abundance, fill = region, alpha=0.5)) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species")+    
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_color_manual(values = palette.region, guide = "none")+ #no legend
    scale_fill_manual(values = palette.region, guide = "none") #no legend
  ggsave(paste("./proportion in pollen/Phylo_Bubble_Region_leg.png", sep = ""), width = 16, height = 16)
  
  # landscape level
  palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
  ggplot(BB22.full.leg.bubble_ordered, aes(x = landscape, y =BB22.full.leg.bubble_ordered$plant.species, color = landscape)) + 
    geom_point(aes(size = Abundance, fill = landscape, alpha=0.5)) + 
    facet_wrap(~bbspecies) +
    labs(y = "plant species")+ 
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_color_manual(values = palette.landscape, guide = "none") +
    scale_fill_manual(values = palette.landscape, guide = "none")
  ggsave(paste("./proportion in pollen/Phylo_Bubble_Landscape_leg.png", sep = ""), width = 16, height = 16)
  setwd(input)
  

############################
#### LEG/BODY POLLEN DIFFERENCES ####
  
  #reset environment
  rm(list=ls())
  
  # set working directory to main repository
  input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
  output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"
  
  # load data
  setwd(input) # set working directory
  BB22.full <- read_csv("BB22_full.csv") # import data set with
  
  # filter for entries only on body pollen per bumblebee species
  # for B.lapidarius
  BB22.lapi.body <- BB22.full %>% 
    filter(str_detect(bborgan, "B") & bbspecies == "B.lapidarius") 
  # for B.pascuorum
  BB22.pasc.body <- BB22.full %>% 
    filter(str_detect(bborgan, "B") & bbspecies == "B.pascuorum")
  
  # create data frame with cumulative abundance per plant species per bumblebee species
  # for B.lapidarius
  BB22.lapi.body <- BB22.lapi.body %>% 
    group_by(plant.species) %>%
    summarise(cum.abund = sum(Abundance))
  # for B.pascuorum
  BB22.pasc.body <- BB22.pasc.body %>% 
    group_by(plant.species) %>%
    summarise(cum.abund = sum(Abundance))
  
  # filter for entries only on leg pollen per bumblebee species
  # for B.lapidarius
  BB22.lapi.leg <- BB22.full %>% 
    filter(str_detect(bborgan, "L") & bbspecies == "B.lapidarius") 
  # for B.pascuorum
  BB22.pasc.leg <- BB22.full %>% 
    filter(str_detect(bborgan, "L") & bbspecies == "B.pascuorum")
  
  # create data frame with cumulative abundance per plant species per bumblebee species
  # for B.lapidarius
  BB22.lapi.leg <- BB22.lapi.leg %>% 
    group_by(plant.species) %>%
    summarise(cum.abund = sum(Abundance))
  # for B.pascuorum
  BB22.pasc.leg <- BB22.pasc.leg %>% 
    group_by(plant.species) %>%
    summarise(cum.abund = sum(Abundance))
  
  # find common visited plant species body
  intersection.pasc <- length(intersect(BB22.pasc.leg$plant.species, BB22.pasc.body$plant.species))
  notshared.1.pasc <- length(BB22.pasc.leg$plant.species[!(BB22.pasc.leg$plant.species %in% intersection.pasc)])
  notshared.2.pasc <- length(BB22.pasc.body$plant.species[!(BB22.pasc.body$plant.species %in% intersection.pasc)])
  
  intersection.lapi <- length(intersect(BB22.lapi.leg$plant.species, BB22.lapi.body$plant.species))
  notshared.1.lapi <- length(BB22.lapi.leg$plant.species[!(BB22.lapi.leg$plant.species %in% intersection.lapi)])
  notshared.2.lapi <- length(BB22.lapi.body$plant.species[!(BB22.lapi.body$plant.species %in% intersection.lapi)])

  # summarize in data frame for plotting
  p_1 <- data.frame(value = c(notshared.1.pasc, intersection.pasc, notshared.2.pasc,
                            notshared.1.lapi, intersection.lapi, notshared.2.lapi),
                  shared = rep(factor(c("only on corbicula", "shared", "only on body"), 
                                  levels = c("only on corbicula", "shared", "only on body")),2),
                  bbspecies = c(rep("B.pascuorum", 3), rep("B.lapidarius", 3)))
  # plotting
  palette.p1 <- c("#FA812F", "#A1BE95", "#88A550")
  p1 <- ggplot(p_1, aes(fill=shared, y=value, x=bbspecies)) + 
    geom_bar(position="fill", stat="identity") + 
    xlab("Bumblebee species") + ylab("percent")+
    theme_classic(base_size=20) +
    guides(fill=guide_legend(title=""))+
    scale_fill_manual(values=palette.p1, name = ""); p1

  # sum of plant species visited: compute and summarize in data frame 
  p_2 <- data.frame(sum = c(length(BB22.lapi.body$plant.species), length(BB22.lapi.leg$plant.species),
                            length(BB22.pasc.body$plant.species), length(BB22.pasc.leg$plant.species)),
                    bborgan = rep(factor(c("body", "corbicula"), levels = c("body", "corbicula")),2),
                    bbspecies = c(rep("B.pascuorum", 2), rep("B.lapidarius", 2)))
  
  palette.p2 <- c("#88A550","#FA812F")
  p2 <- ggplot(p_2, aes(x=bbspecies, y=sum, fill = bborgan)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab("Sum of Plant Species detected") + xlab("Bumblebee species")+
    theme_classic(base_size=20)+
    scale_fill_manual(values=palette.p2, name = ""); p2
  
  setwd(output)
  plot1 <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B")); plot1
  annotate_figure(plot1, top = text_grob("Comparison of corbicula and body pollen", 
                                         face = "bold", size = 22))
  # ggsave("./proportion in pollen/Comparison_Body_Leg_Pollen.png", width = 16, height = 8)
  setwd(input)

############################
#### LEG/BODY POLLEN DIFFERENCES IN LANDSCAPE ####
  #reset environment
  rm(list=ls())
  
  # set working directory to main repository
  input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
  output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"
  
  # load data
  setwd(input) # set working directory
  BB22.full <- read_csv("BB22_full.csv") # import data set with
  
  # define vectors for loops
  landscape <- c("U", "R")
  organ <- c("B", "L")
  species <- c("B.lapidarius", "B.pascuorum")
  bb <- c()
  # produce data frame per landscape, per species, and per body part
  for (i in landscape) {
    BB22.full.loop <- BB22.full %>%
      filter(landscape == i) #filter for landscape
    for (j in organ) {
      BB22.full.loop.2 <- BB22.full.loop %>%
        filter(str_detect(bborgan, j)) #filter for body part
      for (k in species) {
        BB22.full.loop.3 <- BB22.full.loop.2 %>%
          filter(bbspecies == k) #filter for bumblebee species
        df <- BB22.full.loop.3 %>%
          group_by(plant.species, landscape) %>%
          summarise(cum.abund = sum(Abundance)) # calculate cumulative abundance per plant species per data frame
        if (k == "B.lapidarius") { #for naming the data frames more concise replace name with abbreviation
          bb <- "lapi"
        } else {
          bb <- "pasc"
        }
        assign(paste("BB22", bb, j, i, sep = "."), df) # return named data frame
      } #end loop k
    } #end loop j
  }#end loop i
  
  # find common visited plant species body per landscape
  # urban
  intersection.pasc <- length(intersect(BB22.pasc.L.U$plant.species, BB22.pasc.B.U$plant.species))
  notshared.1.pasc <- length(BB22.pasc.L.U$plant.species[!(BB22.pasc.L.U$plant.species %in% intersection.pasc)])
  notshared.2.pasc <- length(BB22.pasc.B.U$plant.species[!(BB22.pasc.B.U$plant.species %in% intersection.pasc)])
  
  intersection.lapi <- length(intersect(BB22.lapi.L.U$plant.species, BB22.lapi.B.U$plant.species))
  notshared.1.lapi <- length(BB22.lapi.L.U$plant.species[!(BB22.lapi.L.U$plant.species %in% intersection.lapi)])
  notshared.2.lapi <- length(BB22.lapi.B.U$plant.species[!(BB22.lapi.B.U$plant.species %in% intersection.lapi)])
  
  # summarize in data frame for plotting
  df.ratio.U <- data.frame(percent = c(notshared.1.pasc, intersection.pasc, notshared.2.pasc,
                              notshared.1.lapi, intersection.lapi, notshared.2.lapi),
                           plant_species = rep(factor(c("only on corbicula", "shared", "only on body"), 
                                        levels = c("only on corbicula", "shared", "only on body")),2),
                           bbspecies = c(rep("B.pascuorum", 3), rep("B.lapidarius", 3)),
                           landscape = rep("urban",6))
  
  # rural
  intersection.pasc <- length(intersect(BB22.pasc.L.R$plant.species, BB22.pasc.B.R$plant.species))
  notshared.1.pasc <- length(BB22.pasc.L.R$plant.species[!(BB22.pasc.L.R$plant.species %in% intersection.pasc)])
  notshared.2.pasc <- length(BB22.pasc.B.R$plant.species[!(BB22.pasc.B.R$plant.species %in% intersection.pasc)])
  
  intersection.lapi <- length(intersect(BB22.lapi.L.R$plant.species, BB22.lapi.B.R$plant.species))
  notshared.1.lapi <- length(BB22.lapi.L.R$plant.species[!(BB22.lapi.L.R$plant.species %in% intersection.lapi)])
  notshared.2.lapi <- length(BB22.lapi.B.R$plant.species[!(BB22.lapi.B.R$plant.species %in% intersection.lapi)])
  
  # summarize in data frame for plotting
  df.ratio.R <- data.frame(percent = c(notshared.1.pasc, intersection.pasc, notshared.2.pasc,
                                       notshared.1.lapi, intersection.lapi, notshared.2.lapi),
                           plant_species = rep(factor(c("only on corbicula", "shared", "only on body"), 
                                                      levels = c("only on corbicula", "shared", "only on body")),2),
                           bbspecies = c(rep("B.pascuorum", 3), rep("B.lapidarius", 3)),
                           landscape = rep("rural",6))
  
  df.ratio <- rbind(df.ratio.U, df.ratio.R)
  
  # plotting
  palette.p1 <- c("#FA812F", "#A1BE95", "#88A550")
  p3 <- ggplot(df.ratio, aes(fill=plant_species, y=percent, x=landscape)) + 
    geom_bar(position="fill", stat="identity") + xlab("") +
    facet_wrap(~bbspecies)+ 
    theme_classic(base_size=20)+
    scale_fill_manual(values=palette.p1, name = ""); p3
  
  
  # sum of plant species visited per landscape: compute and summarize in data frame 
  df.sum.urban <- data.frame(sum = c(length(BB22.pasc.L.U$plant.species), length(BB22.pasc.B.U$plant.species), 
                                     length(BB22.lapi.L.U$plant.species), length(BB22.lapi.B.U$plant.species)),
                             bbspecies = c(rep("B.pascuorum", 2), rep("B.lapidarius", 2)),
                             bborgan = rep(factor(c("body", "corbicula"), levels = c("body", "corbicula")),2),
                             landscape = rep("urban",4))
  df.sum.rural <- data.frame(sum = c(length(BB22.pasc.L.R$plant.species), length(BB22.pasc.B.R$plant.species), 
                                     length(BB22.lapi.L.R$plant.species), length(BB22.lapi.B.R$plant.species)),
                             bbspecies = c(rep("B.pascuorum", 2), rep("B.lapidarius", 2)),
                             bborgan = rep(factor(c("body", "corbicula"), levels = c("body", "corbicula")),2),
                             landscape = rep("rural",4))
  df.sum <- rbind(df.sum.urban, df.sum.rural)
  
  # plotting
  palette.p2 <- c("#88A550","#FA812F")
  
  p4 <- ggplot(df.sum, aes(x=landscape, y=sum, fill = bborgan)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab("Sum of Plant Species detected") + xlab("")+ 
    facet_wrap(~bbspecies)+ 
    theme_classic(base_size=20)+
    scale_fill_manual(values=palette.p2, name = ""); p4
  
  setwd(output)
  plot1 <- ggarrange(p3, p4, ncol = 2, labels = c("A", "B"))
  annotate_figure(plot1, top = text_grob("Comparison Body and Leg Pollen across Landscapes", 
                                         face = "bold", size = 22))
  ggsave("./proportion in pollen/Comparison_landcape_Body_Leg_Pollen.png", width = 16, height = 8)
  setwd(input)
  
  


