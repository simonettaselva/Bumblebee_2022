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

##########################
### POLLEN COMPOSITION ###
##########################

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
  B22.full$family.agg <- BB22.full$family
  # group all rare families into "other families"
  for(h in rare.families$family){
    for(i in 1:nrow(BB22.full)){
      if(BB22.full$family.agg[i] == h){
        BB22.full$family.agg[i] <- "Other families"
      }
    }
  }

# 2. produce data frame with famlilie's abundances per site
    families.site <- BB22.full %>%
      dplyr::group_by(family.agg, site) %>%
      dplyr::summarise(cum.abund = sum(Abundance))
  






