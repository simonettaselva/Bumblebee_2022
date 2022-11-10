################################################
# R data preparation PLANTS first looks
# by Simonetta Selva
#
# Created: October 19, 2022
# Project: Bumblebee 2022
################################################
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)


# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"


# load data
setwd(input)
BB22 <- read.csv("BB22_data_tres_0.01.csv")
BB22 <- BB22%>% 
  mutate(binom.abund = if_else(Abundance == 0, 0, 1),
         ID = substring(Sample, 2),
         site = paste(location, replicate, sep="_"),
         bbspecies = as_factor(bbspecies),
         OTU = as_factor(OTU))
levels(BB22$bbspecies) <- c("B.lapidarius", "B.pascuorum")

BB22.plantspecies <- unique(BB22$OTU)

#import Joans dataset
JC_floral_traits <- read.csv2("floral_traits.csv",sep = ",")%>% 
  mutate(Sp.merge = as_factor(Sp.merge))
JC.plantspecies <- unique(JC_floral_traits$Sp.merge)

# see how many of hte plant species are shared by datasets
shared.species <- intersect(JC.plantspecies,BB22.plantspecies)
# 169 are shared

species.not.covered <- setdiff(BB22.plantspecies, shared.species)


#neu

###############
plants <- read.csv2("phenology_pollen_nectar_sugar_database_copy.csv",sep = ",")
plants.sum <- plants%>% 
  summarize(plant.species = plant.species,
            community = community,
            continent.region = continent.region,
            country = country,
            location.name = location.name,
            year = year.of.the..data.collection,
            mean.yearly.air.temperature = mean.yearly.air.temperature.in.degrees.Celsius,
            total.yearly.pecipitation= total.yearly..pecipitation.in.mm,
            flowering.start =flowering.start,                                                  
            flowering.end = flowering.end,                                                    
            flowering.peak = flowering.peak,
            flowering.peak.2..if.occured. = flowering.peak.2..if.occured.,
            flowering.lenght = flowering.lenght,
            flower.longevity = flower.longevity..days.,
            sugar.concentration...in.nectar = sugar.concentration...in.nectar,
            pollen.flower = total.pollen.per.flower..the.number.of.grains.)
names(plants)






