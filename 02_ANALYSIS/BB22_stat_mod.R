################################################
# R data preparation PLANTS first looks
# by Simonetta Selva
#
# Created: October 19, 2022
# Project: Bumblebee 2022
################################################

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
