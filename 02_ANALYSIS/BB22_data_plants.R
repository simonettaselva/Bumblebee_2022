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

species.not.covered <- setdiff(BB22.plantspecies, shared.species); species.not.covered
# check names on "plants of the world online"; add the missing and adapt synonyms

JC_floral_traits <- read.csv2("floral_traits_ss.csv",sep = ",")%>% 
  mutate(Sp.merge = as_factor(Sp.merge),
         plant_height_m = as.numeric(plant_height_m))
JC.plantspecies <- unique(JC_floral_traits$Sp.merge)

BB22.plantspecies <- recode(BB22.plantspecies,
                            "Phedimus spurius" = "Sedum spurium",
                            "Phedimus aizoon" = "Sedum aizoon",
                            "Clematis sp. YX-2018" = "Clematis sp.",
                            "Salvia amplexicaulis" = "Salvia nemorosa",
                            "Phedimus kamtschaticus" = "Sedum aizoon",
                            "Tilia americana x Tilia x moltkei" = "Tilia americana",
                            "Hylotelephium telephium" = "Sedum telephium",
                            "Petunia sp. LR-2018" = "Petunia sp.",
                            "Onopordum illyricum" = "Onopordum acanthium",
                            "Fabaceae spc" = "Fabaceae sp.",
                            "Potentilla glabra" = "Dasiphora fruticosa",
                            "Linaria spc" = "Linaria sp.",
                            "Begonia spc" = "Begonia sp.",
                            "Lamiaceae spc" = "Lamiaceae sp.",
                            "Allium spc" = "Allium sp.",
                            "Asteraceae spc" = "Asteraceae sp.",
                            "Boraginaceae spc" = "Boraginaceae sp.",
                            "Betonica officinalis" = "Stachys officinalis",
                            "Cyclamen spc" = "Cyclamen sp.",
                            "Crepis spc" = "Crepis sp.",
                            "Eruca pinnatifida" = "Eruca vesicaria",
                            "Sedum montanum" = "Sempervivum montanum",
                            "Trifolium spc" = "Trifolium sp.",
                            "Lotus spc" = "Lotus sp.",
                            "Hypochaeris spc" = "Hypochaeris sp.",
                            "Nepeta spc" = "Nepeta sp.",
                            "x Chitalpa tashkentensis" = "xChitalpa tashkentensis")

# see how many of hte plant species are shared by data sets
shared.species <- intersect(JC.plantspecies,BB22.plantspecies)
length(shared.species)
# 230 are shared

plants_meta <- JC_floral_traits %>% 
  filter(Sp.merge %in% BB22.plantspecies)%>% 
  rename(plant.species = Sp.merge)%>%
  select(-City)%>%
  droplevels()


# add data on nutrients
plants.nutrients <- read.csv2("phenology_pollen_nectar_sugar_database_copy.csv",sep = ",")
str(plants.nutrients)
plants.nutrients <- plants.nutrients%>% 
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
            pollen.flower = total.pollen.per.flower..the.number.of.grains.)%>% 
  mutate(mean.yearly.air.temperature = as.numeric(mean.yearly.air.temperature),
         total.yearly.pecipitation= as.numeric(total.yearly.pecipitation),
         flowering.start =as.numeric(flowering.start),                                              
         flowering.end = as.numeric(flowering.end),                                              
         flowering.peak = as.numeric(flowering.peak),
         flowering.peak.2..if.occured. = as.numeric(flowering.peak.2..if.occured.),
         flowering.lenght = as.numeric(flowering.lenght),
         flower.longevity = as.numeric(flower.longevity),
         sugar.concentration...in.nectar = as.numeric(sugar.concentration...in.nectar),
         pollen.flower = as.numeric(pollen.flower))

# filter for species in BB22 dataset
plants.nutri.filter <- plants.nutrients %>% 
  filter(plant.species %in% BB22.plantspecies)

# summarize entries
plants.nutri.filter <- plants.nutri.filter %>% 
  group_by(plant.species)%>%
  summarize(plant.species = plant.species,
            flowering.start = mean(flowering.start, na.rm=TRUE),
            flowering.end = mean(flowering.end, na.rm=TRUE),
            flowering.peak = mean(flowering.peak, na.rm=TRUE),
            flowering.lenght = mean(flowering.lenght, na.rm=TRUE),
            flower.longevity = mean(flower.longevity, na.rm=TRUE),
            sugar.concentration = mean(sugar.concentration...in.nectar, na.rm=TRUE),
            pollen.flower = mean(pollen.flower, na.rm=TRUE))%>% 
  distinct()
            
BB22_plant_traits <- merge(plants_meta,plants.nutri.filter, by = "plant.species", all=TRUE)
# write_csv(BB22_plant_traits, "BB22_plant_traits.csv")
# data frame has been changes manually -> do not write again
names(BB22_plant_traits)

BB22_plant_traits <- read.csv("BB22_plant_traits.csv", sep=",")

### look at plant traits ###
ggplot(BB22_plant_traits, aes(x = native_exotic)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x = pollination_mode)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x = growth_form_category)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x = nectar)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x = pollen)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x = oil)) +
  geom_bar()
ggplot(BB22_plant_traits, aes(x=plant_height_m)) + # 172 NAs
  geom_histogram()
ggplot(BB22_plant_traits, aes(x=sugar.concentration)) + # 172 NAs
  geom_histogram()

str(BB22_plant_traits)
