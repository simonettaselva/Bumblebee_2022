################################################
# R compute Functional Diversity
# by Simonetta Selva
#
# Created: December 14, 2022
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

setwd(input)

# load data
BB22_full <- read_csv("BB22_full.csv")

# choose only to numeric data
BB22_full.numeric <- BB22_full %>% 
  summarise(ID = ID,
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            plant.species = plant.species,
            binom.abund = binom.abund,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)
BB22_full.numeric$Flowering_duration <- as.numeric(BB22_full.numeric$Flowering_duration)

# Impute missing values and reduce data dimensionality using PCA 
require(caret)
require(vegan)
BB22_full.pascuroum <- BB22_full.numeric[BB22_full.numeric$bbspecies == "B.pascuorum",]
BB22_full.pascuroum.species <- BB22_full.pascuroum %>% 
  select(plant.species, Flowering_duration, Flowering_start, growth_form_numeric, 
         structural_blossom_numeric, sugar.concentration, symmetry_numeric, plant_height_m) %>% 
  filter(plant.species!="Fabaceae sp.") %>% # Remove this uninteresting entry, where no traits are found
  distinct() # remove duplicates

trt.mis.pred <- preProcess(as.data.frame(BB22_full.pascuroum.species[,-c(1)]), "knnImpute")
traits <- predict(trt.mis.pred, BB22_full.pascuroum.species[,-c(1)]); head(trt.mis.pred)
rownames(traits)  = BB22_full.pascuroum.species$plant.species

# PCA -> reduce dimensionality
trt.pca <- prcomp(traits, scale. = T, center = T)
cumsum(trt.pca$sdev/sum(trt.pca$sdev))
trt.scaled <- scores(trt.pca)[,1:2] # adjust number of axes for each group

# bring into wide format (for each BB species) on site level
library(reshape2)
wide.pascuroum <- dcast(BB22_full.pascuroum, site ~ plant.species, value.var="binom.abund")
sp.pa.pascuroum <- decostand(wide.pascuroum[,-1], "pa")
rownames(sp.pa.pascuroum)  = wide.pascuroum$site #re-introduce rownames

# baskets_fruits_weights = sp.pa.pascuroum
# fruits_traits = traits
# fruits_traits_cat = traits_cat

# create data frame with information on the variables (numeric, character, factor etc.)
# for details go to "mFD: General Workflow", by Camille Magneville 2022
traits_cat <- data_frame(trait_name = colnames(traits),
                         trait_type = rep("Q", length(colnames(traits))))

# Species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = traits, 
  stop_if_NA = TRUE)






