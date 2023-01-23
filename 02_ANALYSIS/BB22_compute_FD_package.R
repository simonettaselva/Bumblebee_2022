################################################
# R compute Functional Diversity with FD package
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
BB22_full <- read_csv("BB22_full.csv")%>% 
  mutate(ID = as.character(ID))
BB22_full$ID.short = as.factor(substring(BB22_full$ID,1, nchar(BB22_full$ID)-1)) #remove information on leg and body from ID

# choose only to numeric data
BB22_full.numeric <- BB22_full %>% 
  summarise(ID.short = ID.short, #remove information on leg and body from ID
            ID = as.factor(ID),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            plant.species = plant.species,
            binom.abund = binom.abund,
            abundance = Abundance,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)
BB22_full.numeric$Flowering_duration <- as.numeric(BB22_full.numeric$Flowering_duration)

# correlation analysis of the variables
library(Hmisc)
library(corrplot)
res <- rcorr(as.matrix(BB22_full.numeric[,c(12:18)]),type="pearson")
M <- cor(BB22_full.numeric[,c(12:18)], use = "complete.obs")
corrplot::corrplot(M, type="upper", order="hclust", p.mat = res$P, sig.level = 0.05)

# remove plant height, flowering start and symmetry
res2 <- rcorr(as.matrix(BB22_full.numeric[,c(12, 14, 15, 16)]), type="pearson") # all p values are <0.05
M <- cor(BB22_full.numeric[,c(12, 14, 15, 16)], use = "complete.obs")
corrplot::corrplot(M, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.05)

# remove plant height, flowering start and symmetry from data set
BB22_full.red <- BB22_full.numeric%>% 
  select(-Flowering_start,
         -symmetry_numeric,
         -plant_height_m) 

# for (i in levels(BB22_full.red$bbspecies)) { #loop trough bumblebee species
#   for (j in c("ID.short", "site", "landscape")) {

j <- "ID.short"
i <- "B.pascuorum"


# 1. Impute missing values and reduce data dimensionality using PCA 
require(caret)
require(vegan)

# remove plant species entrences that do not have traits
BB22_full.loop <- BB22_full.red[BB22_full.red$bbspecies == i,]%>%
  filter(plant.species!= "Fabaceae sp.",
         plant.species!= "Cyclamen sp.",
         plant.species!= "Petunia sp.",
         plant.species!= "Mandevilla atroviolacea" )%>%
  mutate(Flowering_duration = as.numeric(Flowering_duration))%>% 
  droplevels()

# select columnes used in computing FDs
BB22_full.loop.species <- BB22_full.loop %>% 
  select(plant.species, Flowering_duration, structural_blossom_numeric, 
         sugar.concentration, growth_form_numeric) %>% 
  distinct() # remove duplicates

# impute missing data
trt.mis.pred <- preProcess(as.data.frame(BB22_full.loop.species[,-c(1)]), "knnImpute")
traits <- predict(trt.mis.pred, BB22_full.loop.species[,-c(1)]); head(trt.mis.pred)
traits <- as.data.frame(traits)
rownames(traits)  <- BB22_full.loop.species$plant.species
traits <- traits[order(row.names(traits)), ] # reorder traits into alphabetical order


# bring plants species per site/ID into wide format (for each BB species)
library(reshape2)
# presence/absence data
wide <- dcast(BB22_full.loop, BB22_full.loop[[j]] ~ plant.species, value.var="binom.abund")[,-1]
rownames(wide)  <- levels(BB22_full.loop[[j]])
sp.pa <- decostand(wide, "pa")
sp.pa <- as.matrix(sp.pa) #turn into matrix

### ACHTUNG KANN NICHT IM LOOP LAUFEN !!!
# relative abundance data
BB22_full.ab <- BB22_full.loop%>%
  group_by(ID.short, plant.species)%>%
  summarise(abundance = sum(abundance))%>% 
  distinct() # remove duplicates

BB22_full.ab.new <- c()
for (h in unique(BB22_full.ab$ID.short)) {
  temp <- BB22_full.ab[BB22_full.ab$ID.short==h,]
  perc <- sum(temp$abundance)
  for (k in 1:nrow(temp)) {
    temp$ab.new[k] <- 100/perc*temp$abundance[k]/100
  }
  BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
}

wide <- dcast(BB22_full.ab.new, ID.short ~ plant.species, value.var="ab.new")[,-1]
rownames(wide)  <- levels(BB22_full.ab.new$ID.short)
# abundance
sp.ab <- as.matrix(wide)
# presence/absence
sp.pa <- decostand(wide, "pa")
sp.pa <- as.matrix(sp.pa) #turn into matrix

  
# compute FD
library(FD)
fd <- FD::dbFD(x = traits , a = sp.ab)





