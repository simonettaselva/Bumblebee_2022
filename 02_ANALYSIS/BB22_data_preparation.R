################################################
# R data preparation
# by Simonetta Selva
#
# Created: November 17, 2022
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

#create binom. abundance variable
BB22 <- BB22%>% 
  mutate(binom.abund = if_else(Abundance == 0, 0, 1),
         ID = substring(Sample, 2),
         site = paste(location, landscape, replicate, sep=""),
         bbspecies = as_factor(bbspecies))%>% 
  select(-Sample, -X, -project, -xID)
levels(BB22$bbspecies) <- c("B.lapidarius", "B.pascuorum") # rename BB species
BB22 <- BB22[, c(15,16,1,3,4,5,6,7,8,9,10,11,12,13,2,14)] # reorder columns
BB22$OTU <- recode(BB22$OTU, # match plant names
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


# 230 species
# 47 famillies

#remove entries with abundance = 0
setwd(input)
BB22.abund <- BB22%>% 
  filter(binom.abund == 1)

# add information on plants species (file BB22_data_plants.R)
BB22_plant_traits <- read_csv("BB22_plant_traits_added.csv")
BB22 <- BB22 %>% rename(plant.species = OTU)

BB22_full <- full_join(BB22, BB22_plant_traits, by = "plant.species")%>%
  select(-Family, -Genus, -species)%>% 
  rename(species = Sp2) %>% 
  filter(binom.abund == 1)%>% 
  mutate(native_exotic = as_factor(native_exotic),
         pollination_mode = as_factor(pollination_mode),
         growth_form_category = as_factor(growth_form_category),
         inflorescence = as_factor(inflorescence),
         structural_blossom_class = as_factor(structural_blossom_class),
         symmetry = as_factor(symmetry),
         Flowering_months_duration = as.numeric(Flowering_months_duration))
BB22_full <- BB22_full[, c(1,2,4,5,6,7,8,9,10,11,12,13,16,3,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30)] # reorder columns

# create numerical values for factors
BB22_full$growth_form_numeric <- BB22_full$growth_form_category
levels(BB22_full$growth_form_numeric) <- 1:5
BB22_full$growth_form_numeric <- as.numeric(BB22_full$growth_form_numeric)
# herb = 1, shrub = 2, climber = 3, tree = 4 , shrub/tree = 5

BB22_full$structural_blossom_numeric <- BB22_full$structural_blossom_class
levels(BB22_full$structural_blossom_numeric) <- 1:7
BB22_full$structural_blossom_numeric <- as.numeric(BB22_full$structural_blossom_numeric)
# flag = 1, gullet = 2, dish_bowl = 3, stalk_disk = 4, bell_trumpet = 5, tube = 6, brush = 7

BB22_full$symmetry_numeric <- BB22_full$symmetry
levels(BB22_full$symmetry_numeric) <- 1:3
BB22_full$symmetry_numeric <- as.numeric(BB22_full$symmetry_numeric)
# sigomorph = 1, actinomorph = 2, no_symmetry = 3

# write.csv(BB22_full, "BB22_full.csv")

# see how many NA in sugar concentration
sum(is.na(BB22_full$sugar.concentration)) #around 1/5 of entries


# summarize data COMMUNITY METRICS

library(codyn)
library(vegan)
library(fundiversity)

## CWM
BB22.metrics <- BB22_full %>% 
  group_by(ID) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            Shannon = diversity(Abundance),
            NrSpecies=n_distinct(species),
            Flowering_duration_cwm = weighted.mean(Flowering_months_duration, Abundance, na.rm=T),
            Flowering_start_cwm = weighted.mean(start_flowering, Abundance, na.rm=T),
            growth_form_cwm = weighted.mean(growth_form_numeric, Abundance, na.rm=T),
            structural_blossom_cwm = weighted.mean(structural_blossom_numeric, Abundance, na.rm=T),
            sugar_concentration_cwm = weighted.mean(sugar.concentration, Abundance, na.rm=T),
            plant_height_cwm = weighted.mean(plant_height_m, Abundance, na.rm=T)
            ) %>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F"),
         landscape = fct_relevel(landscape, "U", "R")) %>%
  distinct()

setwd(input)
# write_csv(BB22.metrics, "BB22.metrics.csv")

# look at distributions of the CWM variables
CWM <- colnames(BB22.metrics)[str_detect(colnames(BB22.metrics), "_cwm", negate = FALSE)]
par(mfrow = c(2, 3))
for (i in CWM){
  hist(BB22.metrics[[i]], main = i)
}
par(mfrow = c(1,1))


### Hahs and Fournier et al. 2022 - Compute functional diversity indices ------------
### Script to calculate various functional diversity metrics
### load useful function ("Toolkit") -------------------------------------------
BB22_full.numeric <- BB22_full %>% 
  summarise(ID = ID,
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            plant.species = plant.species,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)

source("Bertrand_Function_ToolKit.R")

# Imput missing values and reduce data dimensionality using PCA 
# -> enable to calculate FD metrics for sites with low diversity
require(caret)
require(vegan)
BB22_full.mis.model = preProcess(BB22_full.numeric[,8:13], "knnImpute")
# preProcess not possible since to much missing data -> either delete rows aor use different approach
BB22_full.mis.model = predict(BB22_full.mis.model, BB22_full.numeric[,8:13]); head(BB22_full.mis.model)
#SS: does not work

### imputing Simos style (not correct probably)
impute <- function(x){
  x[is.na(x)] <- mean (x, na.rm = TRUE)
  return(x)
}
df = NULL
for (column in colnames(BB22_full.numeric[,9:15])) {
  imputed <- impute(BB22_full.numeric[[column]])
  df <- cbind(df, imputed)
}
colnames(df) <- colnames(BB22_full.numeric[,9:15])
# trt.mis.pred <- cbind(BB22_full.numeric[, 1:8], df)
trt.mis.pred <- df

# PCA -> reduce dimensionality
trt.pca <- prcomp(trt.mis.pred, scale. = T, center = T)
cumsum(trt.pca$sdev/sum(trt.pca$sdev))
trt.scaled <- scores(trt.pca)[,1:2] # adjust number of axes for each group

### FD incices from mFD package of Magneville et al. 2022 ------------------------
### mFD
library(mFD)
FD_mFD <- calc_mFD(trt=trt.mis.pred)







