################################################
# GOAL 3: Diet Consistency
# by Simonetta Selva
#
# Created: January 26, 2023
# Project: Bumblebee 2022
################################################

# Aim: Examine the diet consistency between and within urban and rural landscapes for both species.


# preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rlist)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# set working directory
setwd(input) 

# load data
# import species lists on bee pollen (file BB22_sites_plants_GBIF.R)
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

# B.PACUORUM ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

## site level  ----
# prepare list element to store correlation matrices for species, genus and family
pasc.site <- list()

### species ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.species.bb <- list()

for (i in sitenames) {
  site.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
}

bla  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/", "ZHUC", "_IF_1500.csv", sep = "")) # no genus in this data frame
blad  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", "ZHUC", "_1500.csv", sep = ""))
  
# import data from GBIF and InfoFlora
# InfoFlora
site.list.InfoFlora <- list()
for (i in sitenames) {
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  j <- 0
  site.data.InfoFlora$species <- c()
  for (j in 1:nrow(site.data.InfoFlora)){
    site.data.InfoFlora$species[j] <- paste(unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[1],
                                                 unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[2])
  }
  site.list.InfoFlora[[i]] <- unique(site.data.InfoFlora$species)
}

# GBIF
site.list.gbif <- list()

for (i in sitenames) {
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif[[i]] <- unique(site.data.gbif$species)
}

# combine species list GBIF and InfoFlora (site level)
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora[[i]])
  y <- as.factor(site.list.gbif[[i]])
  site.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
}

# create list with intersection of species lists for all sites
intersection.sites.species <- list()
table.sites.species <- list()
correlation.site.species <- data.frame(species = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    intersection.sites.species[[paste(i,"/",j, sep="")]] <- intersect(site.species.occ[[i]], site.species.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.sites.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.species <- rbind(correlation.site.species, correlation)
    
  } # end loop j
} # end loop i

correlation.site.species$site1 <- as.factor(correlation.site.species$site1)
matrix.site.species <- split(correlation.site.species, f = correlation.site.species$site1)
temp <- list.cbind(matrix.site.species)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.species$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.species <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.species) <- "numeric"

# store in list
pasc.site[["species"]] <- matrix.site.species











### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.genus.bb <- list()

for (i in sitenames) {
  site.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$site == i])
}



### family ----

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.family.bb <- list()

for (i in sitenames) {
  site.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$site == i])
}

## region level ----
### species ----
### genus ----
### family ----

# B.LAPIDARIUS ----
## site level  ----
### species ----
### genus ----
### family ----

## region level ----
### species ----
### genus ----
### family ----



