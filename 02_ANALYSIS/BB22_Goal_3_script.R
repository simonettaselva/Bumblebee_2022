################################################
# GOAL 3: Diet Consistency
# by Simonetta Selva
#
# Created: January 18, 2023
# Project: Bumblebee 2022
################################################

# Aim: Examine the diet consistency between and within urban and rural landscapes for both species.


# xxx ----
## preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

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

## B.lapidarius ----
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.list.bb <- list()

for (i in sitenames) {
  site.list.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
}

# import species lists form GBIF and InfoFlora (file BB22_sites_plants_GBIF.R)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# create list with intersection of species lists for all regions
list.intersections <- list()
bb.intersections <- list()
table.interaction <- list()
list.table <- list()
correlation.tables <- list()

for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    list.intersections[[paste(i,"/",j, sep="")]] <- intersect(site.list.occ[[i]], site.list.occ[[j]])
    
    # filter binary dataframe for the two sites
    list.table[[paste(i,"/",j, sep="")]] <- select(BB22.full.table, i, j)
    # bb.intersections[[paste("bee",i,":",i,"/",j, sep="")]] <- intersect(site.list.bb[[i]], list.intersections[[paste(i,"/",j, sep="")]])
    
    # filter binary data frame for species also found in pollen per site
    table.interaction[[paste("bee",i,":",i,"/",j, sep="")]] <- 
      list.table[[paste(i,"/",j, sep="")]] %>% filter(row.names(list.table[[paste(i,"/",j, sep="")]]) %in% site.list.bb[[i]])
    
  } # end loop j
} # end loop i

cor.test(table.interaction$`beeZHUA:ZHUA/ZHUC`$ZHUA, table.interaction$`beeZHUA:ZHUA/ZHUC`$ZHUC, method=c("pearson"))

typeof(table.interaction$`beeZHUA:ZHUA/ZHUC`$ZHUA)


