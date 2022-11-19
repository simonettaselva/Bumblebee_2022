################################################
# R site plant occurances GBIF
# by Simonetta Selva
#
# Created: Nove,ber 19, 2022
# Project: Bumblebee 2022
################################################
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

## set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

for (i in sitenames) {
  site.data <- read_csv(paste("./sites_plant_list/BB22_", i, ".csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  assign(paste(i, ".list", sep = ""), unique(site.data$species)) 
  site.list[[i]] <- unique(site.data$species)
}






