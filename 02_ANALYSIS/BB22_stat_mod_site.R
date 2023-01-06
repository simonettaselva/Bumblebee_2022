################################################
# Statistical Models script SITE
# by Simonetta Selva
#
# Created: January, 6th, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22.metrics.body <- read_csv("BB22.metrics.csv") %>%
  filter(bborgan == "B")%>%
  mutate(ID = substring(ID,1, nchar(ID)-1))
BB22.bb.traits <- read_csv("BB22.bb.traits.csv")

# replace the names for urban and rural
for (i in 1:nrow(BB22.bb.traits)) {
  if (BB22.bb.traits$landscape[i] == "urban") {
    BB22.bb.traits$landscape[i] = "U"
  } else {
    BB22.bb.traits$landscape[i] = "R"
  }
}
# rename column site to match other dataframes
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(site = paste(location, landscape, replicate, sep = ""))

# 
BB22.bb.traits.site <- BB22.bb.traits %>%
  group_by(site) %>%
  summarise(site = site,
            intertegular_distance = mean(intertegular_distance),
            glossa = mean(glossa),
            prementum = mean(prementum),
            proboscis_length = mean(proboscis_length),
            proboscis_ratio = mean(proboscis_ratio),
            fore_wing_length = mean(fore_wing_length),
            fore_wing_ratio = mean(fore_wing_ratio),
            corbicula_length = mean(corbicula_length),
            corbicula_ratio = mean(corbicula_ratio))

BB22.metrics.traits <- merge(BB22.metrics.body, BB22.bb.traits[, -c(2:8)], by = "ID")
BB22.metrics.traits$site <- as.factor(BB22.metrics.traits$site)




# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")%>%
  mutate(site = as_factor(site))

# add site coordinates to the data frame (in LV95)
BB22.metrics.traits <- merge(BB22.metrics.traits,BB22.sites.meta[, c(1,3,4)], by  = "site", all.x=TRUE) 



