################################################
# R site plant occurances GBIF and InfoFlora
# by Simonetta Selva
#
# Created: November 19, 2022
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
setwd(input)

# my data: create species lists for each site and save in list()
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

BB22.abund <- read_csv("BB22.abund.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species)) 
BB22.abund <- BB22.abund%>%
  mutate(site = fct_recode(site,"ZHUA" = "ZH_A","ZHUB" = "ZH_B","ZHUC" = "ZH_C","ZHRD" = "ZH_D","ZHRE" = "ZH_E","ZHRF" = "ZH_F","BEUA" = "BE_A",
                           "BEUB" = "BE_B","BEUC" = "BE_C","BERD" = "BE_D","BSUA" = "BS_A","BSUB" = "BS_B","BSUC" = "BS_C","BSRD" = "BS_D",
                           "BSRE" = "BS_E","BSRF" = "BS_F"))
site.list <- list()

for (i in sitenames) {
  site.list[[i]]  <- unique(BB22.abund$species[BB22.abund$site == i]) %>% 
    droplevels()
}

#### GBIF ####

# GBIF: create species lists for each site and save in list()
site.list.gbif <- list()

for (i in sitenames) {
  site.data.gbif  <- read_csv(paste("./sites_plant_list/BB22_", i, ".csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif [[i]] <- unique(site.data.gbif$species)
}


# how many species of total occurred species are visited by BB per site
ratios <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.gbif [[i]]))
  gbif.occured <- length(site.list.gbif [[i]])
  first <- c(i, shared/gbif.occured, "visited by bumblebees")
  second<- c(i, (gbif.occured-shared)/gbif.occured, "GBIF")
  ratios <- rbind(ratios, first, second)
  colnames(ratios)<- c("site", "percentage", "group")
  ratios <- as.data.frame(ratios)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("Comp_GBIF_visited_per_site.jpeg")


#### INFO FLORA ####

setwd(input)
site.list.InfoFlora <- list()
# site.data.InfoFlora <- read_csv(paste("./sites_plant_list/01_InfoFlora/BB22_ZHUA_InfoFlora.csv", sep = ""))%>% 
# mutate(species = as_factor(Species)) %>% 
#   filter(occ != 0)
# rename(site.data.InfoFlora, occ = ?)

for (i in sitenames) {
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/01_InfoFlora/BB22_", i, "_InfoFlora.csv", sep = "")) %>%
    mutate(species = as_factor(Species)) %>% 
    filter(occ != 0)
  j <- 0
  site.data.InfoFlora$species <- c()
  for (j in 1:nrow(site.data.InfoFlora)){
    site.data.InfoFlora$species[j] <- paste(unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[1],
                                            unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[2])
  }
  site.list.InfoFlora [[i]] <- unique(site.data.InfoFlora$species)
}

# how many species of total occurred species are visited by BB per site
ratios.InfoFlora <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.InfoFlora [[i]]))
  InfoFlora.occured <- length(site.list.InfoFlora [[i]])
  first <- c(i, shared/InfoFlora.occured, "visited by bumblebees")
  # hoi <- rbind(ratios, first)
  second<- c(i, (InfoFlora.occured-shared)/InfoFlora.occured, "InfoFlora")
  ratios.InfoFlora <- rbind(ratios.InfoFlora, first, second)
  colnames(ratios.InfoFlora)<- c("site", "percentage", "group")
  ratios.InfoFlora <- as.data.frame(ratios.InfoFlora)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.InfoFlora, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("Comp_InfoFlora_visited_per_site.jpeg")



#### Intersection of InfoFlora and GBIF per site ####
ratios.comparison <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list.gbif[[i]], site.list.InfoFlora [[i]]))
  comparison.occured <- length(site.list.InfoFlora [[i]])
  first <- c(i, shared/comparison.occured, "shared by GBIF")
  # hoi <- rbind(ratios, first)
  second<- c(i, (comparison.occured-shared)/comparison.occured, "InfoFlora")
  ratios.comparison <- rbind(ratios.comparison, first, second)
  colnames(ratios.comparison)<- c("site", "percentage", "group")
  ratios.comparison <- as.data.frame(ratios.comparison)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.comparison, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison GBIF and InfoFlora species lists per site") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("Comp_GBIF_InfoFlora.jpeg")
setwd(input)






















