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
locationnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

BB22.abund <- read_csv("BB22.abund.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species)) 
BB22.abund <- BB22.abund%>%
  mutate(site = fct_recode(site,"ZHUA" = "ZH_A","ZHUB" = "ZH_B","ZHUC" = "ZH_C","ZHRD" = "ZH_D","ZHRE" = "ZH_E","ZHRF" = "ZH_F","BEUA" = "BE_A",
                           "BEUB" = "BE_B","BEUC" = "BE_C","BERD" = "BE_D","BSUA" = "BS_A","BSUB" = "BS_B","BSUC" = "BS_C","BSRD" = "BS_D",
                           "BSRE" = "BS_E","BSRF" = "BS_F")) 
BB22.abund <- BB22.abund%>%
  mutate(locationnames = substr(site, 1, 3))

# create species lists per site
site.list <- list()
for (i in sitenames) {
  site.list[[i]]  <- unique(BB22.abund$species[BB22.abund$site == i]) %>%
    droplevels()
}

# create species lists per location
location.list <- list()
for (i in locationnames) {
  location.list[[i]]  <- unique(BB22.abund$species[BB22.abund$locationnames == i]) %>%
    droplevels()
}



#### GBIF ####

# GBIF: create species lists for each site and save in list()
site.list.gbif <- list()

for (i in sitenames) {
  site.data.gbif  <- read_csv(paste("./sites_plant_list/03_GBIF/BB22_", i, ".csv", sep = "")) %>%
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
# ggsave("Comp_GBIF_visited_per_site.jpeg")


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
# ggsave("Comp_InfoFlora_visited_per_site.jpeg")



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
  ggtitle("Comparison GBIF and InfoFlora 5x5 species lists per site") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("Comp_GBIF_InfoFlora.jpeg")
setwd(input)


#### MORE EXACT INFOFLORA DATA ####
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.list.InfoFlora.ex <- list()

for (i in sitenames) {
  site.data.InfoFlora.ex  <- read_csv(paste("./sites_plant_list/02_InfoFlora_exact/",i , "_IF.csv", sep = "")) 
  site.data.InfoFlora.ex <- site.data.InfoFlora.ex %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora.ex)))
  j <- 0
  site.data.InfoFlora.ex$species <- c()
  for (j in 1:nrow(site.data.InfoFlora.ex)){
    site.data.InfoFlora.ex$species[j] <- paste(unlist(strsplit(site.data.InfoFlora.ex$Species[j], split=' ', fixed=TRUE))[1],
                                            unlist(strsplit(site.data.InfoFlora.ex$Species[j], split=' ', fixed=TRUE))[2])
  }
  site.list.InfoFlora.ex [[i]] <- unique(site.data.InfoFlora.ex$species)
}

#### Intersection of the two InfoFlora data sets per site ####
ratios.comparison <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list.InfoFlora [[i]], site.list.InfoFlora.ex[[i]]))
  comparison.occured <- length(site.list.InfoFlora [[i]])
  first <- c(i, shared/comparison.occured, "shared by more precise species lists")
  # hoi <- rbind(ratios, first)
  second<- c(i, (comparison.occured-shared)/comparison.occured, "InfoFlora 5x5")
  ratios.comparison <- rbind(ratios.comparison, first, second)
  colnames(ratios.comparison)<- c("site", "percentage", "group")
  ratios.comparison <- as.data.frame(ratios.comparison)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.comparison, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison InfoFlora 5x5 and exact species lists per site") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("Comp_both_InfoFlora.jpeg")
setwd(input)


#### Intersection of InfoFlora more exact and GBIF per site ####
ratios.comparison <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list.gbif[[i]], site.list.InfoFlora.ex[[i]]))
  comparison.occured <- length(site.list.InfoFlora.ex[[i]])
  first <- c(i, shared/comparison.occured, "shared by GBIF")
  # hoi <- rbind(ratios, first)
  second<- c(i, (comparison.occured-shared)/comparison.occured, "InfoFlora exact")
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
# ggsave("Comp_GBIF_InfoFlora.jpeg")
setwd(input)


# how many species of total occurred species by InfoFlora exact are visited by BB per site
ratios.InfoFlora <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.InfoFlora.ex [[i]]))
  InfoFlora.occured <- length(site.list.InfoFlora.ex [[i]])
  first <- c(i, shared/InfoFlora.occured, "visited by bumblebees")
  # hoi <- rbind(ratios, first)
  second<- c(i, (InfoFlora.occured-shared)/InfoFlora.occured, "InfoFlora exact")
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
# ggsave("Comp_InfoFloraEx_visited_per_site.jpeg")
setwd(input)



#### TRY EVERYTHING ON LOCATION LEVEL (not on site level) ####
locationnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

location.list.InfoFlora.ex <- list()
# temp <- c()
for (i in locationnames) {
  location.data.InfoFlora.ex  <- read_csv(paste("./sites_plant_list/02_InfoFlora_exact/",i, "_comb_IF.csv", sep = "")) 
  location.data.InfoFlora.ex <- location.data.InfoFlora.ex %>%
    summarise(Species = Taxon,
              location = rep(i , nrow(location.data.InfoFlora.ex)))
  j <- 0
  location.data.InfoFlora.ex$species <- NA
  for (j in 1:nrow(location.data.InfoFlora.ex)){
    location.data.InfoFlora.ex$species[j] <- paste(unlist(strsplit(location.data.InfoFlora.ex$Species[j], split=' ', fixed=TRUE))[1],
                                               unlist(strsplit(location.data.InfoFlora.ex$Species[j], split=' ', fixed=TRUE))[2])
  }
  location.list.InfoFlora.ex [[i]] <- unique(location.data.InfoFlora.ex$species)
  
}


# how many species of total occurred species by InfoFlora exact are visited by BB per location
ratios.InfoFlora <- c()

for (i in locationnames) {
  shared <- length(intersect(location.list[[i]], location.list.InfoFlora.ex [[i]]))
  InfoFlora.occured <- length(location.list.InfoFlora.ex [[i]])
  first <- c(i, shared/InfoFlora.occured, "visited by bumblebees")
  # hoi <- rbind(ratios, first)
  second<- c(i, (InfoFlora.occured-shared)/InfoFlora.occured, "InfoFlora exact")
  ratios.InfoFlora <- rbind(ratios.InfoFlora, first, second)
  colnames(ratios.InfoFlora)<- c("location", "percentage", "group")
  ratios.InfoFlora <- as.data.frame(ratios.InfoFlora)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.InfoFlora, aes(fill=group, y=percentage, x=location)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("Comp_InfoFloraEx_visited_per_location.jpeg")
setwd(input)






