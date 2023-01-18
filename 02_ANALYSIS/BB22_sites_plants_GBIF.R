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
region <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

# BB22.abund <- read_csv("BB22.abund.csv") %>%
#   mutate(site = as_factor(site),
#          species = as_factor(species)) 
# BB22.abund <- BB22.abund%>%
#   mutate(site = paste(location, landscape, replicate, sep = ""))
# BB22.abund <- BB22.abund%>%
#   mutate(region = substr(site, 1, 3))

BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

# create species lists per site
site.list <- list()
for (i in sitenames) {
  site.list[[i]]  <- unique(BB22.abund$species[BB22.abund$site == i]) %>%
    droplevels()
}

saveRDS(site.list, file="site_list_bb.RData")

# create species lists per region
location.list <- list()
for (i in region) {
  location.list[[i]]  <- unique(BB22.abund$species[BB22.abund$region == i]) %>%
    droplevels()
}

saveRDS(site.list, file="region_list_bb.RData")


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
# ggsave("./GBIF and InfoFlora/Comp_GBIF_visited_per_site.jpeg")


#### INFO FLORA ####

setwd(input)
site.list.InfoFlora <- list()


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
# ggsave("./GBIF and InfoFlora/Comp_InfoFlora_visited_per_site.jpeg")



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
# ggsave("./GBIF and InfoFlora/Comp_GBIF_InfoFlora.jpeg")
setwd(input)


#### MORE EXACT INFOFLORA DATA ####
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.list.InfoFlora.ex <- list()

for (i in sitenames) {
  site.data.InfoFlora.ex  <- read_csv(paste("./sites_plant_list/03_InfoFlora_800/",i , "_IF.csv", sep = "")) 
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
# ggsave("./GBIF and InfoFlora/Comp_both_InfoFlora.jpeg")
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
# ggsave("./GBIF and InfoFlora/Comp_GBIF_InfoFlora.jpeg")
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
# ggsave("./GBIF and InfoFlora/Comp_InfoFloraEx_visited_per_site.jpeg")
setwd(input)



#### TRY EVERYTHING ON REGION LEVEL (not on site level) ####
region <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

location.list.InfoFlora.ex <- list()
# temp <- c()
for (i in region) {
  location.data.InfoFlora.ex  <- read_csv(paste("./sites_plant_list/03_InfoFlora_800/",i, "_comb_IF.csv", sep = "")) 
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

for (i in region) {
  shared <- length(intersect(location.list[[i]], location.list.InfoFlora.ex [[i]]))
  InfoFlora.occured <- length(location.list.InfoFlora.ex [[i]])
  first <- c(i, shared/InfoFlora.occured, "visited by bumblebees")
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
# ggsave("./GBIF and InfoFlora/Comp_InfoFloraEx_visited_per_location.jpeg")
setwd(input)

# see also if there are plants in pollen not listed in data bases
# on site level INFOFLORA
p <- data.frame()
plot <- data.frame()
for (i in sitenames) {
  intersection <- length(intersect(site.list[[i]], site.list.InfoFlora.ex [[i]]))
  notshared.1 <- length(site.list[[i]][!(site.list[[i]] %in% intersection)])
  notshared.2 <- length(site.list.InfoFlora.ex[[i]][!(site.list.InfoFlora.ex[[i]] %in% intersection)])
  p <- data.frame(value = c(notshared.1, intersection, notshared.2),
                  shared = factor(c("not shared BB", "shared", "not shared InfoFlora"), levels = c("not shared BB", "shared", "not shared InfoFlora")),
                  site = rep(i, 3))
  plot <- rbind(plot,p)
}

setwd(output)
ggplot(plot, aes(fill=shared, y=value, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("./GBIF and InfoFlora/Comp_InfoFloraEx_visited_per_site_2.jpeg")
setwd(input)

# on site level GBIF
p <- data.frame()
plot <- data.frame()
for (i in sitenames) {
  intersection <- length(intersect(site.list[[i]], site.list.gbif [[i]]))
  notshared.1 <- length(site.list[[i]][!(site.list[[i]] %in% intersection)])
  notshared.2 <- length(site.list.gbif[[i]][!(site.list.gbif[[i]] %in% intersection)])
  p <- data.frame(value = c(notshared.1, intersection, notshared.2),
                  shared = factor(c("not shared BB", "shared", "not shared GBIF"), levels = c("not shared BB", "shared", "not shared GBIF")),
                  site = rep(i, 3))
  plot <- rbind(plot,p)
}

setwd(output)
ggplot(plot, aes(fill=shared, y=value, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("./GBIF and InfoFlora/Comp_GBIF_visited_per_site_2.jpeg")
setwd(input)


# on location level
p <- data.frame()
plot <- data.frame()
for (i in region) {
  intersection <- length(intersect(location.list[[i]], location.list.InfoFlora.ex [[i]]))
  notshared.1 <- length(location.list[[i]][!(location.list[[i]] %in% intersection)])
  notshared.2 <- length(location.list.InfoFlora.ex[[i]][!(location.list.InfoFlora.ex[[i]] %in% intersection)])
  p <- data.frame(value = c(notshared.1, intersection, notshared.2),
                     shared = factor(c("not shared BB", "shared", "not shared InfoFlora"), levels = c("not shared BB", "shared", "not shared InfoFlora")),
                     site = rep(i, 3))
  plot <- rbind(plot,p)
}

setwd(output)
ggplot(plot, aes(fill=shared, y=value, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("location") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("./GBIF and InfoFlora/Comp_InfoFloraEx_visited_per_location_2.jpeg")
setwd(input)


## see which plant species are shared
intersection <- list()

for (i in sitenames) {
  intersection[[i]] <- intersect(site.list[[i]], site.list.InfoFlora.ex [[i]])
}

intersection$ZHRF
intersection$ZHRD

#### 1500m Buffer ####
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
region <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

# InfoFlora
site.list.InfoFlora.1500 <- list()
for (i in sitenames) {
  site.data.InfoFlora.1500  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  site.data.InfoFlora.1500 <- site.data.InfoFlora.1500 %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora.1500)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora.1500)))
  j <- 0
  site.data.InfoFlora.1500$species <- c()
  for (j in 1:nrow(site.data.InfoFlora.1500)){
    site.data.InfoFlora.1500$species[j] <- paste(unlist(strsplit(site.data.InfoFlora.1500$Species[j], split=' ', fixed=TRUE))[1],
                                               unlist(strsplit(site.data.InfoFlora.1500$Species[j], split=' ', fixed=TRUE))[2])
  }
  site.list.InfoFlora.1500[[i]] <- unique(site.data.InfoFlora.1500$species)
}


# GBIF
site.list.gbif.1500 <- list()

for (i in sitenames) {
  site.data.gbif.1500  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.1500[[i]] <- unique(site.data.gbif.1500$species)
}

# combine species list GBIF and InfoFlora (site level)
site.list.1500 <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.1500 [[i]])
  y <- as.factor(site.list.gbif.1500 [[i]])
  site.list.1500[[i]] <- unique(c(x,y))%>% 
    droplevels()
}
saveRDS(site.list.1500, file="sp_list_gbif_infoflora.RData")


# combine species list GBIF and InfoFlora (region level)
site.list.1500 <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.1500 [[i]])
  y <- as.factor(site.list.gbif.1500 [[i]])
  site.list.1500[[i]] <- unique(c(x,y))%>% 
    droplevels()
}
saveRDS(site.list.1500, file="sp_list_gbif_infoflora.RData")



ratios.1500 <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.1500[[i]]))
  occured <- length(site.list.1500[[i]])
  first <- c(i, shared/occured, "percentage of plants visited by bumblebees")
  second<- c(i, (occured-shared)/occured, "GBIF and InfoFlora (1500m)")
  ratios.1500 <- rbind(ratios.1500, first, second)
  colnames(ratios.1500)<- c("site", "percentage", "group")
  ratios.1500 <- as.data.frame(ratios.1500)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.1500, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("./GBIF and InfoFlora/Comp_1500_visited_per_site.jpeg", width = 8, height = 8)
setwd(input)


#### 2500m Buffer ####
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
# InfoFlora
site.list.InfoFlora.2500 <- list()
for (i in sitenames) {
  site.data.InfoFlora.2500  <- read_csv(paste("./sites_plant_list/05_InfoFlora_2500/",i , "_IF_2500.csv", sep = "")) 
  site.data.InfoFlora.2500 <- site.data.InfoFlora.2500 %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora.2500)))
  j <- 0
  site.data.InfoFlora.2500$species <- c()
  for (j in 1:nrow(site.data.InfoFlora.2500)){
    site.data.InfoFlora.2500$species[j] <- paste(unlist(strsplit(site.data.InfoFlora.2500$Species[j], split=' ', fixed=TRUE))[1],
                                                 unlist(strsplit(site.data.InfoFlora.2500$Species[j], split=' ', fixed=TRUE))[2])
  }
  site.list.InfoFlora.2500[[i]] <- unique(site.data.InfoFlora.2500$species)
}

# GBIF
site.list.gbif.2500 <- list()

for (i in sitenames) {
  site.data.gbif.2500  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_2500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.2500[[i]] <- unique(site.data.gbif.2500$species)
}


# combine species list GBIF and InfoFlora
site.list.2500 <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.2500 [[i]])
  y <- as.factor(site.list.gbif.2500 [[i]])
  site.list.2500[[i]] <- unique(c(x,y))%>% 
    droplevels()
}

ratios.2500 <- c()

for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.2500[[i]]))
  occured <- length(site.list.2500[[i]])
  first <- c(i, shared/occured, "percentage of plants visited by bumblebees")
  second<- c(i, (occured-shared)/occured, "GBIF and InfoFlora (2500m)")
  ratios.2500 <- rbind(ratios.2500, first, second)
  colnames(ratios.2500)<- c("site", "percentage", "group")
  ratios.2500 <- as.data.frame(ratios.2500)%>% 
    mutate(percentage = as.numeric(percentage))
}

setwd(output)
ggplot(ratios.2500, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
# ggsave("./GBIF and InfoFlora/Comp_2500_visited_per_site.jpeg", width = 8, height = 8)
setwd(input)















