################################################
# Consistency of Diet IN a site
# by Simonetta Selva
#
# Created: January, 31th, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())

# preparation and load data ----

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
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

#remove information on leg and body from ID
BB22.full$ID.short = as.factor(substring(BB22.full$ID,1, nchar(BB22.full$ID)-1)) 

# sites' names for looping
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

# B.PACUORUM ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

# prepare list element to store correlation matrices for species, genus and family
pasc.ID <- list()
pasc.ID.mean <- list()

## calculate distance ----
for (i in sitenames) {
  # select a site
  # i <- "BSRE"
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site%>%
    group_by(ID.short, family)%>%
    summarise(abundance = sum(Abundance))%>% 
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
  
  # ON FAMILY
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var="ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method="bray")
  pasc.ID[[i]] <- as.numeric(dist[1:1000])
  pasc.ID.mean[[i]] <- mean(dist)
}


## plotting ----
# transform data
pasc.ID.df <- do.call(cbind, pasc.ID)
pasc.ID.df <- melt(pasc.ID.df)[, -1]
colnames(pasc.ID.df) <- c("site", "dist")
pasc.ID.df$landscape <- substr(pasc.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

ggplot(pasc.ID.df, aes(x=site, y = dist, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") + ylab("distance") +
  ggtitle("B.pascuorum: distance between individuals per site") +
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

setwd(output)
ggsave("./01_Goal 0/distances/pasc_dist_family.png", width = 10, height = 10)
setwd(input)

# B.LAPIDARIUS ----
# only use B.B.lapidarius
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

# prepare list element to store correlation matrices for species, genus and family
lapi.ID <- list()
lapi.ID.mean <- list()

## calculate distance ----
for (i in sitenames) {
  # select a site
  # i <- "BSRE"
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site%>%
    group_by(ID.short, family)%>%
    summarise(abundance = sum(Abundance))%>% 
    distinct() # remove duplicates
  
  BB22_full.ab.new <- c()
  for (h in unique(BB22_full.ab$ID.short)) {
    # h <- "36_ZHUA_p"
    temp <- BB22_full.ab[BB22_full.ab$ID.short==h,]
    perc <- sum(temp$abundance)
    for (k in 1:nrow(temp)) {
      temp$ab.new[k] <- 100/perc*temp$abundance[k]/100
    }
    BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
  }
  
  # ON FAMILY
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var="ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method="bray")
  lapi.ID[[i]] <- as.numeric(dist[1:1000])
  lapi.ID.mean[[i]] <- mean(dist)
}

## plotting ----
# transform data
lapi.ID.df <- do.call(cbind, lapi.ID)
lapi.ID.df <- melt(lapi.ID.df)[, -1]
colnames(lapi.ID.df) <- c("site", "dist")
lapi.ID.df$landscape <- substr(lapi.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

ggplot(lapi.ID.df, aes(x=site, y = dist, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") + ylab("distance") +
  ggtitle("B.lapidarius: distance between individuals per site") +
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

setwd(output)
ggsave("./01_Goal 0/distances/lapi_dist_family.png", width = 10, height = 10)
setwd(input)
