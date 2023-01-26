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
# convert into binary data frame for comaring and later calculate correlations between sites
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.species.bb <- list()

for (i in sitenames) {
  site.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  j <- 0
  site.data.InfoFlora$species <- c()
  for (j in 1:nrow(site.data.InfoFlora)){
    # species names are different than in GBIF -> need to be adapted
    site.data.InfoFlora$species[j] <- paste(unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[1],
                                                 unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[2])
  } # end loop j
  
  # store species per site in a list
  site.list.InfoFlora[[i]] <- unique(site.data.InfoFlora$species)
  
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif[[i]] <- unique(site.data.gbif$species)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora[[i]])
  y <- as.factor(site.list.gbif[[i]])
  site.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

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
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
library(stringr)
#prepare list to store species lists
site.list.InfoFlora <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(genus = word(Taxon, 1),
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora[[i]] <- unique(site.data.InfoFlora$genus)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(genus = as_factor(genus)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif[[i]] <- unique(site.data.gbif$genus)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.genus.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora[[i]])
  y <- as.factor(site.list.gbif[[i]])
  site.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all sites
intersection.sites.genus <- list()
table.sites.genus <- list()
correlation.site.genus <- data.frame(genus = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common genus of two sites
    intersection.sites.genus[[paste(i,"/",j, sep="")]] <- intersect(site.genus.occ[[i]], site.genus.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.sites.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.genus <- rbind(correlation.site.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.site.genus$site1 <- as.factor(correlation.site.genus$site1)
matrix.site.genus <- split(correlation.site.genus, f = correlation.site.genus$site1)
temp <- list.cbind(matrix.site.genus)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.genus$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.genus <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.genus) <- "numeric"

# store in list
pasc.site[["genus"]] <- matrix.site.genus


### family ----

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.family.bb <- list()

for (i in sitenames) {
  site.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(family = Familie,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora[[i]] <- unique(site.data.InfoFlora$family)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(family = as_factor(family)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif[[i]] <- unique(site.data.gbif$family)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.family.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora[[i]])
  y <- as.factor(site.list.gbif[[i]])
  site.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all sites
intersection.sites.family <- list()
table.sites.family <- list()
correlation.site.family <- data.frame(family = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common family of two sites
    intersection.sites.family[[paste(i,"/",j, sep="")]] <- intersect(site.family.occ[[i]], site.family.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.sites.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.family <- rbind(correlation.site.family, correlation)
    
  } # end loop j
} # end loop i

correlation.site.family$site1 <- as.factor(correlation.site.family$site1)
matrix.site.family <- split(correlation.site.family, f = correlation.site.family$site1)
temp <- list.cbind(matrix.site.family)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.family$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.family <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.family) <- "numeric"

# store in list
pasc.site[["family"]] <- matrix.site.family

### plotting ----
library(ggcorrplot)
a <- ggcorrplot(pasc.site[["species"]], hc.order = F, type = "lower", #method = "circle",
                outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(-0.02,1), low = "white", high =  "#07575B")

b <- ggcorrplot(pasc.site[["genus"]], hc.order = F, type = "lower", #method = "circle",
                outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(-0.02,1), low = "white", high =  "#07575B")

c <- ggcorrplot(pasc.site[["family"]], hc.order = F, type = "lower", #method = "circle",
                outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(-0.02,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
plot1 <- ggarrange(a, b, c, ncol = 1, nrow = 3,
          labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("B.pascuroum", face = "bold", size = 22))
# 
ggsave("./04_Goal_3/corr_pasc_site.png", width = 10, height = 30)
setwd(input)

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



