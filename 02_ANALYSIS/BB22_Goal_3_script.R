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
site.list.InfoFlora.species <- list()
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
  site.list.InfoFlora.species[[i]] <- unique(site.data.InfoFlora$species)
  
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.species <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.species[[i]] <- unique(site.data.gbif$species)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.species[[i]])
  y <- as.factor(site.list.gbif.species[[i]])
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
pasc.site[["species"]] <- abs(matrix.site.species)


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
site.list.InfoFlora.genus <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(genus = word(Taxon, 1),
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.genus[[i]] <- unique(site.data.InfoFlora$genus)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.genus <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(genus = as_factor(genus)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.genus[[i]] <- unique(site.data.gbif$genus)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.genus.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.genus[[i]])
  y <- as.factor(site.list.gbif.genus[[i]])
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
pasc.site[["genus"]] <- abs(matrix.site.genus)


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
site.list.InfoFlora.family <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(family = Familie,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.family[[i]] <- unique(site.data.InfoFlora$family)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.family <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(family = as_factor(family)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.family[[i]] <- unique(site.data.gbif$family)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.family.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.family[[i]])
  y <- as.factor(site.list.gbif.family[[i]])
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
pasc.site[["family"]] <- abs(matrix.site.family)

### plotting ----
library(ggcorrplot)
axis.order.site <- c("BERD", 
                     "BSRD", "BSRE", "BSRF",
                     "ZHRD", "ZHRE", "ZHRF",
                     "BEUA", "BEUB", "BEUC",
                     "BSUA", "BSUB", "BSUC", 
                     "ZHUA", "ZHUB", "ZHUC") # to consistently order the axis

a1 <- ggcorrplot(pasc.site[["species"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

a2 <- ggcorrplot(pasc.site[["genus"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

a3 <- ggcorrplot(pasc.site[["family"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
plot1 <- ggarrange(a1, a2, a3, ncol = 1, nrow = 3,
          labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("B.pascuroum", face = "bold", size = 22))
ggsave("./04_Goal_3/corr_pasc_site.png", width = 10, height = 30)
setwd(input)

## region level ----
# prepare list element to store correlation matrices for species, genus and family
pasc.region <- list()

### species ----

# convert into binary data frame for comaring and later calculate correlations between regions
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per region
regionnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")
region.species.bb <- list()

for (i in regionnames) {
  region.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.species$ZHUA, site.list.InfoFlora.species$ZHUB, site.list.InfoFlora.species$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.species$ZHRD, site.list.InfoFlora.species$ZHRE, site.list.InfoFlora.species$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.species$BSUA, site.list.InfoFlora.species$BSUB, site.list.InfoFlora.species$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.species$BSRD, site.list.InfoFlora.species$BSRE, site.list.InfoFlora.species$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.species$BEUA, site.list.InfoFlora.species$BEUB, site.list.InfoFlora.species$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.species$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.species$ZHUA, site.list.gbif.species$ZHUB, site.list.gbif.species$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.species$ZHRD, site.list.gbif.species$ZHRE, site.list.gbif.species$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.species$BSUA, site.list.gbif.species$BSUB, site.list.gbif.species$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.species$BSRD, site.list.gbif.species$BSRE, site.list.gbif.species$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.species$BEUA, site.list.gbif.species$BEUB, site.list.gbif.species$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.species$BERD))


# combine species list GBIF and InfoFlora (region level) in new list
region.species.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

# create list with intersection of species lists for all regions
intersection.regions.species <- list()
table.regions.species <- list()
correlation.region.species <- data.frame(species = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common species of two regions
    intersection.regions.species[[paste(i,"/",j, sep="")]] <- intersect(region.species.occ[[i]], region.species.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.regions.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.species <- rbind(correlation.region.species, correlation)
    
  } # end loop j
} # end loop i

correlation.region.species$region1 <- as.factor(correlation.region.species$region1)
matrix.region.species <- split(correlation.region.species, f = correlation.region.species$region1)
temp <- list.cbind(matrix.region.species)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.species$region1)
rownames(temp) <- regionnames
matrix.region.species <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.species) <- "numeric"

# store in list
pasc.region[["species"]] <- abs(matrix.region.species)


### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per region
region.genus.bb <- list()

for (i in regionnames) {
  region.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.genus$ZHUA, site.list.InfoFlora.genus$ZHUB, site.list.InfoFlora.genus$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.genus$ZHRD, site.list.InfoFlora.genus$ZHRE, site.list.InfoFlora.genus$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.genus$BSUA, site.list.InfoFlora.genus$BSUB, site.list.InfoFlora.genus$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.genus$BSRD, site.list.InfoFlora.genus$BSRE, site.list.InfoFlora.genus$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.genus$BEUA, site.list.InfoFlora.genus$BEUB, site.list.InfoFlora.genus$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.genus$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.genus$ZHUA, site.list.gbif.genus$ZHUB, site.list.gbif.genus$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.genus$ZHRD, site.list.gbif.genus$ZHRE, site.list.gbif.genus$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.genus$BSUA, site.list.gbif.genus$BSUB, site.list.gbif.genus$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.genus$BSRD, site.list.gbif.genus$BSRE, site.list.gbif.genus$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.genus$BEUA, site.list.gbif.genus$BEUB, site.list.gbif.genus$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.genus$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.genus.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all regions
intersection.regions.genus <- list()
table.regions.genus <- list()
correlation.region.genus <- data.frame(genus = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common genus of two regions
    intersection.regions.genus[[paste(i,"/",j, sep="")]] <- intersect(region.genus.occ[[i]], region.genus.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.regions.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.genus <- rbind(correlation.region.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.region.genus$region1 <- as.factor(correlation.region.genus$region1)
matrix.region.genus <- split(correlation.region.genus, f = correlation.region.genus$region1)
temp <- list.cbind(matrix.region.genus)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.genus$region1)
rownames(temp) <- regionnames
matrix.region.genus <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.genus) <- "numeric"

# store in list
pasc.region[["genus"]] <- abs(matrix.region.genus)



### family ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
region.family.bb <- list()

for (i in regionnames) {
  region.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$region == i])
} # end loop i

# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.family$ZHUA, site.list.InfoFlora.family$ZHUB, site.list.InfoFlora.family$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.family$ZHRD, site.list.InfoFlora.family$ZHRE, site.list.InfoFlora.family$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.family$BSUA, site.list.InfoFlora.family$BSUB, site.list.InfoFlora.family$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.family$BSRD, site.list.InfoFlora.family$BSRE, site.list.InfoFlora.family$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.family$BEUA, site.list.InfoFlora.family$BEUB, site.list.InfoFlora.family$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.family$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.family$ZHUA, site.list.gbif.family$ZHUB, site.list.gbif.family$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.family$ZHRD, site.list.gbif.family$ZHRE, site.list.gbif.family$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.family$BSUA, site.list.gbif.family$BSUB, site.list.gbif.family$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.family$BSRD, site.list.gbif.family$BSRE, site.list.gbif.family$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.family$BEUA, site.list.gbif.family$BEUB, site.list.gbif.family$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.family$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.family.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all regions
intersection.regions.family <- list()
table.regions.family <- list()
correlation.region.family <- data.frame(family = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common family of two regions
    intersection.regions.family[[paste(i,"/",j, sep="")]] <- intersect(region.family.occ[[i]], region.family.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.regions.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.family <- rbind(correlation.region.family, correlation)
    
  } # end loop j
} # end loop i

correlation.region.family$region1 <- as.factor(correlation.region.family$region1)
matrix.region.family <- split(correlation.region.family, f = correlation.region.family$region1)
temp <- list.cbind(matrix.region.family)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.family$region1)
rownames(temp) <- regionnames
matrix.region.family <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.family) <- "numeric"

# store in list
pasc.region[["family"]] <- abs(matrix.region.family)

### plotting ----
library(ggcorrplot)
axis.order.reg <- c("BER", "BSR", "ZHR", "BEU", "BSU", "ZHU") # to consistenly order the axis

b1 <- ggcorrplot(pasc.region[["species"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower",
                outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")


b2 <- ggcorrplot(pasc.region[["genus"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

b3 <- ggcorrplot(pasc.region[["family"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
plot2 <- ggarrange(b1,b2,b3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot2, top = text_grob("B.pascuroum", face = "bold", size = 22))
ggsave("./04_Goal_3/corr_pasc_region.png", width = 10, height = 30)
setwd(input)

# B.LAPIDARIUS ----
# only use B.lapidarius
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

## site level  ----
# prepare list element to store correlation matrices for species, genus and family
lapi.site <- list()

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
site.list.InfoFlora.species <- list()
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
  site.list.InfoFlora.species[[i]] <- unique(site.data.InfoFlora$species)
  
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.species <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.species[[i]] <- unique(site.data.gbif$species)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.species[[i]])
  y <- as.factor(site.list.gbif.species[[i]])
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
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
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
lapi.site[["species"]] <- abs(matrix.site.species)

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
site.list.InfoFlora.genus <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(genus = word(Taxon, 1),
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.genus[[i]] <- unique(site.data.InfoFlora$genus)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.genus <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(genus = as_factor(genus)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.genus[[i]] <- unique(site.data.gbif$genus)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.genus.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.genus[[i]])
  y <- as.factor(site.list.gbif.genus[[i]])
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
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
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
lapi.site[["genus"]] <- abs(matrix.site.genus)


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
site.list.InfoFlora.family <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./sites_plant_list/04_InfoFlora_1500/",i , "_IF_1500.csv", sep = "")) 
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(family = Familie,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.family[[i]] <- unique(site.data.InfoFlora$family)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.family <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./sites_plant_list/01_GBIF/BB22_", i, "_1500.csv", sep = "")) %>%
    mutate(family = as_factor(family)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.family[[i]] <- unique(site.data.gbif$family)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.family.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.family[[i]])
  y <- as.factor(site.list.gbif.family[[i]])
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
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
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
lapi.site[["family"]] <- abs(matrix.site.family)

### plotting ----
library(ggcorrplot)
axis.order.site <- c("BERD", 
                     "BSRD", "BSRE", "BSRF",
                     "ZHRD", "ZHRE", "ZHRF",
                     "BEUA", "BEUB", "BEUC",
                     "BSUA", "BSUB", "BSUC", 
                     "ZHUA", "ZHUB", "ZHUC") # to consistenly order the axis

c1 <- ggcorrplot(lapi.site[["species"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

c2 <- ggcorrplot(lapi.site[["genus"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

c3 <- ggcorrplot(lapi.site[["family"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
plot3 <- ggarrange(c1,c2,c3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot3, top = text_grob("B.lapidarius", face = "bold", size = 22))
ggsave("./04_Goal_3/corr_lapi_site.png", width = 10, height = 30)
setwd(input)

## region level ----
# prepare list element to store correlation matrices for species, genus and family
lapi.region <- list()

### species ----

# convert into binary data frame for comaring and later calculate correlations between regions
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per region
regionnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")
region.species.bb <- list()

for (i in regionnames) {
  region.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.species$ZHUA, site.list.InfoFlora.species$ZHUB, site.list.InfoFlora.species$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.species$ZHRD, site.list.InfoFlora.species$ZHRE, site.list.InfoFlora.species$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.species$BSUA, site.list.InfoFlora.species$BSUB, site.list.InfoFlora.species$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.species$BSRD, site.list.InfoFlora.species$BSRE, site.list.InfoFlora.species$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.species$BEUA, site.list.InfoFlora.species$BEUB, site.list.InfoFlora.species$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.species$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.species$ZHUA, site.list.gbif.species$ZHUB, site.list.gbif.species$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.species$ZHRD, site.list.gbif.species$ZHRE, site.list.gbif.species$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.species$BSUA, site.list.gbif.species$BSUB, site.list.gbif.species$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.species$BSRD, site.list.gbif.species$BSRE, site.list.gbif.species$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.species$BEUA, site.list.gbif.species$BEUB, site.list.gbif.species$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.species$BERD))


# combine species list GBIF and InfoFlora (region level) in new list
region.species.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

# create list with intersection of species lists for all regions
intersection.regions.species <- list()
table.regions.species <- list()
correlation.region.species <- data.frame(species = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common species of two regions
    intersection.regions.species[[paste(i,"/",j, sep="")]] <- intersect(region.species.occ[[i]], region.species.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.regions.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.species <- rbind(correlation.region.species, correlation)
    
  } # end loop j
} # end loop i

correlation.region.species$region1 <- as.factor(correlation.region.species$region1)
matrix.region.species <- split(correlation.region.species, f = correlation.region.species$region1)
temp <- list.cbind(matrix.region.species)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.species$region1)
rownames(temp) <- regionnames
matrix.region.species <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.species) <- "numeric"

# store in list
lapi.region[["species"]] <- abs(matrix.region.species)


### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per region
region.genus.bb <- list()

for (i in regionnames) {
  region.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.genus$ZHUA, site.list.InfoFlora.genus$ZHUB, site.list.InfoFlora.genus$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.genus$ZHRD, site.list.InfoFlora.genus$ZHRE, site.list.InfoFlora.genus$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.genus$BSUA, site.list.InfoFlora.genus$BSUB, site.list.InfoFlora.genus$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.genus$BSRD, site.list.InfoFlora.genus$BSRE, site.list.InfoFlora.genus$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.genus$BEUA, site.list.InfoFlora.genus$BEUB, site.list.InfoFlora.genus$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.genus$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.genus$ZHUA, site.list.gbif.genus$ZHUB, site.list.gbif.genus$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.genus$ZHRD, site.list.gbif.genus$ZHRE, site.list.gbif.genus$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.genus$BSUA, site.list.gbif.genus$BSUB, site.list.gbif.genus$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.genus$BSRD, site.list.gbif.genus$BSRE, site.list.gbif.genus$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.genus$BEUA, site.list.gbif.genus$BEUB, site.list.gbif.genus$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.genus$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.genus.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all regions
intersection.regions.genus <- list()
table.regions.genus <- list()
correlation.region.genus <- data.frame(genus = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common genus of two regions
    intersection.regions.genus[[paste(i,"/",j, sep="")]] <- intersect(region.genus.occ[[i]], region.genus.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.regions.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.genus <- rbind(correlation.region.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.region.genus$region1 <- as.factor(correlation.region.genus$region1)
matrix.region.genus <- split(correlation.region.genus, f = correlation.region.genus$region1)
temp <- list.cbind(matrix.region.genus)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.genus$region1)
rownames(temp) <- regionnames
matrix.region.genus <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.genus) <- "numeric"

# store in list
lapi.region[["genus"]] <- abs(matrix.region.genus)

### family ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
region.family.bb <- list()

for (i in regionnames) {
  region.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$region == i])
} # end loop i

# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.family$ZHUA, site.list.InfoFlora.family$ZHUB, site.list.InfoFlora.family$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.family$ZHRD, site.list.InfoFlora.family$ZHRE, site.list.InfoFlora.family$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.family$BSUA, site.list.InfoFlora.family$BSUB, site.list.InfoFlora.family$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.family$BSRD, site.list.InfoFlora.family$BSRE, site.list.InfoFlora.family$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.family$BEUA, site.list.InfoFlora.family$BEUB, site.list.InfoFlora.family$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.family$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.family$ZHUA, site.list.gbif.family$ZHUB, site.list.gbif.family$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.family$ZHRD, site.list.gbif.family$ZHRE, site.list.gbif.family$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.family$BSUA, site.list.gbif.family$BSUB, site.list.gbif.family$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.family$BSRD, site.list.gbif.family$BSRE, site.list.gbif.family$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.family$BEUA, site.list.gbif.family$BEUB, site.list.gbif.family$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.family$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.family.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all regions
intersection.regions.family <- list()
table.regions.family <- list()
correlation.region.family <- data.frame(family = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common family of two regions
    intersection.regions.family[[paste(i,"/",j, sep="")]] <- intersect(region.family.occ[[i]], region.family.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.regions.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.family <- rbind(correlation.region.family, correlation)
    
  } # end loop j
} # end loop i

correlation.region.family$region1 <- as.factor(correlation.region.family$region1)
matrix.region.family <- split(correlation.region.family, f = correlation.region.family$region1)
temp <- list.cbind(matrix.region.family)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.family$region1)
rownames(temp) <- regionnames
matrix.region.family <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.family) <- "numeric"

# store in list
lapi.region[["family"]] <- abs(matrix.region.family)

### plotting ----
library(ggcorrplot)
axis.order.reg <- c("BER", "BSR", "ZHR", "BEU", "BSU", "ZHU") # to consistenly order the axis

d1 <- ggcorrplot(lapi.region[["species"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")


d2 <- ggcorrplot(lapi.region[["genus"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

d3 <- ggcorrplot(lapi.region[["family"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

# arrange them into one file to export
setwd(output)
plot4 <- ggarrange(d1,d2,d3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot4, top = text_grob("B.lapidarius", face = "bold", size = 22))
ggsave("./04_Goal_3/corr_lapi_region.png", width = 10, height = 30)
setwd(input)




