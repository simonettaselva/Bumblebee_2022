################################################
# R compute Functional Diversity Bumblebees
# by Simonetta Selva
#
# Created: February 4th, 2022
# Project: Bumblebee 2022
################################################
rm(list=ls())

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
BB22.bb.traits <- read_csv("BB22_traits.csv")

# data preparation -----
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
  mutate(site = paste(location, landscape, replicate, sep = ""),
         region = paste(location, landscape, sep = "")) %>%
  dplyr::select(-NrSpecies, -Shannon)
BB22.bb.traits$site <- as.factor(BB22.bb.traits$site)
BB22.bb.traits$ID <- as.factor(BB22.bb.traits$ID)

# B.pascuorum ----
# select only B.pascuorum
BB22.bb.traits.species <- BB22.bb.traits%>%
  filter(bbspecies == "B.pascuorum")

# prepare functional matrix
traits <- scale(BB22.bb.traits.species[, c(11,13,15)],center=TRUE,scale=TRUE)
rownames(traits)  <- BB22.bb.traits.species$ID

# prepare ecological matrix
library(reshape2)
sp.pa <- dcast(BB22.bb.traits.species, site ~ ID, length)[,-1]
rownames(sp.pa)  <- levels(BB22.bb.traits.species$site)
sp.pa[is.na(sp.pa)] <- 0
sp.pa <- as.matrix(sp.pa)

# compute FDs
library(FD)
fd.list <- FD::dbFD(x = traits , a = sp.pa) # not weighted
fd.bb <- data_frame(site = names(fd.list$nbsp),
                    nbsp.bb = fd.list$nbsp,
                    FRic.bb = fd.list$FRic,
                    FEve.bb = fd.list$FEve,
                    FDiv.bb = fd.list$FDiv) %>% 
  mutate(landscape = substring(site, 3,3))

# import FD of plants
fd.plants.pasc <-  read_csv(paste("./FD/FD_package_B.pascuorum_site.csv", sep = "")) %>% 
  rename(site = ...1,
         nbsp = nbsp.w,
         FRic = FRic.w,
         FEve = FEve.w,
         FDiv = FDiv.w) 

fd.all <- merge(fd.bb, fd.plants, by  = "site", all.x=TRUE)

# relationship bumblebee FD and plant FD
fds <- c("FRic","FEve","FDiv")
plot_list <- list()
library(nlme)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (i in fds) {
  i.bb <- paste(i, ".bb", sep = "")
  f <- formula(paste(i,"~", i.bb))
  fit <- lme(f, random=~1|landscape, data = fd.all, na.action=na.omit)
  p <- ggplot(fd.all, aes_string(i.bb, i, colour = "landscape")) + 
           geom_point() + 
           theme_classic(base_size = 20) + 
           theme(aspect.ratio=1) + 
           geom_smooth(method="lm", se = FALSE) +
           scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
           stat_cor(aes(color = landscape), size = 5)
  if (i == "FRic") {
    p <- p + ylim(0, 82) + ylab("funtional richness plants")+
      xlim(0, 80) + xlab("funtional richness bumblebees")
  } else if (i == "FDiv") {
    p <- p + ylim(0.6, 1) + ylab("funtional divergence plants")+
      xlim(0.6, 0.9) + xlab("funtional divergence bumblebees")
  } else {
    p <- p + ylim(0.2, 0.75) + ylab("funtional evenness plants")+
      xlim(0.6, 0.9) + xlab("funtional evenness bumblebees")
  }
  plot_list[[i]] <- p
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],
                  ncol = 1, nrow = 3,
                  labels = c("A", "B", "C", "D"),
                  common.legend = T)
annotate_figure(plot, top = text_grob("B.pascuorum: bumblebee FDs vs. plant FDs", 
                                      face = "bold", size = 22))
ggsave("./FD bumblebees/FD_B.pascuorum_bb_plants.png", width = 6, height = 24)
setwd(input)

# B.lapidarius ----
# select only B.lapidarius
BB22.bb.traits.species <- BB22.bb.traits%>%
  filter(bbspecies == "B.lapidarius")

# prepare functional matrix
traits <- scale(BB22.bb.traits.species[, c(11,13,15)],center=TRUE,scale=TRUE)
rownames(traits)  <- BB22.bb.traits.species$ID
colnames(BB22.bb.traits.species)

# prepare ecological matrix
library(reshape2)
sp.pa <- dcast(BB22.bb.traits.species, site ~ ID, length)[,-1]
rownames(sp.pa)  <- levels(BB22.bb.traits.species$site)
sp.pa[is.na(sp.pa)] <- 0
sp.pa <- as.matrix(sp.pa)

# compute FDs
library(FD)
fd.list <- FD::dbFD(x = traits , a = sp.pa) # not weighted
fd.bb <- data_frame(site = names(fd.list$nbsp),
                    nbsp.bb = fd.list$nbsp,
                    FRic.bb = fd.list$FRic,
                    FEve.bb = fd.list$FEve,
                    FDiv.bb = fd.list$FDiv) %>% 
  mutate(landscape = substring(site, 3,3))

# import FD of plants
fd.plants <-  read_csv(paste("./FD/FD_package_B.lapidarius_site.csv", sep = "")) %>% 
  rename(site = ...1,
         nbsp = nbsp.w,
         FRic = FRic.w,
         FEve = FEve.w,
         FDiv = FDiv.w) 
fd.all <- merge(fd.bb, fd.plants, by  = "site", all.x=TRUE)

# relationship bumblebee FD and plant FD
fds <- c("FRic","FEve","FDiv")
plot_list <- list()
library(nlme)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (i in fds) {
  i.bb <- paste(i, ".bb", sep = "")
  f <- formula(paste(i,"~", i.bb))
  fit <- lme(f, random=~1|landscape, data = fd.all, na.action=na.omit)
  p <- ggplot(fd.all, aes_string(i.bb, i, colour = "landscape")) + 
    geom_point() + 
    theme_classic(base_size = 20) + 
    theme(aspect.ratio=1) + 
    geom_smooth(method="lm", se = FALSE) +
    scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
    stat_cor(aes(color = landscape), size = 5)
  if (i == "FRic") {
    p <- p + ylim(0, 82) + ylab("funtional richness plants")+
      xlim(0, 80) + xlab("funtional richness bumblebees")
  } else if (i == "FDiv") {
    p <- p + ylim(0.6, 1) + ylab("funtional divergence plants")+
      xlim(0.6, 0.9) + xlab("funtional divergence bumblebees")
  } else {
    p <- p + ylim(0.2, 0.75) + ylab("funtional evenness plants")+
      xlim(0.6, 0.9) + xlab("funtional evenness bumblebees")
  }
  plot_list[[i]] <- p
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],
                  ncol = 1, nrow = 3,
                  labels = c("A", "B", "C", "D"),
                  common.legend = T)
annotate_figure(plot, top = text_grob("B.lapidarius: bumblebee FDs vs. plant FDs", 
                                      face = "bold", size = 22))
ggsave("./FD bumblebees/FD_B.lapidarius_bb_plants.png", width = 6, height = 24)
setwd(input)

# species combined no landscape ----
