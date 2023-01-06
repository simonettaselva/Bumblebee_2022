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

#only B.pascuroum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with funtional diversity of plants per site
BB22.fun.site <- read_csv("./FD/fd_B.pascuorum_site.csv")%>% 
  rename_with(.cols = 1, ~"site")

# compute means of bumblebee traits to compare with functional metrics of plant traits on site level
BB22.bb.traits.site <- BB22.bb.traits.sp %>%
  group_by(site) %>%
  summarise(site = site,
            location = location,
            landscape = landscape,
            intertegular_distance = mean(intertegular_distance),
            glossa = mean(glossa),
            prementum = mean(prementum),
            proboscis_length = mean(proboscis_length),
            proboscis_ratio = mean(proboscis_ratio),
            fore_wing_length = mean(fore_wing_length),
            fore_wing_ratio = mean(fore_wing_ratio),
            corbicula_length = mean(corbicula_length),
            corbicula_ratio = mean(corbicula_ratio)) %>%
  distinct()

# add site coordinates to the trait data frame (in LV95)
BB22.sites <- merge(BB22.bb.traits.site, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 
BB22.sites <- merge(BB22.sites, BB22.fun.site, by  = "site", all.x=TRUE)

# look at data
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.sites[, -c(1:14)])

# Boxplots for all the variables we want to look at
resp <- c("sp_richn","fdis", "fric", "fdiv", "feve", "fspe")
library(psych)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
for (j in resp) {
  gg.data <- data.frame(landscape=BB22.sites$landscape,value=BB22.sites[[j]])
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    ylab(j) + xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none")
  plot_list[[j]] <- p
  describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                   plot_list[[3]],plot_list[[4]],
                   plot_list[[5]],plot_list[[6]], 
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D", "E", "F"))
annotate_figure(plot, top = text_grob("B.pascuorum: species richness and funtional diversity across landscapes", 
                                       face = "bold", size = 22))
ggsave("./functional diversity/FD_B.pascuorum.png", width = 4, height = 24)
setwd(input)


#### only B.lapidarius ####
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.lapidarius",]


# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with funtional diversity of plants per site
BB22.fun.site <- read_csv("./FD/fd_B.lapidarius_site_.csv")%>% 
  rename_with(.cols = 1, ~"site")

# compute means of bumblebee traits to compare with functional metrics of plant traits on site level
BB22.bb.traits.site <- BB22.bb.traits.sp %>%
  group_by(site) %>%
  summarise(site = site,
            location = location,
            landscape = landscape,
            intertegular_distance = mean(intertegular_distance),
            glossa = mean(glossa),
            prementum = mean(prementum),
            proboscis_length = mean(proboscis_length),
            proboscis_ratio = mean(proboscis_ratio),
            fore_wing_length = mean(fore_wing_length),
            fore_wing_ratio = mean(fore_wing_ratio),
            corbicula_length = mean(corbicula_length),
            corbicula_ratio = mean(corbicula_ratio)) %>%
  distinct()

# add site coordinates to the trait data frame (in LV95)
BB22.sites <- merge(BB22.bb.traits.site, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 
BB22.sites <- merge(BB22.sites, BB22.fun.site, by  = "site", all.x=TRUE)

# look at data
# FDis, FRic, FDiv, FEve, FSpe
library(Hmisc)
hist.data.frame(BB22.sites[, -c(1:14)])

# Boxplots for all the variables we want to look at
resp <- c("sp_richn","fdis", "fric", "fdiv", "feve", "fspe")
library(psych)
plot_list  <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  gg.data <- data.frame(landscape=BB22.sites$landscape,value=BB22.sites[[j]])
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    ylab(j) + xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) + guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none")
  plot_list[[j]] <- p
  describeBy(gg.data$value, gg.data$landscape)
}

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  plot_list[[5]],plot_list[[6]], 
                  ncol = 1, nrow = 6,
                  labels = c("A", "B", "C", "D", "E", "F"))
annotate_figure(plot, top = text_grob("B.lapidarius: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
ggsave("./functional diversity/FD_B.lapidarius.png", width = 4, height = 24)
setwd(input)



# comparison between species



