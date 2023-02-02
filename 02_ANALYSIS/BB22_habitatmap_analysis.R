################################################
# Site Analysis
# by Simonetta Selva
#
# Created: February 1st, 2023
# Project: Bumblebee 2022
################################################

rm(list=ls())
# preparation and load data ----

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
BB22.habitat <- read_csv("BB22_habitat_buffer.csv") %>%
  mutate(type = substring(as.character(TypoCH_NUM), 1, 1))
BB22.habitat$type <- as.factor(BB22.habitat$type)
levels(BB22.habitat$type) <- c("Water body", "Riparian and wetland", "Glacier, rock, rubble and scree", 
                               "Greenland", "Bush and shrubbery", "Forest", "Turf and ruderal areas", 
                               "Tree and field crops", "Building")
BB22.habitat$site <- factor(BB22.habitat$site, levels = c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF"))

head(BB22.habitat)

BB22.habitat.area <- BB22.habitat %>%
  group_by(site, type)%>%
  summarise(location = Location,
            landscape = Habitat,
            replicate = Code,
            area = sum(area)) %>%
  distinct()

# plotting
ggplot(BB22.habitat.area, aes(fill=type, y=area, x=site)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("area [m2]") + xlab("")+
  scale_fill_manual(values = c("#3470f3", "#d4f7fe", "#d8e7f1", "#f8e3cd", "#f0f04f", "#b5e674", "#c09853", "#d779de", "#a0a0a0"))+
  theme_classic(base_size = 20) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
setwd(output)
ggsave(paste("./habitat map/sites_habitat_map.png", sep = ""), width = 16, height = 10)
setwd(input)




