################################################
# R data preparation PLANTS first looks
# by Simonetta Selva
#
# Created: October 19, 2022
# Project: Bumblebee 2022
################################################

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)


# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
plants <- read.csv2("phenology_pollen_nectar_sugar_database_copy.csv",sep = ",")
plants <- read.delim("phenology_pollen_nectar_sugar_database_copy.txt")
