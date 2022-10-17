################################################
# R data preparation
# by Simonetta Selva
#
# Created: October 17, 2022
# Project: Bumblebee 2022
################################################

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22 <- read.csv("BB22_data_tres_0.01.csv")

#create binom. abundance variable
BB22 <- BB22%>% 
  mutate(binom.abund = if_else(Abundance == 0, 0, 1))

#remove entries with abundance = 0
BB22.abund <- BB22%>% 
  filter(binom.abund == 1)

