################################################
# R compute Functional Diversity
# by Simonetta Selva
#
# Created: December 14, 2022
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
BB22.full <- read_csv("BB22_full.csv")















