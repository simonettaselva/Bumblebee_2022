setwd("~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA/sites_plant_list/05_InfoFlora_2500")
library(tidyverse)
library(dplyr)


df <- read.csv(paste("ZHUA", "_IF_2500.csv", sep = ""))

for (i in unique(BB22_full.numeric$site)) {
df <- read.csv(paste(i, "_IF_2500.csv", sep = ""))
df <- df%>% 
  filter(In.Buffer.1500 == "x")
assign(i, df)
}

write_csv(ZHUC, paste("~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA/sites_plant_list/04_InfoFlora_1500/", "ZHUB", "_IF_1500", ".csv", sep = ""))

# BERD BEUA BEUB BEUC BSRD BSRE BSRF BSUA BSUB BSUC ZHRD ZHRE ZHRF ZHUA ZHUB ZHUC