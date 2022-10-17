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
library(ggpubr)


# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

# load data
setwd(input)
BB22 <- read.csv("BB22_data_tres_0.01.csv")

#create binom. abundance variable
BB22 <- BB22%>% 
  mutate(binom.abund = if_else(Abundance == 0, 0, 1),
         ID = substring(Sample, 2))

#remove entries with abundance = 0
BB22.abund <- BB22%>% 
  filter(binom.abund == 1)

# calculate Shannon and Number of Species by individual
library(vegan)

BB22.shannon <- BB22.abund %>% 
  group_by(ID) %>%
  summarise(Shannon = diversity(Abundance),
            NrSpecies=n_distinct(species),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, sep="_")))%>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F"),
         landscape = fct_relevel(landscape, "U", "R")) %>%
  distinct()

#### SHANNON AND S across locations ####
# plot mean shannon and S vs. landscape and location
sh1 <- ggplot(BB22.shannon, aes(Shannon, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("BODY and CUTICULA")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

s1<- ggplot(BB22.shannon, aes(NrSpecies, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("BODY and CUTICULA")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

# only body pollen
BB22.shannon.body <- BB22.shannon%>%
  filter(bborgan == "B")

sh2 <- ggplot(BB22.shannon.body, aes(Shannon, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("BODY")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

s2 <- ggplot(BB22.shannon.body, aes(NrSpecies, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("BODY")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

# only leg pollen
BB22.shannon.leg <- BB22.shannon%>%
  filter(bborgan == "L")

sh3 <- ggplot(BB22.shannon.leg, aes(Shannon, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("CUTICULA")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

s3 <- ggplot(BB22.shannon.leg, aes(NrSpecies, landscape, fill=location)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("CUTICULA")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies)

setwd(output)
plot1 <- ggarrange(sh1, sh2, sh3, nrow = 3, labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("Alpha Diversity across Landscapes and Locations", 
                                      face = "bold", size = 14))
ggsave("Shannon_Landscapes_Locations.png", width = 8, height = 14)

plot2 <- ggarrange(s1, s2, s3, nrow = 3, labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot2, top = text_grob("Number of Species across Landscapes and Locations", 
                                       face = "bold", size = 14))
ggsave("NrSpecies_Landscapes_Locations.png", width = 8, height = 14)
setwd(input)

#### SHANNON AND S in locations ####
# plot mean shannon and S vs. location
colors_fill <- scale_fill_manual(values = c("A" = "darkslategray3",
                             "B" = "darkslategray4",
                             "C" = "darkslategray",
                             "D" = "cornsilk",
                             "E" = "cornsilk2",
                             "F" = "cornsilk4"))
# BERN
BB22.shannon.bern <- BB22.shannon%>%
  filter(location == "BE")

bern1 <- ggplot(BB22.shannon.bern, aes(Shannon, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("Alpha Diversity")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

bern2 <- ggplot(BB22.shannon.bern, aes(NrSpecies, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("Number of Species")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

setwd(output)
plot1 <- ggarrange(bern1, bern2, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("Bern", 
                                       face = "bold", size = 14))
ggsave("Bern_overview.png", width = 8, height = 14)
setwd(input)


# BASEL
BB22.shannon.basel <- BB22.shannon%>%
  filter(location == "BS")

basel1 <- ggplot(BB22.shannon.basel, aes(Shannon, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("Alpha Diversity")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

basel2 <- ggplot(BB22.shannon.basel, aes(NrSpecies, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("Number of Species")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

setwd(output)
plot1 <- ggarrange(basel1, basel2, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("Basel", 
                                       face = "bold", size = 14))
ggsave("Basel_overview.png", width = 8, height = 14)
setwd(input)



# ZUERICH
BB22.shannon.zurich <- BB22.shannon%>%
  filter(location == "ZH")

zurich1 <- ggplot(BB22.shannon.zurich, aes(Shannon, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Shannon Index") +
  ylab("Landscape") + ggtitle("Alpha Diversity")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

zurich2 <- ggplot(BB22.shannon.zurich, aes(NrSpecies, location, fill=replicate)) +
  geom_boxplot(notch = T) + xlab("Number of Species") +
  ylab("Landscape") + ggtitle("Number of Species")+
  coord_flip() + theme_bw() + facet_wrap(~bbspecies) + colors_fill

setwd(output)
plot1 <- ggarrange(zurich1, zurich2, nrow = 2, labels = c("A", "B"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("Zurich", 
                                       face = "bold", size = 14))
ggsave("Zurich_overview.png", width = 8, height = 14)
setwd(input)


#### Compare and look at leg/body differences ####

# filter for entries only on body per bbspecies
BB22.lapi.body <- BB22.abund %>% 
  filter(str_detect(bborgan, "B") & bbspecies == "l")
BB22.lapi.body.OTU <- BB22.lapi.body %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.body <- BB22.abund %>% 
  filter(str_detect(bborgan, "B") & bbspecies == "p")
BB22.pasc.body.OTU <- BB22.pasc.body %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

# filter for entries only on legs
BB22.lapi.leg <- BB22.abund %>% 
  filter(str_detect(bborgan, "L") & bbspecies == "l")
BB22.lapi.leg.OTU <- BB22.lapi.leg %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.leg <- BB22.abund %>% 
  filter(str_detect(bborgan, "L") & bbspecies == "p")
BB22.pasc.leg.OTU <- BB22.pasc.leg %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

# find common visited plant species body
shared.pasc <- length(intersect(BB22.pasc.leg.OTU$OTU, BB22.pasc.body.OTU$OTU)) #shared species B.p.
visited.pasc <- length(unique(as_factor(c(BB22.pasc.leg.OTU$OTU, BB22.pasc.body.OTU$OTU)))) #all species B.p.

shared.lapi <- length(intersect(BB22.lapi.leg.OTU$OTU, BB22.lapi.body.OTU$OTU)) #shared species B.l.
visited.lapi <- length(unique(as_factor(c(BB22.lapi.leg.OTU$OTU, BB22.lapi.body.OTU$OTU)))) #all species B.l.

df.ratio <- data.frame (percent = c(shared.pasc, visited.pasc, shared.lapi, visited.lapi),
                  bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                  plant_species = c("shared", "all (body+leg)", "shared", "all (body+leg)"))
# Stacked + percent
p1 <- ggplot(df.ratio, aes(fill=plant_species, y=percent, x=bbspecies)) + 
  geom_bar(position="fill", stat="identity") + xlab("") +
  ggtitle("Body and Leg: Ratio of Shared Plant Species detected")

# plot sum of plant species visited
df.sum <- data.frame(sum = c(length(BB22.pasc.leg.OTU$OTU), length(BB22.pasc.body.OTU$OTU), length(BB22.lapi.leg.OTU$OTU), length(BB22.lapi.body.OTU$OTU)),
                     bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                     bborgan = c("leg", "body","leg", "body"))
p2 <- ggplot(df.sum, aes(x=bbspecies, y=sum, fill = bborgan)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Sum of Plant Species detected") + xlab("")+
  ggtitle("Body and Leg: Sum of Plant Species detected")

setwd(output)
plot1 <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
annotate_figure(plot1, top = text_grob("Comparison Body and Leg Pollen", 
                                       face = "bold", size = 14))
ggsave("Comparison_Body_Leg_Pollen.png", width = 16, height = 8)
setwd(input)


#### Compare and look at leg/body differences and landscapes ####





