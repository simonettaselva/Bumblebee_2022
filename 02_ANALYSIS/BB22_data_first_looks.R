################################################
# R data first looks
# by Simonetta Selva
#
# Created: Nove,ber 17, 2022
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


# load data
setwd(input)
BB22.abund <- read.csv("BB22.abund.csv") #from BB22_data_preparation

#### SHANNON AND S across locations ####
# calculate Shannon and Number of Species by individual
library(vegan)

BB22.shannon <- BB22.abund %>% 
  group_by(ID) %>%
  summarise(site = site,
            Shannon = diversity(Abundance),
            NrSpecies=n_distinct(species),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, sep="_"))) %>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F"),
         landscape = fct_relevel(landscape, "U", "R")) %>%
  distinct()

# write.csv(BB22.shannon, "BB22_s_h.csv")

#Densitiy plots of Shannon and Nr Species
density.NrSpecies <- ggdensity(BB22.shannon, x = "NrSpecies", 
                               fill = "bbspecies", palette = "jco")
density.Shannon <- ggdensity(BB22.shannon, x = "Shannon", 
                             fill = "bbspecies", palette = "jco")

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
##ggsave("Shannon_Landscapes_Locations.png", width = 8, height = 14)

plot2 <- ggarrange(s1, s2, s3, nrow = 3, labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot2, top = text_grob("Number of Species across Landscapes and Locations", 
                                       face = "bold", size = 14))
##ggsave("NrSpecies_Landscapes_Locations.png", width = 8, height = 14)
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
##ggsave("Bern_overview.png", width = 8, height = 14)
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
##ggsave("Basel_overview.png", width = 8, height = 14)
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
##ggsave("Zurich_overview.png", width = 8, height = 14)
setwd(input)


#### Compare and look at leg/body differences ####

# filter for entries only on body per bbspecies
BB22.lapi.body <- BB22.abund %>% 
  filter(str_detect(bborgan, "B") & bbspecies == "B.lapidarius")
BB22.lapi.body.OTU <- BB22.lapi.body %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.body <- BB22.abund %>% 
  filter(str_detect(bborgan, "B") & bbspecies == "B.pascuorum")
BB22.pasc.body.OTU <- BB22.pasc.body %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

# filter for entries only on legs
BB22.lapi.leg <- BB22.abund %>% 
  filter(str_detect(bborgan, "L") & bbspecies == "B.lapidarius")
BB22.lapi.leg.OTU <- BB22.lapi.leg %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.leg <- BB22.abund %>% 
  filter(str_detect(bborgan, "L") & bbspecies == "B.pascuorum")
BB22.pasc.leg.OTU <- BB22.pasc.leg %>% 
  group_by(OTU) %>%
  summarise(cum.abund = sum(Abundance))

# find common visited plant species body
shared.pasc <- length(intersect(BB22.pasc.leg.OTU$OTU, BB22.pasc.body.OTU$OTU)) #shared species B.p.
visited.pasc <- length(unique(as_factor(c(BB22.pasc.leg.OTU$OTU, BB22.pasc.body.OTU$OTU)))) #all species B.p.

shared.lapi <- length(intersect(BB22.lapi.leg.OTU$OTU, BB22.lapi.body.OTU$OTU)) #shared species B.l.
visited.lapi <- length(unique(as_factor(c(BB22.lapi.leg.OTU$OTU, BB22.lapi.body.OTU$OTU)))) #all species B.l.

df.ratio <- data.frame (percent = c(shared.pasc/visited.pasc, (visited.pasc-shared.pasc)/visited.pasc, 
                                    shared.lapi/visited.lapi, (visited.lapi-shared.lapi)/visited.lapi),
                        bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                        plant_species = c("shared", "all (body+leg)", "shared", "all (body+leg)"))

# Stacked + percent
p1 <- ggplot(df.ratio, aes(fill=plant_species, y=percent, x=bbspecies)) + 
  geom_bar(position="fill", stat="identity") + xlab("") +
  ggtitle("Body and Leg: Ratio of Shared Plant Species detected")+ theme_classic(base_size=20)

# plot sum of plant species visited
df.sum <- data.frame(sum = c(length(BB22.pasc.leg.OTU$OTU), length(BB22.pasc.body.OTU$OTU), length(BB22.lapi.leg.OTU$OTU), length(BB22.lapi.body.OTU$OTU)),
                     bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                     bborgan = c("leg", "body","leg", "body"))
p2 <- ggplot(df.sum, aes(x=bbspecies, y=sum, fill = bborgan)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Sum of Plant Species detected") + xlab("")+
  ggtitle("Body and Leg: Sum of Plant Species detected")+ theme_classic(base_size=20)

setwd(output)
plot1 <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
annotate_figure(plot1, top = text_grob("Comparison Body and Leg Pollen", 
                                       face = "bold", size = 14))
##ggsave("Comparison_Body_Leg_Pollen.png", width = 16, height = 8)
setwd(input)


#### Compare and look at leg/body differences and landscapes ####
# BB22.OTU <- BB22.abund %>% 
#   group_by(OTU, landscape, bbspecies, bborgan) %>%
#   summarise(cum.abund = sum(Abundance))


# URBAN
BB22.abund.urban <- BB22.abund %>%
  filter(landscape == "U")

# filter for entries only on body per bbspecies
BB22.lapi.body <- BB22.abund.urban %>%
  filter(str_detect(bborgan, "B") & bbspecies == "B.lapidarius")
BB22.lapi.body.OTU.urban <- BB22.lapi.body %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.body <- BB22.abund.urban %>%
  filter(str_detect(bborgan, "B") & bbspecies == "B.pascuorum")
BB22.pasc.body.OTU.urban <- BB22.pasc.body %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

# filter for entries only on legs
BB22.lapi.leg <- BB22.abund.urban %>%
  filter(str_detect(bborgan, "L") & bbspecies == "B.lapidarius")
BB22.lapi.leg.OTU.urban <- BB22.lapi.leg %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.leg <- BB22.abund.urban %>%
  filter(str_detect(bborgan, "L") & bbspecies == "B.pascuorum")
BB22.pasc.leg.OTU.urban <- BB22.pasc.leg %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

# RURAL
BB22.abund.rural <- BB22.abund %>%
  filter(landscape == "R")

# filter for entries only on body per bbspecies
BB22.lapi.body <- BB22.abund.rural %>%
  filter(str_detect(bborgan, "B") & bbspecies == "B.lapidarius")
BB22.lapi.body.OTU.rural <- BB22.lapi.body %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.body <- BB22.abund.rural %>%
  filter(str_detect(bborgan, "B") & bbspecies == "B.pascuorum")
BB22.pasc.body.OTU.rural <- BB22.pasc.body %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

# filter for entries only on legs
BB22.lapi.leg <- BB22.abund.rural %>%
  filter(str_detect(bborgan, "L") & bbspecies == "B.lapidarius")
BB22.lapi.leg.OTU.rural <- BB22.lapi.leg %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

BB22.pasc.leg <- BB22.abund.rural %>%
  filter(str_detect(bborgan, "L") & bbspecies == "B.pascuorum")
BB22.pasc.leg.OTU.rural <- BB22.pasc.leg %>%
  group_by(OTU, landscape) %>%
  summarise(cum.abund = sum(Abundance))

# find common visited plant species body per landscape
shared.pasc <- length(intersect(BB22.pasc.leg.OTU.urban$OTU, BB22.pasc.body.OTU.urban$OTU)) #shared species B.p.
visited.pasc <- length(unique(as_factor(c(BB22.pasc.leg.OTU.urban$OTU, BB22.pasc.body.OTU.urban$OTU)))) #all species B.p.

shared.lapi <- length(intersect(BB22.lapi.leg.OTU.urban$OTU, BB22.lapi.body.OTU.urban$OTU)) #shared species B.l.
visited.lapi <- length(unique(as_factor(c(BB22.lapi.leg.OTU.urban$OTU, BB22.lapi.body.OTU.urban$OTU)))) #all species B.l.

df.ratio.urban <- data.frame (percent = c(shared.pasc/visited.pasc, (visited.pasc-shared.pasc)/visited.pasc, 
                                          shared.lapi/visited.lapi, (visited.lapi-shared.lapi)/visited.lapi),
                              bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                              plant_species = c("shared", "all (body+leg)", "shared", "all (body+leg)"),
                              landscape = rep("U",4))

shared.pasc <- length(intersect(BB22.pasc.leg.OTU.rural$OTU, BB22.pasc.body.OTU.rural$OTU)) #shared species B.p.
visited.pasc <- length(unique(as_factor(c(BB22.pasc.leg.OTU.rural$OTU, BB22.pasc.body.OTU.rural$OTU)))) #all species B.p.

shared.lapi <- length(intersect(BB22.lapi.leg.OTU.rural$OTU, BB22.lapi.body.OTU.rural$OTU)) #shared species B.l.
visited.lapi <- length(unique(as_factor(c(BB22.lapi.leg.OTU.rural$OTU, BB22.lapi.body.OTU.rural$OTU)))) #all species B.l.

df.ratio.rural <- data.frame (percent = c(shared.pasc/visited.pasc, (visited.pasc-shared.pasc)/visited.pasc, 
                                          shared.lapi/visited.lapi, (visited.lapi-shared.lapi)/visited.lapi),
                              bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                              plant_species = c("shared", "all (body+leg)", "shared", "all (body+leg)"),
                              landscape = rep("R",4))

df.ratio <- rbind(df.ratio.urban, df.ratio.rural)

# Stacked + percent
p3 <- ggplot(df.ratio, aes(fill=plant_species, y=percent, x=landscape)) + 
  geom_bar(position="fill", stat="identity") + xlab("") +
  ggtitle("Body and Leg: Ratio of Shared Plant Species detected") + facet_wrap(~bbspecies)+ theme_classic(base_size=20)


# plot sum of plant species visited
df.sum.urban <- data.frame(sum = c(length(BB22.pasc.leg.OTU.urban$OTU), length(BB22.pasc.body.OTU.urban$OTU), 
                                   length(BB22.lapi.leg.OTU.urban$OTU), length(BB22.lapi.body.OTU.urban$OTU)),
                           bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                           bborgan = c("leg", "body","leg", "body"),
                           landscape = rep("U",4))
df.sum.rural <- data.frame(sum = c(length(BB22.pasc.leg.OTU$OTU), length(BB22.pasc.body.OTU$OTU), length(BB22.lapi.leg.OTU$OTU), length(BB22.lapi.body.OTU$OTU)),
                           bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                           bborgan = c("leg", "body","leg", "body"),
                           landscape = rep("R",4))

df.sum <- rbind(df.sum.urban, df.sum.rural)

p4 <- ggplot(df.sum, aes(x=landscape, y=sum, fill = bborgan)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Sum of Plant Species detected") + xlab("")+
  ggtitle("Body and Leg: Sum of Plant Species detected")+ facet_wrap(~bbspecies)+ theme_classic(base_size=20)

setwd(output)
plot1 <- ggarrange(p3, p4, ncol = 2, labels = c("A", "B"))
annotate_figure(plot1, top = text_grob("Comparison Body and Leg Pollen across Landscapes", 
                                       face = "bold", size = 14))
##ggsave("Comparison_landcape_Body_Leg_Pollen.png", width = 16, height = 8)
setwd(input)


#### Compare and look at pollen differences in landscapes #### 
BB22.abund.land <- BB22.abund %>%
  group_by(OTU, bbspecies, landscape) %>%
  summarise(cum.abund = sum(Abundance),
            landscape = landscape)%>%
  distinct()


shared.pasc.land <- length(intersect(BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.pascuorum"&BB22.abund.land$landscape=="U"], 
                                     BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.pascuorum"&BB22.abund.land$landscape=="R"])) #shared species B.p.
visited.pasc.land <- length(unique(as_factor(c(BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.pascuorum"&BB22.abund.land$landscape=="U"], 
                                               BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.pascuorum"&BB22.abund.land$landscape=="R"])))) #all species B.p.

shared.lapi.land <- length(intersect(BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.lapidarius"&BB22.abund.land$landscape=="U"], 
                                     BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.lapidarius"&BB22.abund.land$landscape=="R"])) #shared species B.p.
visited.lapi.land <- length(unique(as_factor(c(BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.lapidarius"&BB22.abund.land$landscape=="U"], 
                                               BB22.abund.land$OTU[BB22.abund.land$bbspecies=="B.lapidarius"&BB22.abund.land$landscape=="R"])))) #all species B.p.

df.ratio <- data.frame (percent = c(shared.pasc.land/visited.pasc.land, (visited.pasc.land-shared.pasc.land)/visited.pasc.land, 
                                    shared.lapi.land/visited.lapi.land, (visited.lapi.land-shared.lapi.land)/visited.lapi.land),
                        bbspecies = c("B.pascuorum", "B.pascuorum", "B.lapidarius", "B.lapidarius"),
                        plant_species = c("shared", "all (U & R)", "shared", "all (U & R)"))
# Stacked + percent
ggplot(df.ratio, aes(fill=plant_species, y=percent, x=bbspecies)) + 
  geom_bar(position="fill", stat="identity") + xlab("") +
  ggtitle("Landscapes: Ratio of Shared Plant Species detected") + theme_classic(base_size=20)

setwd(output)
##ggsave("Comparison_Landscpe_Pollen.png", width = 8, height = 8)
setwd(input)

#### Relationship Shannon/NrSpecies and phenotype traits BB #### 
BB16 <- read.csv("Bumblebee_Data_Eggenberger_et.al_JAE_2019.csv")
BB22.shannon.body <- BB22.shannon %>%
  filter(str_detect(bborgan, "B"))%>%
  mutate(ID = substring(ID,1, nchar(ID)-1))

BB22.bb.traits <- merge(BB22.shannon.body,BB16, by ="ID") %>%
  summarize(ID = ID,
            site = site,
            Shannon = Shannon,
            NrSpecies = NrSpecies,
            location = as.factor(location),
            landscape = as.factor(Habitat),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            intertegular_distance = Intertegular_distance,
            glossa = Glossa,
            prementum = Prementum,
            proboscis_length = Proboscis_length,
            proboscis_ratio = proboscis_ratio,
            fore_wing_length = fore_wing_length,
            fore_wing_ratio = fore_wing_ratio,
            corbicula_length = corbicula_length,
            corbicula_ratio = corbicula_ratio)
names(BB22.bb.traits)
levels(BB22.bb.traits$bbspecies) <- c("B.lapidarius", "B.pascuorum")

write_csv(BB22.bb.traits, "BB22.bb.traits.csv")

# plot relation Shannon/NrSpecies and phenitypic traits BB
#NrSpecies
a1 <- ggplot(BB22.bb.traits, aes(glossa, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a2 <- ggplot(BB22.bb.traits, aes(prementum, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a3 <- ggplot(BB22.bb.traits, aes(proboscis_length, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a4 <- ggplot(BB22.bb.traits, aes(proboscis_ratio, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a5 <- ggplot(BB22.bb.traits, aes(fore_wing_length, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a6 <- ggplot(BB22.bb.traits, aes(fore_wing_ratio, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a7 <- ggplot(BB22.bb.traits, aes(corbicula_length, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a8 <- ggplot(BB22.bb.traits, aes(corbicula_ratio, NrSpecies, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))

setwd(output)
plot3 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, ncol = 4, nrow=2, labels = c(LETTERS[1:8]),   common.legend = TRUE)
annotate_figure(plot3, top = text_grob("Comparison Traits and Number of Species across Landscapes", 
                                       face = "bold", size = 14))
#ggsave("Comparison_Traits_NrSpecies_Landscapes.png", width = 16, height = 8)
setwd(input)


# SHANNON
a1 <- ggplot(BB22.bb.traits, aes(glossa, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a2 <- ggplot(BB22.bb.traits, aes(prementum, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a3 <- ggplot(BB22.bb.traits, aes(proboscis_length, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a4 <- ggplot(BB22.bb.traits, aes(proboscis_ratio, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a5 <- ggplot(BB22.bb.traits, aes(fore_wing_length, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a6 <- ggplot(BB22.bb.traits, aes(fore_wing_ratio, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a7 <- ggplot(BB22.bb.traits, aes(corbicula_length, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))
a8 <- ggplot(BB22.bb.traits, aes(corbicula_ratio, Shannon, colour = landscape, shape = bbspecies, linetype = bbspecies)) + 
  geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
  scale_linetype_manual(values=c("solid", "dotted"))

setwd(output)
plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, ncol = 4, nrow=2, labels = c(LETTERS[1:8]),   common.legend = TRUE)
annotate_figure(plot4, top = text_grob("Comparison Traits and Shannon across Landscapes", 
                                       face = "bold", size = 14))
#ggsave("Comparison_Traits_Shannon_Landscapes.png", width = 16, height = 8)
setwd(input)


#### Look at metrics relationships ###
rm(list=ls())
# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/Bumblebee_2022/03_OUTPUT"

setwd(input)
BB22.metrics.body <- read_csv("BB22.metrics.csv")%>%
  filter(bborgan == "B")%>%
  mutate(ID = substring(ID,1, nchar(ID)-1))
BB22.bb.traits <- read_csv("BB22.bb.traits.csv")

BB22.metrics.traits <- merge(BB22.metrics.body, BB22.bb.traits[, -c(2:8)], by = "ID")


traits <- colnames(BB22.metrics.traits[, 16:24])
metrics <- colnames(BB22.metrics.traits[, 8:15])

for (i in metrics) {
  x <- 1
  for (j in traits) {
    assign(paste("a", x, sep=""),
           ggplot(BB22.metrics.traits, aes_string(j, i, colour = "landscape", shape = "bbspecies", linetype = "bbspecies")) + 
             geom_point(alpha = 0.4 ) + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam")+
             scale_linetype_manual(values=c("solid", "dotted")))
    x <- x+1
  }
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, ncol = 4, nrow=2, labels = c(LETTERS[1:8]),   common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("Comparison of", i, "and traits across Landscapes", sep = ""),
                                         face = "bold", size = 14))
  ggsave(paste("Correlation of", i, "and traits across Landscapes.png", sep = ""), width = 16, height = 8)
  setwd(input)
}

#with no differentiation of urban rural and species 
for (i in metrics) {
  x <- 1
  for (j in traits) {
    assign(paste("a", x, sep=""),
           ggplot(BB22.metrics.traits, aes_string(j, i)) + 
             geom_point() + theme_bw() + theme(aspect.ratio=1) + geom_smooth(method="gam"))
    x <- x+1
  }
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, ncol = 4, nrow=2, labels = c(LETTERS[1:8]),   common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("Comparison of", i, "and traits", sep = ""),
                                         face = "bold", size = 14))
  ggsave(paste("Correlation of", i, "and traits.png", sep = ""), width = 16, height = 8)
  setwd(input)
}

#Linear regression R2 and slope
for (i in metrics) {
  x <- 1
  for (j in traits) {
    f <- (paste(i,"~", j))
    fit <- lm(f, data = BB22.metrics.traits)
    assign(paste("a", x, sep=""),
           ggplot(BB22.metrics.traits, aes_string(j, i)) + 
             geom_point() + theme_bw() + theme(aspect.ratio=1) + 
             geom_smooth(method="lm", col = "red")+
             labs(subtitle = paste(
               paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "P =",signif(summary(fit)$coef[2,4], 5)),
               paste("Intercept =",signif(fit$coef[[1]],5 ),
                     "Slope =",signif(fit$coef[[2]], 5)), sep="\n"))
           )
    x <- x+1
  }
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8, ncol = 4, nrow=2, labels = c(LETTERS[1:8]),   common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("Comparison of", i, "and traits", sep = ""),
                                         face = "bold", size = 14))
  # ggsave(paste("Correlation of", i, "and traits.png", sep = ""), width = 16, height = 8)
  setwd(input)
}

#### reproduce Joans plots ####
rm(list=ls())

BB22.full <- read_csv("BB22_full.csv")
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <-paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep="")
}

#### plot plant families per site and species
setwd(output)
ggplot(BB22.full, aes(fill=family, y=Abundance, x=site)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic() + facet_wrap(~bbspecies)+ 
  ggtitle("Plant Families per Site") +
  theme(axis.text.x = element_text(angle = 90))
# ggsave(paste("PlantFamilies_per_Site.png", sep = ""), width = 16, height = 8)
setwd(input)

families.overview <- BB22.full%>%
  group_by(family)%>%
  summarise(cum.abund = sum(Abundance))
ggplot(families.overview, aes(y=reorder(family, -cum.abund), x=cum.abund)) + 
  geom_bar(stat="identity")+ theme_classic()

rare.families <- families.overview %>% top_n(nrow(families.overview)-30, -cum.abund)
BB22.full$family.agg <- BB22.full$family
for(h in rare.families$family){
  for(i in 1:nrow(BB22.full)){
    if(BB22.full$family.agg[i] == h){
      BB22.full$family.agg[i] <- "Other families"
    }
  }
}

ggplot(BB22.full, aes(fill=family.agg, y=Abundance, x=site)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic() + facet_wrap(~bbspecies)+ 
  ggtitle("Plant Families per Site") +
  theme(axis.text.x = element_text(angle = 90))
                                                                           

#### plot plant growth form per site and species
setwd(output)
ggplot(BB22.full, aes(fill=native_exotic, y=binom.abund, x=site)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic() + facet_wrap(~bbspecies)+ 
  ggtitle("Origin status") + ylab("Proportion of species in the pollen")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("OriginStatus_per_Site.png", sep = ""), width = 16, height = 8)
setwd(input)

#### growth type per site and species
setwd(output)
ggplot(BB22.full, aes(fill=growth_form_category, y=binom.abund, x=site)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic() + facet_wrap(~bbspecies)+ 
  ggtitle("Growth Form") + ylab("Proportion of species in the pollen")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("GrowthForm__per_Site.png", sep = ""), width = 16, height = 8)
setwd(input)


#### blossom class per site and species
setwd(output)
ggplot(BB22.full, aes(fill=structural_blossom_class, y=binom.abund, x=site)) + 
  geom_bar(position="fill", stat="identity")+ theme_classic() + facet_wrap(~bbspecies)+ 
  ggtitle("Blossom Class") + ylab("Proportion of species in the pollen")+
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("BlossonClass__per_Site.png", sep = ""), width = 16, height = 8)
setwd(input)


### phylogenetic tree

library(V.PhyloMaker)
species <- BB22.full$plant.species
phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

# run the function (load data since it takes a long time)
result <- phylo.maker(phylo, scenarios=c("S1","S2","S3"))
# saveRDS(example, file="phylo.RData")
# result <- readRDS("phylo.RData")

# plot the phylogenies with node ages displayed.
library(ape)
tree <- plot.phylo(result$scenario.3, cex = 0.5, main = "Phylogenetic tree of Species in Pollen") #Joan used S3
tree
# write.tree(result$scenario.3, "tree.tre")

 # bubble plot relative abundances
BB22.full.site <- BB22.full%>%
  group_by(site, plant.species, bbspecies)%>%
  summarise(Abundance = sum(Abundance))
p2 <- ggplot(BB22.full.site, aes(x = site, y =BB22.full.site$plant.species, color = site)) + 
  geom_point(aes(size = Abundance, fill = site, alpha=0.5)) + facet_wrap(~bbspecies) +theme_classic()



library(patchwork)
library(ggtree)

atree | p2

p2 %>% insert_left(tree)




p <- ggtree(result)











