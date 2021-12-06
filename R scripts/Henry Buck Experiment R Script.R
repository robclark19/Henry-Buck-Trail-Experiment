# Henry Buck Experiment and Associated Data

# Libraries ####
library("tidyverse")
library("lme4")
library("car")
library("blmeco")
library("ggplot2")
library("gridExtra")
library("readxl")
library("vegan")
library("tidyverse")
library("emmeans")
library("multcomp")
library("doBy")
library("ggpubr")
library("cowplot")
library("ggrepel")
library("viridis")
library("BiodiversityR")

# To do ###
# import updated 2017 ant data
# import updated 2017 plant data
# clean plant data and make rows for myrmecochore abundance
# run non-parametric tests of treatment effects on ants
# run nmds for plant community response to experiment
# run non-parametric tests of treatment effects on myrmecochore abundance
# complete visualization for all three analyses


# 2017 Ant Data #####
ant_dat <- read.csv("./data/alchemy cookie.csv", header=TRUE)
# look for any extra data from june 24 ant baits


# Aphaenogaster colony recruitment
# simple non-parametric test
kruskal.test(A_picea ~ Treatment, data=ant_dat)

# Kruskal-Wallis chi-squared = 1.158, df = 2, p-value = 0.5605



# make a data.frame with mean and SE
ant_plot_dat <- summaryBy(A_picea ~ Treatment, data=ant_dat, FUN=c(length,mean,sd))
ant_plot_dat

# Rename column for treatment length (values) to just N
names(ant_plot_dat)[names(ant_plot_dat)=="A_picea.length"] <- "N"

#calculate standard error of the mean manually
ant_plot_dat$SE <- ant_plot_dat$A_picea.sd / sqrt(ant_plot_dat$N)
ant_plot_dat


# Fig 3 ####
# Violin plot (?) of ant abundance for Apheanogaster
ant_violin <- ggplot(ant_plot_dat, aes(x=Treatment, y=A_picea.mean)) +
  geom_point(stat="identity") +
  geom_errorbar(aes(ymin=A_picea.mean-SE, ymax=A_picea.mean+SE), position=position_dodge(0.5), width=0.2) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Treatment", y=expression('# of'~italic(Aphaenogaster)~'colonies at baits')) + 
  # scale_x_discrete(labels=c("Present", "Excluded")) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_violin(data=ant_dat, aes(x=Treatment, y=A_picea), alpha=0, adjust = 1)
ant_violin

# I'd call this figure 3 because #1 is plot design, #2 is myrmecochore % cover
# ggsave(filename = "./Figures/Fig3.svg", plot = ant_violin , device = "svg",
#       width = 4, height = 4, units = "in")


# All non-apheanogster counts ("nac")
ant_dat$nac <- ant_dat$Colony_total - ant_dat$A_picea

# no impact on the abundance of other ant species showing up at baits
kruskal.test(nac  ~ Treatment, data=ant_dat)

# Kruskal-Wallis chi-squared = 0.78385, df = 2, p-value = 0.6758










# 2017 Plant Data #######
buck_dat <- read_excel("Data/2017 henry buck trail plant community plots.xlsx") %>% 
  replace(is.na(.), 0) %>% 
  as.data.frame()



# For NMDS and specaccum make an incidence based dataset
# change to incidence based data rather than worker counts
incidence_dat <- buck_dat %>% mutate_if(is.numeric, ~1 * (. != 0))

matrix_dat <- incidence_dat
matrix_dat$Treatment_2017 <- NULL
matrix_dat$Feet <- NULL
matrix_dat$Block <- NULL

# NMDS ####

env_dat <- as.data.frame(buck_dat)
names(env_dat)[names(env_dat) == 'Treatment_2017'] <- 'sites'

# change back to matrix.dat to use just abundance data

Accum.3 <- accumcomp(matrix_dat, y=env_dat, factor='sites', 
                     method="exact", conditioned=F, plotit=F)



# possibly need for extrafont::loadfonts(device="win") to have Arial
# as alternative, use library(ggThemeAssist)
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())


#this command prepares data for ggplot
accum.long3 <- accumcomp.long(Accum.3, ci=NA, label.freq=5)

plotgg1 <- ggplot(data=accum.long3, aes(x = Sites, y = Richness, ymax =  UPR, ymin= LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(accum.long3, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "# of meters sampled", y = "Understory plant species richness", colour = "Treatment", shape = "Treatment")


# Figure S4 #####
ggsave(filename = "./Figures/FigS4.svg", plot = plotgg1 , device = "svg",
       width = 6, height = 4.5, units = "in")

# Fig S4 basic stats
str(incidence_dat)

# basic richness
# Create variable for total richness
incidence_dat_2 <- incidence_dat %>%
  rowwise() %>%
  mutate(richness = sum(across(Spring_beauty:Fuzzy_plant), na.rm = T))   %>%
  as.data.frame()

hist(incidence_dat_2$richness)

# huh, not what i expected.
plant_rich_glm <- lmer(richness ~ Treatment_2017 + (1|Feet), data=incidence_dat_2)
Anova(plant_rich_glm)
cld(emmeans(plant_rich_glm, ~ Treatment_2017))

# yeah still sig with kruskal test.
kruskal.test(richness  ~ Treatment_2017, data=incidence_dat_2)


# what about if you compare the overlap of the rarefaction curves









# Fig 2 #### 
# change to % coverage of each group among treatments
# Create variable for ant-dispersed plants
buck_dat$ant_plants <- 
  (buck_dat$Spring_beauty +
  buck_dat$Mayflower +
  buck_dat$Trillium +
  buck_dat$D_breeches +
  buck_dat$Anemone +
  buck_dat$Hepatica +
  buck_dat$Yellow_violet)

# Create variable for total coverage
buck_dat <- buck_dat %>%
rowwise() %>%
  mutate(total_plants = sum(across(Spring_beauty:Fuzzy_plant), na.rm = T))   %>%
 as.data.frame()

# Create variable for non-ant-dispersed plants
buck_dat$non_ant_plants <- buck_dat$total_plants - buck_dat$ant_plants


# Change to ant plant coverage across treatments

ant_plant_glm <- glm(ant_plants ~ Treatment_2017, data=buck_dat)
Anova(ant_plant_glm)

# preliminary analysis indicates no change to ant plant coverage
# so they recuperated over the duration of the study, 
# or seed removal didnt impact pops
plot(emmeans(ant_plant_glm,  ~ Treatment_2017, type="response"))

# plant glm
plant_glm <- glm(non_ant_plants ~ Treatment_2017, data=buck_dat)
Anova(plant_glm)

# slight increase in non-ant plants
# perhaps indicates some change in competition historically?
plot(emmeans(plant_glm,  ~ Treatment_2017, type="response"))



# percent ant plants (pap) of total sampled
# this model is great. it has the proportional data, weighted by total, etc.
# removing elaiosome bearing plants reduced their prop in the plant community
# supplementation had no real effect, control has the highest

buck_dat$pap <- buck_dat$ant_plants / buck_dat$total_plants

pap_glm <- glmer(pap ~ Treatment_2017 + (1|Feet), weights=total_plants, family=binomial, data=buck_dat)

Anova(pap_glm)
plot(emmeans(pap_glm,  ~ Treatment_2017, type="response"))

cld(emmeans(pap_glm,  ~ Treatment_2017, type="response"))

# feet = 10 meter block transect in which % cover was measured
# robust even with that random effect (edge to center)


## make a data.frame with mean and SE
# plant_plot_dat <- buck_dat %>% drop_na() 
# plant_plot_dat<- summaryBy(pap ~ Treatment_2017, data=plant_plot_dat, FUN=c(length,mean,sd))
# names(plant_plot_dat)[names(plant_plot_dat)=="pap.length"] <- "N"
# calculate standard error of the mean manually
# plant_plot_dat$SE <- plant_plot_dat$pap.sd / sqrt(plant_plot_dat$N)


# try with model estimates
plant_plot_dat <- emmeans(pap_glm,  ~ Treatment_2017) %>%
                        cld(type="response", Letters=c("abcd"))

plant_plot_dat$pap.mean <- plant_plot_dat$prob

#fix group labels so posthoc test looks nice
plant_plot_dat$.group=gsub(" ", "", plant_plot_dat$.group)

plant_jitter <- ggplot(plant_plot_dat, aes(x=Treatment_2017, y=pap.mean)) +
  geom_point(stat="identity", size=2) +
  geom_errorbar(aes(ymin=pap.mean-SE, ymax=pap.mean+SE), position=position_dodge(0.5), width=0.05) +
  theme_bw(base_size = 12) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Treatment", y="Proportion covered by ant-dispersed plants") + 
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(limits=c(0.65, 0.8)) +
  geom_text(aes(x = Treatment_2017, y = (pap.mean+SE), label = .group, hjust=-.9))
plant_jitter


# I'd call this figure 2 because you want to show the treatment evaluation first
ggsave(filename = "./Figures/Fig2.svg", plot = plant_jitter, device = "svg",
      width = 4, height = 4, units = "in")






# Richness Curves ####

# Figure S1a ####
#Comparing the diversity across all sites from 2009 to 2010 transect data

# Henry buck 2010 April 19th transect
hbt_dat <- read.csv("./data/henry buck 2010 transect.csv", header=TRUE) %>% 
  replace(is.na(.), 0)

# What is the proportional cover of ant-dispersed plants?
#

# What is the rarefaction curve at this site?
hbt_species <- hbt_dat %>% dplyr::select(spring_beauty:partridgeberry) %>% 
  mutate_if(is.numeric, ~1 * (. != 0))


hbt_curve <- specaccum(hbt_species)

FigS1a <- plot(hbt_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main="Henry Buck Trail")

#richness estimates (various approaches)
specpool(hbt_species)

poolaccum(hbt_species)

plot(poolaccum(hbt_species), display = "jack1")

estimateR(hbt_species) %>% plot()






# Figure S1b #####

# Galcoe preserve
gft_dat <- read.csv("./data/galko farm wallingford transect.csv", header=TRUE) %>% 
  replace(is.na(.), 0)

# need to drop woody plants first line-by-line
gft_dat$cherry.seedling <- NULL
gft_dat$hornbeam.seedling <- NULL
gft_dat$hornbeam.seedling <- NULL
gft_dat$hickory.seedling <- NULL
gft_dat$chestnut.oak.seedling <- NULL
gft_dat$virginia.creeper <- NULL
gft_dat$poison.ivy <- NULL
gft_dat$asiatic.bittersweet <- NULL
gft_dat$raspberry <- NULL
gft_dat$multiflora.rose.seedling <- NULL
gft_dat$burning.bush.seedling <- NULL



# What is the proportional cover of ant-dispersed plants?
#






# What is the rarefaction curve at this site?
gft_species <- gft_dat %>% dplyr::select(solomon.seal:goldenrod) %>% 
  mutate_if(is.numeric, ~1 * (. != 0))


gft_curve <- specaccum(gft_species)

FigS1b <- plot(gft_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main = "Galcoe Farm Preserve")



# Figure S1c ####
# Ragged mountain upslope
rmu_dat <- read.csv("./data/ragged mountain april 12 2010 transect.csv", header=TRUE) %>% 
  replace(is.na(.), 0)


# What is the proportional cover of ant-dispersed plants?
#


# What is the rarefaction curve at this site?
rmu_species <- rmu_dat %>% dplyr::select(trout.lily:red.trillium) %>% 
  mutate_if(is.numeric, ~1 * (. != 0))


rmu_curve <- specaccum(rmu_species)

FigS1c <- plot(rmu_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main="Ragged Mountain upslope")





# Figure S1d #####
# ragged mountain downslope
rmd_dat <- read.csv("./data/ragged mountain april 14 2010 transect.csv", header=TRUE) %>% 
  replace(is.na(.), 0)



# What is the proportional cover of ant-dispersed plants?
#


# What is the rarefaction curve at this site?
rmd_species <- rmd_dat %>% dplyr::select(trout.lily:canada.mayflower) %>% 
  mutate_if(is.numeric, ~1 * (. != 0))


rmd_curve <- specaccum(rmd_species)

FigS1d <-  plot(rmd_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main="Ragged Mountain downslope")





# Combine Fig S1 ####
# FigS1abcd <- ggarrange(FigS1a, FigS1b, FigS1c, FigS1d, labels = c("A", "B", "C", "D"), nrow = 2)


ggsave(filename= "./Figures/FigS1a.jpeg", plot = FigS1a , device = "jpeg", 
       width=3, height=3, units = "in")

