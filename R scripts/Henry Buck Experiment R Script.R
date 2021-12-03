###############################################################################

# Henry Buck Experiment and Associated Data

###############################################################################

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


# Ant figure ####
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

# Figure 2 ####
# Impact of treatment on host plant diversity #######
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
  labs(x = "# of Meter Transects", y = "Plant Species Richness", colour = "Treatment", shape = "Treatment")

plotgg1


# Figure 3 ### 
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



# Figure S1 ####
#Comparing the diversity across all sites including non-henry-buck sites

# buck_dat will have to be merged with other sites and have the same plant names
# used across all, including morphospecies names
## or each can be analyzed separately as a unique own matrix


