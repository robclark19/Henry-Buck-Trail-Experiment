###############################################################################

# Henry Buck Experiment and Associated Data

###############################################################################

library("tidyverse")
library("lme4")
library("car")
library("blmeco")
library("ggplot2")
library("gridExtra")
library("piecewiseSEM")
library("readxl")

# To do ###
# import updated 2017 ant data
# import updated 2017 plant data
# clean plant data and make rows for myrmecochore abundance
# run non-parametric tests of treatment effects on ants
# run nmds for plant community response to experiment
# run non-parametric tests of treatment effects on myrmecochore abundance
# complete visualization for all three analyses


# 2017 Ant Data #####
ant.dat <- read.csv("./data/alchemy cookie.csv", header=TRUE)
# add data from june 24 ant baits

str(ant.data)
# check format of ant.data

# Aphaenogaster counts
ap.model <- glm(A_picea ~ Treatment, data=ant.data)
summary(ap.model)
Anova(ap.model)

kruskal.test(A_picea ~ Treatment, data=ant.data)

# honestly why not use the simplest test?
kruskal.test(sum.ap ~ Treatment, data=ant.data)

# All colony counts
all.model <- glm(Colony_total ~ Treatment, data=ant.data)
Anova(all.model)

# carpenter ants
cp.model <- glm(C_pennsylvanicus ~ Treatment, data=ant.data)
Anova(cp.model)


# 2017 Plant Data #######
buck_plants_dat <- read_excel("Data/2017 henry buck trail plant community plots.xlsx") %>% replace(is.na(.), 0)



