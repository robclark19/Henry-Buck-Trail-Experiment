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
ant_plot_dat <- 
  summaryBy(A_picea ~ Treatment, data=ant_dat, FUN=c(length,mean,sd))
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

# get summary statistics for nac (non apheaenogaster colonies)

# make a data.frame with mean and SE
nac_dat <- 
  summaryBy(nac ~ Treatment, data=ant_dat, FUN=c(length,mean,sd))
nac_dat

# Rename column for treatment length (values) to just N
names(nac_dat)[names(nac_dat)=="nac.length"] <- "N"

#calculate standard error of the mean manually
nac_dat$SE <- nac_dat$nac.sd / sqrt(nac_dat$N)
nac_dat








# 2017 Plant Data #######
buck_dat <- read_excel("Data/2017 henry buck trail plant community plots.xlsx") %>% 
  replace(is.na(.), 0) %>% 
  as.data.frame()



# For NMDS and specaccum make an incidence based dataset
# change to incidence based data rather than plant cm coverage
incidence_dat <- buck_dat %>% mutate_if(is.numeric, ~1 * (. != 0))

# make a drop column
incidence_dat <- incidence_dat %>%
  rowwise() %>%
  mutate(total_plants = sum(across(Spring_beauty:Fuzzy_plant), na.rm = T))   %>%
  as.data.frame()

incidence_dat <- subset(incidence_dat, total_plants != 0)
incidence_dat$total_plants <- NULL

# Make species matrix
matrix_dat <- incidence_dat
matrix_dat$Treatment_2017 <- NULL
matrix_dat$Feet <- NULL
matrix_dat$Block <- NULL

# Must remove some transect rows that are all zeroes
matrix_dat <- matrix_dat %>% 
  filter_all(any_vars(. != 0))

# Might need to do some pooling or removals

# remove species columns with less than 5
#matrix_dat$Bloodroot <- NULL
#matrix_dat$Fern_plant <- NULL
#matrix_dat$Cinquefoil <- NULL
#matrix_dat$Hepatica <- NULL
#matrix_dat$Fern_plant <- NULL
#matrix_dat$Anemone <- NULL
#matrix_dat$Wood_aster <- NULL
#matrix_dat$Starflower <- NULL
#matrix_dat$Hepatica <- NULL
#matrix_dat$Jumpseed <- NULL
#matrix_dat$Cleavers <- NULL
#matrix_dat$Carrot_sp <- NULL
#matrix_dat$Yellow_violet <- NULL
#matrix_dat$house_plant_flower <- NULL
#matrix_dat$Allium<- NULL
#matrix_dat$Stinging_nettle <- NULL
#matrix_dat$Wide_carex <- NULL
#matrix_dat$Fuzzy_plant <- NULL
#matrix_dat$Indian_cucumber <- NULL

# NMDS
hb_nmds <- metaMDS(matrix_dat, k=3, plot=TRUE)
hb_nmds


# Make environmental matrix
# env_dat will be the matrix with env variables added back in to match length
attach(incidence_dat)

# or figure out what ten rows are dropped


#just plot points (sites and species)
plot(hb_nmds)

#ordination plot (starts with blank ordination that layers can be added on)
ordiplot(hb_nmds, type="none")
#species plot (add insect species first)
orditorp(hb_nmds,display="species",col="black",air=0.2,cex=1)
# draws a shape around it based on the environmental variable of interest
ordihull(hb_nmds, groups=Treatment_2017, draw="polygon",col="grey90",label=T)


# Hypothesis test
# actual statistical test to for urban vs. rural has different community structure
hb_fit <- envfit(hb_nmds ~ Treatment_2017, perm=999)
hb_fit

# p = 0.312 r2 = 0.0218 ok cool.








# Fig 1 HB Pie Charts ####
str(buck_dat)

# Create new dataframe for calculating total coverage of ant plants
pie_dat <- buck_dat

pie_dat$ant_plants <- 
  (pie_dat$Spring_beauty +
     pie_dat$Trout_lily +
     pie_dat$Mayflower +
     pie_dat$Trillium +
     pie_dat$D_breeches +
     pie_dat$Anemone +
     pie_dat$Hepatica +
     pie_dat$Yellow_violet)

# most common myrmecochores are:
# Spring beauty, trout lily, dutchman's, trillium



# Non-myrmecochores
# Create variable for total coverage
pie_dat <- pie_dat %>%
  rowwise() %>%
  mutate(total_plants = sum(across(Spring_beauty:Fuzzy_plant), na.rm = T))   %>%
  as.data.frame()

# Create variable for non-ant-dispersed plants
pie_dat$non_ant_plants <- pie_dat$total_plants - pie_dat$ant_plants

# Other non-common myrmecochores
# create variable for common myrmecochores
pie_dat <- pie_dat %>%
  rowwise() %>%
  mutate(total_commons =sum(across(c("Trillium", "Trout_lily","Spring_beauty","D_breeches")), na.rm = T))

# total non common ant plants
pie_dat$non_common_ant_plants <- pie_dat$ant_plants - pie_dat$total_commons





# pie chart should have 6 variables:
# Spring beauty, trout lily, dutchman's, trillium, non-myrmecochores, rare myrmecochores

pie2_dat <- subset(pie_dat, select=c("Block","Trillium", "Trout_lily","Spring_beauty","D_breeches", "non_common_ant_plants", "non_ant_plants"))


# then group into 3 to 9 circles based on these counts sorted by treatment

# summarize by each treatment combo
pie2_dat$Block <- as.factor(pie2_dat$Block)


pie3_dat <- pie2_dat %>%
  group_by(Block) %>%
  summarize(Trillium = sum(Trillium), 
            Trout_lily = sum(Trout_lily), 
            Dutchmans_breeches = sum(D_breeches), 
            Spring_beauty = sum(Spring_beauty),
            Other_myrmecohores = sum(non_common_ant_plants),
            Non_myrmecochores = sum(non_ant_plants))

# pie 1
# jesus christ i would have been done with this in an hour if i used excel

p1 <- colnames(pie3_dat) %>% as.data.frame()
p1$cover <- t(pie3_dat[1,])
p1 <- p1[-c(1), ]
names(p1)[1] <- "plant"
p1$plant <- as.factor(p1$plant)
p1$cover <- as.numeric(p1$cover)



# Make a list of the names to use in the figure:
pie_list <- c("Dutchman's breeches", "Non-myrmecochores", "Other myrmecochores",
              "Spring beauty", "Red trillium", "Trout Lily")

pie_order <- c("Dutchmans_breeches", "Trillium", "Trout_lily", "Spring_beauty", "Other_myrmecohores", "Non-myrmecochores")

pie_list_2 <- c("Dutchman's breeches", "Red trillium", "Trout lily", "Spring beauty", 
                "Other myrmecochores", "Non-myrmecochores")

#ggplot pie chart
# this treatment should be: control
p1_plot <- ggplot(p1, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
           position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name = "Plant Category", labels=pie_list_2, palette = "Accent")
p1_plot


# ok, now do the same thing but for blocks 2-9
p2 <- colnames(pie3_dat) %>% as.data.frame()
p2$cover <- t(pie3_dat[2,])
p2 <- p2[-c(1), ]
names(p2)[1] <- "plant"
p2$plant <- as.factor(p2$plant)
p2$cover <- as.numeric(p2$cover)


# ggplot pie chart
# this treatment should be: remove
p2_plot <- ggplot(p2, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p2_plot

#p3
p3 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p3$cover <- t(pie3_dat[3,])
# drop block var
p3 <- p3[-c(1), ]
# rename column header
names(p3)[1] <- "plant"
p3$plant <- as.factor(p3$plant)
p3$cover <- as.numeric(p3$cover)


# ggplot pie chart
# this treatment should be: remove
p3_plot <- ggplot(p3, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p3_plot

# p4
p4 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p4$cover <- t(pie3_dat[4,])
# drop block var
p4 <- p4[-c(1), ]
# rename column header
names(p4)[1] <- "plant"
p4$plant <- as.factor(p4$plant)
p4$cover <- as.numeric(p4$cover)


# ggplot pie chart
# this treatment should be: remove
p4_plot <- ggplot(p4, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1, check.overlap = TRUE) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p4_plot


#p5
p5 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p5$cover <- t(pie3_dat[5,])
# drop block var
p5 <- p5[-c(1), ]
# rename column header
names(p5)[1] <- "plant"
p5$plant <- as.factor(p5$plant)
p5$cover <- as.numeric(p5$cover)

# this treatment should be: control
p5_plot <- ggplot(p5, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p5_plot



#p6
p6 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p6$cover <- t(pie3_dat[6,])
# drop block var
p6 <- p6[-c(1), ]
# rename column header
names(p6)[1] <- "plant"
p6$plant <- as.factor(p6$plant)
p6$cover <- as.numeric(p6$cover)

# this treatment should be: remove
p6_plot <- ggplot(p6, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p6_plot



#p7
p7 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p7$cover <- t(pie3_dat[7,])
# drop block var
p7 <- p7[-c(1), ]
# rename column header
names(p7)[1] <- "plant"
p7$plant <- as.factor(p7$plant)
p7$cover <- as.numeric(p7$cover)

# this treatment should be: remove
p7_plot <- ggplot(p7, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p7_plot



#p8
p8 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p8$cover <- t(pie3_dat[8,])
# drop block var
p8 <- p8[-c(1), ]
# rename column header
names(p8)[1] <- "plant"
p8$plant <- as.factor(p8$plant)
p8$cover <- as.numeric(p8$cover)

# this treatment should be: add
p8_plot <- ggplot(p8, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p8_plot


#p9
p9 <- colnames(pie3_dat) %>% as.data.frame()
# change [x] to the block number
p9$cover <- t(pie3_dat[9,])
# drop block var
p9 <- p9[-c(1), ]
# rename column header
names(p9)[1] <- "plant"
p9$plant <- as.factor(p9$plant)
p9$cover <- as.numeric(p9$cover)

# this treatment should be: control
p9_plot <- ggplot(p9, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=3.1) +
  coord_polar(theta="y") +
  theme_void(base_size = 18) +
  scale_fill_brewer(name="Plant Category",labels=pie_list_2, palette="Accent")
p9_plot






# Fig 1 Export ####
p_all_labels <- c("1 | Control", 
                  "2 | Remove", 
                  "3 | Add",
                  "4 | Add",
                  "5 | Control",
                  "6 | Remove",
                  "7 | Remove",
                  "8 | Add",
                  "9 | Control")

#
pie_fig_all <- ggarrange(p1_plot, p2_plot, p3_plot, p4_plot, p5_plot, p6_plot, p7_plot, p8_plot, p9_plot,
                         labels = p_all_labels, 
                         nrow = 3, ncol = 3,
                          common.legend = TRUE, 
                         legend = "right",
                         vjust = 1.5,
                         label.x = c(0, 0.05, 0.2, 
                                     0.2, 0, 0.05,
                                     0.05, 0.2, 0))
pie_fig_all

ggsave(filename = "./Figures/Fig1.svg", plot = pie_fig_all, device = "svg",
      width = 9, height = 8, units = "in")





















# Rarefaction ####

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
# plant_plot_dat<- 
maryBy(pap ~ Treatment_2017, data=plant_plot_dat, FUN=c(length,mean,sd))
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
# ggsave(filename = "./Figures/Fig2.svg", plot = plant_jitter, device = "svg",width = 4,height = 4,units = "in")


# Fig 2 stacked #######
# ok we need to add the same groupings as Fig. 1 to this figure
# drop geom point, then draw a new geom_bar pulling from a different pooled data which contains the %

# find total coverage
stack_dat <- buck_dat %>% group_by(Treatment_2017) %>%
  summarise(ant_plants = sum(ant_plants), total_plants=sum(total_plants), non_ant_plants=sum(non_ant_plants))

# find proportional coverage
stack_dat$prop_ant <- stack_dat$ant_plants / stack_dat$total_plants
stack_dat$prop_non_ant <- stack_dat$non_ant_plants / stack_dat$total_plants

# drop values used for calcs of proportions
stack_dat = subset(stack_dat, select = c("prop_ant","prop_non_ant","Treatment_2017"))

# pivot so there is a column containing the values among plant groups and treatments
# https://tidyr.tidyverse.org/reference/pivot_longer.html
stack_dat <- stack_dat %>% pivot_longer(!Treatment_2017, names_to = "plant_group")

plant_jitter_2 <- ggplot() +
  geom_bar(data=stack_dat, stat="identity", position = position_stack(reverse = TRUE), aes(x=Treatment_2017, y=value, fill=plant_group)) +
  geom_errorbar(data=plant_plot_dat, stat="identity", aes(x=Treatment_2017, ymin=pap.mean-SE, ymax=pap.mean+SE), position=position_dodge(0.5), width=0.05) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme_bw(base_size = 12) +
  geom_text(data=plant_plot_dat, aes(x = Treatment_2017, y = (pap.mean+SE), label = .group, hjust=-.9)) +
  labs(x="Treatment", y="Proportion of plant coverage out of the entire community", fill="Plant Group") +
  scale_fill_manual(values=c("grey45", "grey"), labels=c("Ant plants","Non-ant-plants"))

plant_jitter_2


# Trillium GLMM ####
# tadd = trillium add
buck_dat$Block <- as.factor(buck_dat$Block)

trillium_count <- buck_dat


trillium_count <- trillium_count %>%
  group_by(Block, Treatment_2017) %>%
  summarize(Trillium = sum(Trillium))


tadd_glm <- glmer.nb(Trillium ~ Treatment_2017 + (1|Block), data=trillium_count)

Anova(tadd_glm)
plot(emmeans(tadd_glm,  ~ Treatment_2017, type="response"))

pairs(emmeans(tadd_glm,  ~ Treatment_2017, type="response"))



# ok lets find the total counts per block, yeah no difference
# check that the blocks are right. these totals for the add group make little sense
# they could be right, but just triple check




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
specpool(hbt_species)

FigS1a <- plot(hbt_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main="Henry Buck Trail")

#richness estimates (various approaches)
specpool(hbt_species)

poolaccum(hbt_species)

plot(poolaccum(hbt_species), display = "jack1")

estimateR(hbt_species) %>% plot()






# Figure S1b #####

# Galko preserve
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
specpool(gft_species)

FigS1b <- plot(gft_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main = "Galko Farm Preserve")



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

specpool(rmu_species)

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
specpool(rmd_species)

FigS1d <-  plot(rmd_curve, ci.type = c("line"), ylab="Understory Plant Richness", 
     xlab = "Meters Sampled", ylim = c(0,20), main="Ragged Mountain downslope")





# Combine Fig S1 ####
# FigS1abcd <- ggarrange(FigS1a, FigS1b, FigS1c, FigS1d, labels = c("A", "B", "C", "D"), nrow = 2)
# ggsave(filename= "./Figures/FigS1a.jpeg", plot = FigS1a , device = "jpeg", width=3, height=3, units = "in")





# New Figure S2 ######
# A pie chart for 2010 henry buck data
hbt_dat


# Create new dataframe for calculating total coverage of ant plants
pie_dat <- hbt_dat

pie_dat$ant_plants <- 
  (pie_dat$spring_beauty +
     pie_dat$trout_lily +
     pie_dat$red_trillium +
     pie_dat$d_breeches +
     pie_dat$anenome_sp +
     pie_dat$viola_sp_2 +
     pie_dat$white_violet)

# most common myrmecochores are:
# Spring beauty, trout lily, dutchman's, trillium




# Non-myrmecochores
# Create variable for total coverage
pie_dat <- pie_dat %>%
  rowwise() %>%
  mutate(total_plants = sum(across(spring_beauty:partridgeberry), na.rm = T))   %>%
  as.data.frame()

# Create variable for non-ant-dispersed plants
pie_dat$non_ant_plants <- pie_dat$total_plants - pie_dat$ant_plants

# curious what the proportion is

pie3_dat$total_plants / pie3_dat$ant_plants


# Other non-common myrmecochores
# create variable for common myrmecochores
pie_dat <- pie_dat %>%
  rowwise() %>%
  mutate(total_commons = sum(across(c("red_trillium", "trout_lily","spring_beauty","d_breeches")), na.rm = T))

# total non common ant plants
pie_dat$non_common_ant_plants <- pie_dat$ant_plants - pie_dat$total_commons




# pie chart should have 6 variables:
# Spring beauty, trout lily, dutchman's, trillium, non-myrmecochores, rare myrmecochores

pie2_dat <- subset(pie_dat, select=c("section","red_trillium", "trout_lily","spring_beauty","d_breeches", "non_common_ant_plants", "non_ant_plants"))


# then group into 3 to 9 circles based on these counts sorted by treatment

# summarize by each treatment combo
pie2_dat$Block <- as.factor(pie2_dat$section)


pie3_dat <- pie2_dat %>%
  group_by(Block) %>%
  summarize(Trillium = sum(red_trillium), 
            Trout_lily = sum(trout_lily), 
            Dutchmans_breeches = sum(d_breeches), 
            Spring_beauty = sum(spring_beauty),
            Other_myrmecohores = sum(non_common_ant_plants),
            Non_myrmecochores = sum(non_ant_plants))

# proportions (approx 15%)
1 - (pie3_dat$Trillium + pie3_dat$Trout_lily + pie3_dat$Dutchmans_breeches + pie3_dat$Spring_beauty + pie3_dat$Other_myrmecohores) /
(pie3_dat$Trillium + pie3_dat$Trout_lily + pie3_dat$Dutchmans_breeches + pie3_dat$Spring_beauty + pie3_dat$Other_myrmecohores + pie3_dat$Non_myrmecochores)

# pie 1
# jesus christ i would have been done with this in an hour if i used excel

p1 <- colnames(pie3_dat) %>% as.data.frame()
p1$cover <- t(pie3_dat[1,])
p1 <- p1[-c(1), ]
names(p1)[1] <- "plant"
p1$plant <- as.factor(p1$plant)
p1$cover <- as.numeric(p1$cover)



# Make a list of the names to use in the figure:
pie_list <- c("Dutchman's breeches", "Non-myrmecochores", "Other myrmecochores",
              "Spring beauty", "Red trillium", "Trout Lily")

pie_order <- c("Dutchmans_breeches", "Trillium", "Trout_lily", "Spring_beauty", "Other_myrmecohores", "Non-myrmecochores")

pie_list_2 <- c("Dutchman's breeches", "Red trillium", "Trout lily", "Spring beauty", 
                "Other myrmecochores", "Non-myrmecochores")

#ggplot pie chart
# this treatment should be: control
ps2_plot <- ggplot(p1, aes(x="", y=cover, fill=factor(plant, levels=pie_order))) +
  geom_bar(stat="identity", width=1, color="black") +
  # translate inches to cm for plotting the coverage values
  geom_text(aes(x=1.7, label = signif((cover*2.54), digits=2)),
            position = position_stack(vjust = 0.5), size=5) +
  coord_polar(theta="y") +
  theme_void(base_size = 24) +
  scale_fill_brewer(name = "Plant Category", labels=pie_list_2, palette = "Accent")
ps2_plot


