### Code and data for Doropoulos & Roff (2021) 
### "Colouring coral larvae for tracking dispersal"

###### Load libraries ###### 
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(emmeans)
library(multcomp)
library(car)

###### Experiment 1. load data ###### 

# A. spathulata larvae were placed with 10 mL stain in individual scintillation vials for 1, 6, 12,
# and 24 hr incubation periods at three different stain concentrations (1, 10, 100 mg l-1). For 
# each of the vital stains, three replicates were conducted per each treatment (incubation time 
# * stain concentration), with 20 larvae assigned to each replicate (36 treatments, 720 larvae 
# total). Across all stains, this resulted in 144 treatments, totalling 2880 larvae. Staining 
# was conducted in independent glass scintillation vials (10mL total volume). Larvae were added 
# to each treatment so that the end point of the staining was the same across all treatments (i.e.
# larvae were added to the 1hr treatment at hr 23, 6hr treatment at hr 18, 12 hr treatment at hr 
# 12), at which point larvae were 6 days old. The intensity of staining for each replicate was 
# scored ordinally by a single observer (CD) into four categories: 1) no stain, 2) light staining,
# 3) medium staining, 4) strong staining. The proportion of alive larvae was counted under a 
# dissecting microscope to determine the effects of staining on larval survival. To determine the 
# effects of stain treatments on larval settlement, surviving larvae from each treatment were 
# added to an individual conditioned settlement tile in 250 mL of filtered seawater. Water changes 
# were conducted after 2 days (8 days after spawning), and larval settlement was scored 3 days 
# after tiles were introduced (9 days after spawning).

### Load Acropora spathulata data
acropora.data <- read.csv("1stExpAcropora.csv", header=T) %>% 
      filter(Species=="spathulata") %>%
      mutate(across(LarvaeAlive, .fns = ~replace_na(.,0))) %>%  
      mutate(across(LarvaeAlive, as.integer)) %>% 
      mutate(across(Incubation, as.numeric)) %>%
      mutate(IncubationF = Incubation) %>%
      mutate(across(IncubationF, as.factor)) %>%
      mutate(across(Concentration, as.factor)) %>% 
      mutate(across(Treatment, as.factor)) %>% 
      mutate(across(Stain, as.factor)) %>% 
      filter(Treatment!= "Control") %>% droplevels() %>%
      mutate(Proportion.survived = LarvaeAlive/LarvaeIn) %>%
      mutate(Proportion.settled = Settled/LarvaeIn) %>%
      mutate(Treatment = fct_relevel(Treatment, c("NeutralRed",  "NileBlue", "AlizarinRed", "CalceinBlue"))) %>%
      unite("StainDye",c("Treatment","Stain"),remove = FALSE) %>%
      mutate(across(StainDye, as.factor)) %>%
      mutate(StainDye = fct_relevel(StainDye, c("NeutralRed_None",  "NeutralRed_Light", "NeutralRed_Medium", "NeutralRed_Heavy", "NileBlue_Light", "NileBlue_Medium", "NileBlue_Heavy", "CalceinBlue_None", "AlizarinRed_None", "AlizarinRed_Light"))) %>%
      unite("ConcIncF",c("Concentration", "IncubationF"),remove = FALSE) %>%
      mutate(across(ConcIncF, as.factor))

# Calculate mean settlement (Acropora spathulata)
acropora.mean.data <- acropora.data %>% 
      group_by(Treatment, Incubation, Concentration, Stain, StainDye) %>% 
      summarise(Proportion.survived.mean = mean(Proportion.survived), Proportion.survived.se = sd(Proportion.survived)/sqrt(n()), Proportion.settled.mean = mean(Proportion.settled), Proportion.settled.se = sd(Proportion.settled)/sqrt(n()),n_samples = n()) %>% 
      mutate(Stain = factor(Stain, levels=c("None","Light","Medium","Heavy"))) %>%
      mutate(Treatment = fct_relevel(Treatment, c("NeutralRed",  "NileBlue", "AlizarinRed", "CalceinBlue")))# %>%

# Re-import Acropora spathulata dataset and filter to only timepoint 12 to compare larval survival 
# and larval settlement for the four stains (neutral red, Nile blue, alizarin red, calcein blue) 
# at different concentration levels after 12 hours of incubation and control (unstained) larvae (Figure S1).
# Main analysis excludes control [ filter(Treatment!= "Control") ] as unstained larvae were only
# compared at time-point 12 in the experiment

acropora.data.control <- read.csv("1stExpAcropora.csv", header=T) %>% 
  filter(Species=="spathulata") %>%
  mutate(across(LarvaeAlive, .fns = ~replace_na(.,0))) %>%  
  mutate(across(LarvaeAlive, as.integer)) %>% 
  mutate(across(Incubation, as.numeric)) %>%
  mutate(IncubationF = as.factor(Incubation)) %>%
  mutate(across(Concentration, as.factor)) %>% 
  mutate(across(Treatment, as.factor)) %>% 
  mutate(across(Stain, as.factor)) %>% 
  #filter(Treatment!= "Control") %>% droplevels() %>%
  mutate(Proportion.survived = LarvaeAlive/LarvaeIn) %>%
  mutate(Proportion.settled = Settled/LarvaeIn) %>%
  unite("StainDye",c("Treatment","Stain"),remove = FALSE) %>%
  mutate(across(StainDye, as.factor)) %>%
  mutate(StainDye = fct_relevel(StainDye, c("NeutralRed_None",  "NeutralRed_Light", "NeutralRed_Medium", "NeutralRed_Heavy", "NileBlue_Light", "NileBlue_Medium", "NileBlue_Heavy", "CalceinBlue_None", "AlizarinRed_None", "AlizarinRed_Light"))) %>%
  unite("alltreatments",c("Treatment", "Concentration", "IncubationF"),remove = FALSE) %>%
  mutate(across(alltreatments, as.factor)) %>%
  filter(IncubationF=="12")

# Calculate mean settlement across treatments at timepoint 12 (Acropora spathulata)
acropora.mean.data.control <- acropora.data.control %>%  
  group_by(Treatment, Incubation, Concentration, Stain, StainDye, alltreatments) %>% 
  summarise(Proportion.survived.mean = mean(Proportion.survived), Proportion.survived.se = sd(Proportion.survived)/sqrt(n()), Proportion.settled.mean = mean(Proportion.settled), Proportion.settled.se = sd(Proportion.settled)/sqrt(n()),n_samples = n()) %>% 
  mutate(Stain = factor(Stain, levels=c("None","Light","Medium","Heavy"))) %>%
  mutate(Treatment = fct_relevel(Treatment, c("Control","AlizarinRed", "CalceinBlue",  "NileBlue", "NeutralRed")))
  
### Load Platygyra daedalea survival data
platygyra.data.survival <- read.csv("1stExpPlatygyra.csv", header=T) %>% # 
      filter(Lifestage=="Survival") %>%
      mutate(across(LarvaeAlive, .fns = ~replace_na(.,0))) %>%  # na
      mutate(across(LarvaeAlive, as.integer)) %>% 
      mutate(across(Incubation, as.numeric)) %>%
      mutate(across(Concentration, as.factor)) %>% 
      mutate(across(Treatment, as.factor)) %>% 
      mutate(across(Stain, as.factor)) %>% 
      filter(Treatment!= "Control") %>% droplevels() %>%
      mutate(Proportion.survived = LarvaeAlive/LarvaeIn) %>%
      mutate(Treatment = fct_relevel(Treatment, c("NeutralRed",  "NileBlue", "AlizarinRed", "CalceinBlue"))) %>%
      unite("StainDye",c("Treatment","Stain"),remove = FALSE) %>%
      mutate(across(StainDye, as.factor)) %>%
      mutate(StainDye = fct_relevel(StainDye, c("NeutralRed_Heavy",  "NileBlue_Heavy", "CalceinBlue_None", "AlizarinRed_None")))

### Load Platygyra daedalea settlement data
platygyra.data.settlment <- read.csv("1stExpPlatygyra.csv", header=T) %>% 
      filter(Lifestage=="Settlement") %>%
      mutate(across(LarvaeAlive, .fns = ~replace_na(.,0))) %>%  # na
      mutate(across(LarvaeAlive, as.integer)) %>% 
      mutate(across(Incubation, as.numeric)) %>%
      mutate(across(Concentration, as.factor)) %>% 
      mutate(across(Treatment, as.factor)) %>% 
      mutate(across(Stain, as.factor)) %>% 
      filter(Treatment!= "Control") %>% droplevels() %>%
      mutate(Proportion.settled = Settled/LarvaeIn) %>%
      mutate(Treatment = fct_relevel(Treatment, c("NeutralRed",  "NileBlue", "AlizarinRed", "CalceinBlue"))) %>%
      unite("StainDye",c("Treatment","Stain"),remove = FALSE) %>%
      mutate(across(StainDye, as.factor)) %>%
      mutate(StainDye = fct_relevel(StainDye, c("NeutralRed_Heavy",  "NileBlue_Heavy", "CalceinBlue_None", "AlizarinRed_None")))
  
# Calculate mean settlement for Platygyra daedalea
platygyra.data.mean.settlment <- platygyra.data.settlment %>%
      group_by(Treatment, Incubation, Concentration, Stain, StainDye) %>% 
      summarise(Proportion.settled.mean = mean(Proportion.settled), Proportion.settled.se = sd(Proportion.settled)/sqrt(n()), Proportion.settled.mean = mean(Proportion.settled), Proportion.settled.se = sd(Proportion.settled)/sqrt(n()),n_samples = n()) %>% 
      mutate(Stain = factor(Stain, levels=c("None","Light","Medium","Heavy"))) %>%
      mutate(Treatment = fct_relevel(Treatment, c("NeutralRed",  "NileBlue", "AlizarinRed", "CalceinBlue")))# %>%

###### Experiment 1. plot results ###### 

# Figure 2
# Probability of larval survival for Acropora spathulata exposed to four stains at different 
# concentration levels and incubation times
ggplot() + theme_bw() + facet_wrap(Treatment~Concentration, ncol=3) + 
      geom_bar(data = acropora.mean.data, aes(x = as.factor(Incubation), y = Proportion.survived.mean, group = Concentration, fill=StainDye), position ='dodge', stat = "identity",colour = 'black', size = 0.3, width=0.9) + # add fill later
      geom_errorbar(data = acropora.mean.data, aes(x = as.factor(Incubation), ymin = Proportion.survived.mean - Proportion.survived.se, ymax = Proportion.survived.mean + Proportion.survived.se, group = Incubation), size = 0.3, width = 0.0, position = position_dodge(0.9)) +
      xlab("Time") + ylab("Probability of survival") + ggtitle("") +
      scale_y_continuous(limits = c(0, 1.4), breaks = seq(0, 1.15, by = 0.2)) + theme(legend.position = "none") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) 

ggsave("plots/Figure 2 Acropora Survival.pdf", last_plot(), width=5, height=8) 

# Probability of larval settlement Acropora spathulata exposed to four stains at different 
# concentration levels and incubation times
ggplot() + facet_wrap(Treatment ~ Concentration + Incubation, ncol=12) + coord_polar("y") +
      geom_bar(data=acropora.mean.data, aes(x="", y = 1, group=Incubation), fill="#F4F4F5", width = 1, stat="identity", colour="#000000") +
      geom_bar(data=acropora.mean.data, aes(x="", y = Proportion.settled.mean, group=Incubation, fill=StainDye), width = 1, stat="identity", colour="#000000") +
      scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) +
      theme_bw() + theme(legend.position = "none") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),  axis.text = element_blank(), panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.border = element_blank()) 
  
ggsave("plots/Figure 2 Acropora settlement.pdf", last_plot(), width=5, height=8) 

# Probability of larval survival for Platygyra daedalea exposed to four stains at different 
# concentration levels and incubation times
ggplot() + theme_bw() + facet_wrap(Treatment~Concentration, ncol=2) + 
  geom_bar(data = platygyra.data.survival, aes(x = as.factor(Incubation), y = Proportion.survived, group = Concentration, fill=StainDye), position ='dodge', stat = "identity",colour = 'black', size = 0.3, width=0.9) + # add fill later
  xlab("Time") + ylab("Probability of survival") + ggtitle("") +
  scale_y_continuous(limits = c(0, 1.4), breaks = seq(0, 1.15, by = 0.2)) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) 

ggsave("plots/Figure 2 Platygra Survival.pdf", last_plot(), width=3, height=8) 

# Probability of larval settlement for Platygyra daedalea exposed to four stains at different 
# concentration levels and incubation times
ggplot() + facet_wrap(Treatment ~ Concentration, ncol=2) + coord_polar("y") +
  geom_bar(data=platygyra.data.mean.settlment, aes(x="", y = 1, group=Concentration), fill="#F4F4F5", width = 1, stat="identity", colour="#000000") +
  geom_bar(data=platygyra.data.mean.settlment, aes(x="", y = Proportion.settled.mean, group=Treatment, fill=StainDye), width = 1, stat="identity", colour="#000000") +
  scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) +
  theme_bw() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),  axis.text = element_blank(), panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(), panel.border = element_blank()) 

ggsave("plots/Figure 2 Platygra settlement.pdf", last_plot(), width=3, height=8) 

# Figure S2
# Probability of larval survival for Acropora spathulata exposed to four stains 
# (neutral red, Nile blue, alizarin red, calcein blue) at different concentration levels after 12 hours of 
# incubation and control (unstained) larvae. 

ggplot() + theme_bw() + facet_wrap(~Treatment, ncol=5, scales = "free_x") + 
  geom_bar(data = acropora.mean.data.control, aes(x = Concentration, y = Proportion.survived.mean, group = Concentration, fill=StainDye), position ='dodge', stat = "identity",colour = 'black', size = 0.3, width=0.9) + # add fill later
  geom_errorbar(data = acropora.mean.data.control, aes(x = Concentration, ymin = Proportion.survived.mean - Proportion.survived.se, ymax = Proportion.survived.mean + Proportion.survived.se, group = Incubation), size = 0.3, width = 0.0, position = position_dodge(0.9)) +
  xlab("Concentration") + ylab("Probability of survival") + ggtitle("") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.15, by = 0.2)) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) 

ggsave("plots/Figure S1 Control survival.pdf", last_plot(), width=8, height=4) 

# Probability of larval settlement for Acropora spathulata exposed to four stains 
# (neutral red, Nile blue, alizarin red, calcein blue) at different concentration levels after 12 hours of 
# incubation and control (unstained) larvae. 
ggplot() + theme_bw() + facet_wrap(~Treatment, ncol=5, scales = "free_x") + 
  geom_bar(data = acropora.mean.data.control, aes(x = Concentration, y = Proportion.settled.mean, group = Concentration, fill=StainDye), position ='dodge', stat = "identity",colour = 'black', size = 0.3, width=0.9) + # add fill later
  geom_errorbar(data = acropora.mean.data.control, aes(x = Concentration, ymin = Proportion.settled.mean - Proportion.settled.se, ymax = Proportion.settled.mean + Proportion.settled.se, group = Incubation), size = 0.3, width = 0.0, position = position_dodge(0.9)) +
  xlab("Concentration") + ylab("Probability of settlement") + ggtitle("") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.15, by = 0.2)) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = c("NeutralRed_None"= "#ffffff", "NeutralRed_Light"= "#f69799","NeutralRed_Medium"="#cc3433","NeutralRed_Heavy"="#650e0d", "NileBlue_Light" = "#9ec9eb", "NileBlue_Medium" = "#3c5fb5", "NileBlue_Heavy"= "#04267a", "AlizarinRed_None" = "#ffffff", "AlizarinRed_Light" = "#ffe3e4", "CalceinBlue_None" = "#ffffff")) 

ggsave("plots/Figure S1 Control settlement.pdf", last_plot(), width=8, height=4) 


###### Experiment 1. logistic models ###### 

### Acropora spathulata larval survival models

# all stains
fm1 <- glm(cbind((LarvaeIn-LarvaeDead),LarvaeIn) ~  Concentration*IncubationF*Treatment, weights=LarvaeIn, data = acropora.data, family=binomial)
  Anova(fm1) # full model

# Neutral Red
fm1.neutral.red <- glm(cbind((LarvaeIn-LarvaeDead),LarvaeIn) ~  Concentration*IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="NeutralRed"), family=binomial)
    Anova(fm1.neutral.red) # significant interaction, combine fixed factors for posthoc

fm1.neutral.red.a <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~ ConcIncF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="NeutralRed"), family=binomial)
    Anova(fm1.neutral.red) 
    cld(glht(fm1.neutral.red.a, linfct=mcp(ConcIncF="Tukey")))

#Nile Blue
fm1.nile.blue <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~  Concentration*IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="NileBlue"), family=binomial)
    Anova(fm1.nile.blue) # ns interaction, concentration ns

fm1.nile.blue.a <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~  IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="NileBlue"), family=binomial)
    Anova(fm1.nile.blue.a)
    cld(glht(fm1.nile.blue.a, linfct=mcp(IncubationF="Tukey")))

# Alizarin Red
fm1.alizarin.red <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~  Concentration*IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="AlizarinRed"), family=binomial)
    Anova(fm1.alizarin.red) # ns interaction, concentration ns

fm1.alizarin.red.a <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~  IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="AlizarinRed"), family=binomial)
    Anova(fm1.alizarin.red.a)
    cld(glht(fm1.alizarin.red, linfct=mcp(IncubationF="Tukey")))

# Calcein Blue
fm1.calcein.blue <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~  Concentration*IncubationF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="CalceinBlue"), family=binomial)
    Anova(fm1.calcein.blue) # significant interaction, combine fixed factors for posthoc

fm1.calcein.blue.a <- glm(cbind(LarvaeIn-LarvaeDead,LarvaeIn) ~ ConcIncF, weights=LarvaeIn, data = acropora.data %>% filter(Treatment=="CalceinBlue") , family=binomial)
    Anova(fm1.calcein.blue.a)
    cld(glht(fm1.calcein.blue.a, linfct=mcp(ConcIncF="Tukey")))


# Compare survival of stained larval treatments with controls (available at 12hr time point only)  
fm1.control.survival <- glm(cbind((LarvaeIn-LarvaeDead),LarvaeIn) ~  alltreatments, weights=LarvaeIn, data = acropora.data.control, family=binomial)
Anova(fm1.control.survival) # full model
fm1.control.survival.pairwise <- glht(fm1.control.survival, linfct=mcp(alltreatments="Tukey")) 
summary(fm1.control.survival.pairwise)


# Compare settlement of stained larval treatments with controls (available at 12hr time point only)  
fm1.control.settlement <- glm(cbind((LarvaeIn-LarvaeDead),LarvaeIn) ~  alltreatments, weights=LarvaeIn, data = acropora.data.control, family=binomial)
Anova(fm1.control.settlement) # full model
fm1.control.settlement.pairwise <- glht(fm1.control.settlement, linfct=mcp(alltreatments="Tukey")) 
summary(fm1.control.settlement.pairwise)


### Acropora spathulata larval settlement model

# all stains
fm1.settlement <- glm(cbind(Settled,LarvaeAlive) ~  Concentration*IncubationF*Treatment, weights=LarvaeAlive, data = acropora.data, family=binomial)
Anova(fm1.settlement) # significant interaction, combine fixed factors for posthoc

###### Experiment 2. load data ###### 

# Load all data
larvae.data <- read.csv("2ndExp.csv", header=T) %>%
  mutate(PropAlive = Dead/LarvaeIn) %>%
  mutate(across(Concentration, as.factor)) %>% 
  mutate(across(Species, as.factor)) %>% 
  mutate(across(Incubation, as.factor)) %>% 
  mutate(across(Treatment, as.factor)) %>%
  mutate(across(Stain, as.factor)) %>% 
  #filter(Treatment!= "Control") %>% droplevels() %>%
  unite("TreatmentConcentration",c("Treatment","Concentration"),remove = FALSE) %>%
  mutate(across(TreatmentConcentration, as.factor)) %>%
  unite("IncubationConcentration",c("Concentration","Incubation"),remove = FALSE) %>%
  mutate(across(IncubationConcentration, as.factor)) %>%
  unite("IncubationConcentrationSpecies",c("Species","Concentration","Incubation"),remove = FALSE) %>%
  mutate(across(IncubationConcentrationSpecies, as.factor)) %>%
  mutate(Proportion.survived = (LarvaeIn-Dead)/LarvaeIn) %>%
  mutate(Proportion.settled = Settled/LarvaeIn) 
  
# Calculate mean settlement across species and treamtents  
larvae.data.mean <- larvae.data %>% 
  group_by(Treatment, Incubation, Concentration, Species, IncubationConcentration) %>% 
  summarise(Proportion.survived.mean = mean(Proportion.survived), Proportion.survived.se = sd(Proportion.survived)/sqrt(n()),
            Proportion.settled.mean = mean(Proportion.settled), Proportion.settled.se = sd(Proportion.settled)/sqrt(n()),
            n_samples = n()) %>% as.data.frame() 


###### Experiment 2.  plot results ###### 

### Figure 3
# Probability of larval survival for four species of coral exposed to neutral red 
# and Nile blue stains at different concentration levels and incubation times. 

{Figure3 <- ggplot() + theme_bw() + facet_wrap(Treatment~Species, ncol=4,scales = "free_x") + 
    geom_bar(data = larvae.data.mean, aes(x = as.factor(IncubationConcentration) , y = Proportion.survived.mean, group = Concentration, fill=Treatment), alpha=0.95, position ='dodge', stat = "identity",colour = 'black', size = 0.3, width=0.9) + # add fill later
#    geom_errorbar(data = larvae.data.mean, aes(x = as.factor(IncubationConcentration), ymin = Proportion.survived.mean - Proportion.survived.se, ymax = Proportion.survived.mean + Proportion.survived.se, group = Incubation), size = 0.3, width = 0.0, position = position_dodge(0.9)) +
    xlab("Time") + ylab("Probability of survival") + ggtitle("") +
    scale_y_continuous(limits=c(0,1.4), breaks = seq(0, 1, by = 0.2)) + theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = c("NeutralRed"="#650e0d","NileBlue"= "#04267a")) 
  
  
  Figure3 <-  ggplot_gtable(ggplot_build(Figure3))
  #Figure3$widths
  Figure3$widths[5] = 2.2*Figure3$widths[5]
  #gtable::gtable_show_layout(Figure3)
  grid::grid.draw(Figure3)
  ggsave("plots/Figure 3 Larval Survival.pdf", Figure3, width=8, height=8)   
            }

# Probability of larval survival for four species of coral exposed to neutral red 
# and Nile blue stains at different concentration levels and incubation times. 
ggplot() + facet_wrap(Treatment~Species+IncubationConcentration, ncol=11) + coord_polar("y") +
  geom_bar(data=larvae.data.mean %>% filter(!Treatment=="Control"), aes(x="", y = 1, group=IncubationConcentration), fill="#d3d3d3", width = 1, stat="identity", colour="#000000") +
  geom_bar(data=larvae.data.mean %>% filter(!Treatment=="Control"), aes(x="", y = Proportion.settled.mean, group=IncubationConcentration, fill=Treatment), width = 1, stat="identity", colour="#000000") +
  scale_fill_manual(values = c("NeutralRed"="#650e0d","NileBlue"= "#04267a")) +
  theme_bw() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(),  axis.text = element_blank(), panel.grid=element_blank(), panel.background = element_blank(), axis.line = element_blank(), strip.background = element_blank(), panel.border = element_blank()) 

ggsave("plots/Figure 3 Larval Settlement.pdf", last_plot(), width=5, height=8) 

###### Experiment 2.  logistic models  ###### 

### larval survival model
fm2.survival <- glm(cbind(LarvaeIn-Dead,LarvaeIn) ~  IncubationConcentrationSpecies, weights=LarvaeIn, data = larvae.data %>% filter(!Treatment=="Control"), family=binomial) # no significant difference in survival among species or treatments (Neutral/NileBlue)

Anova(fm2.survival)

### larval settlement model
fm2.settlement <- glm(cbind(Settled,LarvaeIn) ~  IncubationConcentrationSpecies, weights=LarvaeIn, data = larvae.data, family=binomial) #Significant difference in settlement between species*treatment interaction
Anova(fm2.settlement)

# Extract posthoc contrasts for differences in settlement between treatments and control: 
# Acropora anthocercis neutral red
Anova(fm2.settlement.a <- glm(cbind(Settled,LarvaeIn) ~  IncubationConcentration, weights=LarvaeIn, data = larvae.data %>% filter(Species=="anthocersis") %>% filter(Treatment=="NeutralRed" | Treatment=="Control"), family=binomial))
emmeans(fm2.settlement.a, pairwise ~ IncubationConcentration, type="response")

# Platygyra sinensis neutral red
Anova(fm4a <- glm(cbind(Settled,LarvaeIn) ~  IncubationConcentration, weights=LarvaeIn, data = larvae.data %>% filter(Species=="sinensis") %>% filter(Treatment=="NeutralRed" | Treatment=="Control"), family=binomial))
emmeans(fm4a, pairwise ~ IncubationConcentration, type="response")

# Acropora anthocercis nile blue
Anova(fm4a <- glm(cbind(Settled,LarvaeIn) ~  IncubationConcentration, weights=LarvaeIn, data = larvae.data %>% filter(Species=="anthocersis") %>% filter(Treatment=="NileBlue" | Treatment=="Control"), family=binomial))
emmeans(fm4a, pairwise ~ IncubationConcentration, type="response")

# Platygyra sinensis nile blue
Anova(fm4a <- glm(cbind(Settled,LarvaeIn) ~  IncubationConcentration, weights=LarvaeIn, data = larvae.data %>% filter(Species=="sinensis") %>% filter(Treatment=="NileBlue" | Treatment=="Control"), family=binomial))
emmeans(fm4a, pairwise ~ IncubationConcentration, type="response")
