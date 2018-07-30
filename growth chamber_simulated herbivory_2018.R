library(tidyverse)
library(lme4)
library(lmerTest)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


###########################################################################
###########################################################################

###read in data
pots <- read.csv('Inside_datafinal.csv')%>%
  select(block, pot, treatment)%>%
  unique()
trt <- read.csv('Inside_treatments.csv')%>%
  mutate(treatment=number)%>%
  select(-number)%>%
  left_join(pots)
growth <- read.csv('Inside_datafinal.csv')%>%
  left_join(trt)
nodules <- read.csv('Inside_nodule number.csv')%>%
  left_join(trt)
biomass <- read.csv('Inside_biomass.csv')%>%
  left_join(trt)%>%
  mutate(total_biomass=shoot+root)

###repeated measures ANOVA for height
summary(heightRManova <- lmer(height ~ diversity_treatment*herbivory*date + (1+pot|date), data=subset(growth, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(growth, date=='7/23/18' & diversity_treatment!=0), variable="height", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Height (cm)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


###repeated measures ANOVA for leaf number
summary(leavesRManova <- lmer(num_leaves ~ diversity_treatment*herbivory*date + (1+pot|date), data=subset(growth, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(growth, date=='7/23/18' & diversity_treatment!=0), variable="num_leaves", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Leaf Number') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


###ANOVA for nodule number
summary(nodulesGLM <- aov(nodules ~ diversity_treatment*herbivory, data=subset(nodules, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(nodules, diversity_treatment!=0), variable="nodules", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Nodule Number') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


###ANOVA for biomass
summary(shootGLM <- aov(shoot ~ diversity_treatment*herbivory, data=subset(biomass, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="shoot", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Shoot Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


summary(rootGLM <- aov(root ~ diversity_treatment*herbivory, data=subset(biomass, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="root", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Root Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))

summary(totbioGLM <- aov(total_biomass ~ diversity_treatment*herbivory, data=subset(biomass, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="total_biomass", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Number') +
  ylab('Total Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


###monocultures only
###repeated measures ANOVA for height
summary(heightRManova <- lmer(height ~ diversity_description*herbivory*date + (1+pot|date), data=subset(growth, diversity_treatment==1)))

ggplot(data=barGraphStats(data=subset(growth, date=='7/23/18' & diversity_treatment==1), variable="height", byFactorNames=c("diversity_description", "herbivory")), aes(x=as.factor(diversity_description), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Identity') +
  ylab('Height (cm)') +
  scale_x_discrete(labels=c('110', '76', '136', '138')) +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))
#export at 800x800

###repeated measures ANOVA for leaf number
summary(leavesRManova <- lmer(num_leaves ~ diversity_description*herbivory*date + (1+pot|date), data=subset(growth, diversity_treatment==1)))

ggplot(data=barGraphStats(data=subset(growth, date=='7/23/18' & diversity_treatment==1), variable="num_leaves", byFactorNames=c("diversity_description", "herbivory")), aes(x=as.factor(diversity_description), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Identity') +
  ylab('Leaf Number') +
  scale_x_discrete(labels=c('110', '76', '136', '138')) +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))
#export at 800x800

###repeated measures ANOVA for nodule number
summary(nodulesGLM <- aov(nodules ~ diversity_description*herbivory, data=subset(nodules, diversity_treatment!=0)))

ggplot(data=barGraphStats(data=subset(nodules, diversity_treatment==1), variable="nodules", byFactorNames=c("diversity_description", "herbivory")), aes(x=as.factor(diversity_description), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Strain Identity') +
  ylab('Nodule Number') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))
