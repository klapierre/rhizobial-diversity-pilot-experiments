################################################################################
##  growth chamber_simulated herbivory_2018.R: Examining the effects of rhizobial diversity on soybean resilience following simulated herbivory.
##
##  Author: Kimberly Komatsu
##  Date created: July 29, 2018
################################################################################

library(nlme)
library(lsmeans)
library(car)
library(performance)
library(tidyverse)
library(ggpubr)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch\\final data_komatsu approved')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

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
  select(-number)%>%mutate(strains_code=ifelse(diversity_description=='1', 'A',
                                               ifelse(diversity_description=='2', 'B',
                                                      ifelse(diversity_description=='3', 'C',
                                                             ifelse(diversity_description=='4', 'D',
                                                                    ifelse(diversity_description=='1,2', 'AB',
                                                                           ifelse(diversity_description=='1,3', 'AC',
                                                                                  ifelse(diversity_description=='1,4', 'AD',
                                                                                         ifelse(diversity_description=='2,3', 'BC',
                                                                                                ifelse(diversity_description=='2,4', 'BD',
                                                                                                       ifelse(diversity_description=='3,4', 'CD',
                                                                                                              ifelse(diversity_description=='1,2,3', 'ABC',
                                                                                                                     ifelse(diversity_description=='1,2,4', 'ABD',
                                                                                                                            ifelse(diversity_description=='1,3,4', 'ACD',
                                                                                                                                   ifelse(diversity_description=='2,3,4', 'BCD', 'control')))))))))))))))%>%
  left_join(pots)
growth <- read.csv('Inside_datafinal.csv')%>%
  left_join(trt)%>%
  mutate(date=as.Date(date, format='%m/%d/%y', origin='1899-12-30'))
nodules <- read.csv('Inside_nodule number.csv')%>%
  left_join(trt)
biomass <- read.csv('Inside_biomass.csv')%>%
  left_join(trt)%>%
  mutate(total_biomass=shoot+root)

###repeated measures ANOVA for height
summary(heightModel <- lme(log10(height)~as.factor(diversity_treatment)*as.factor(herbivory)*date,
                           data=subset(growth, diversity_treatment!=0 & num_leaves!=0), 
                           random=~1|block,
                           correlation=corAR1(form=~date|block/pot),
                           control=lmeControl(returnObject=T)))
check_model(heightModel)
anova.lme(heightModel, type='marginal') #significant effect of time

heightPlot <- ggplot(data=barGraphStats(data=subset(growth, diversity_treatment!=0 & num_leaves!=0), variable="height", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Height (cm)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'),
                    name='Simulated\nHerbivory')


###repeated measures ANOVA for leaf number
summary(leavesRManova <- lme(log10(num_leaves)~as.factor(diversity_treatment)*as.factor(herbivory)*date,
                           data=subset(growth, diversity_treatment!=0 & num_leaves!=0), 
                           random=~1|block,
                           correlation=corAR1(form=~date|block/pot),
                           control=lmeControl(returnObject=T)))
check_model(leavesRManova)
anova.lme(leavesRManova, type='marginal') #significant effect of time

ggplot(data=barGraphStats(data=subset(growth, diversity_treatment!=0 & num_leaves!=0), variable="num_leaves", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Leaf Number') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))


###ANOVA for nodule number
summary(nodulesANOVA <- lme(sqrt(nodules)~as.factor(diversity_treatment)*as.factor(herbivory),
                             data=subset(nodules, diversity_treatment!=0), 
                             random=~1|block))
anova.lme(nodulesANOVA, type='marginal') #no effect

nodulePlot <- ggplot(data=barGraphStats(data=subset(nodules, diversity_treatment!=0), variable="nodules", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Nodule Number') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'),
                    name='Simulated\nHerbivory')


###ANOVA for biomass
summary(shootANOVA <- lme(shoot~as.factor(diversity_treatment)*as.factor(herbivory),
                            data=subset(biomass, diversity_treatment!=0), 
                            random=~1|block))
anova.lme(shootANOVA, type='marginal') #no effect

shootPlot <- ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="shoot", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Shoot Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'),
                    name='Simulated\nHerbivory')


summary(rootANOVA <- lme(root~as.factor(diversity_treatment)*as.factor(herbivory),
                          data=subset(biomass, diversity_treatment!=0), 
                          random=~1|block))
anova.lme(rootANOVA, type='marginal') #no effect

rootPlot <- ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="root", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Root Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'),
                    name='Simulated\nHerbivory')



summary(totbioANOVA <- lme(total_biomass~as.factor(diversity_treatment)*as.factor(herbivory),
                         data=subset(biomass, diversity_treatment!=0), 
                         random=~1|block))
anova.lme(totbioANOVA, type='marginal') #no effect

ggplot(data=barGraphStats(data=subset(biomass, diversity_treatment!=0), variable="total_biomass", byFactorNames=c("diversity_treatment", "herbivory")), aes(x=as.factor(diversity_treatment), y=mean, fill=herbivory)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity') +
  ylab('Total Biomass (g)') +
  scale_fill_manual(values=c('dark orange', 'dark green', 'dark blue'))
#export at 665 x 610



#resilience figure
ggarrange(heightPlot, shootPlot, rootPlot, nodulePlot,
          ncol = 2, nrow = 2)
#export at 800x1800