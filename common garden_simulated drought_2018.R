library(tidyverse)
library(lme4)
library(lmerTest)

#laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch')

#desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch')

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
pots <- read.csv('Outside_datafinal.csv')%>%
  mutate(bed=bed.., pot=plant..)%>%
  select(bed, pot)%>%
  unique()
trt <- read.csv('Outside_treatments.csv')%>%
  mutate(pot=plant..)%>%
  select(-plant..)%>%
  left_join(pots)%>%
  left_join(read.csv('garden_warming treatments.csv'))
growth <- read.csv('Outside_datafinal.csv')%>%
  filter(height!='NA')%>%
  mutate(bed=bed.., pot=plant.., height=height, leaf_num=X..leaves, herbivorized_leaves=X..leaves.herb., percent_herbivory=avg..herb...leaf)%>%
  select(bed, pot, date, height, leaf_num, herbivorized_leaves, percent_herbivory)%>%
  left_join(trt)%>%
  #select only the last date (two weeks after rabbit herbivory)
  filter(date=='7/27/2018', height!=999, leaf_num!=999, herbivorized_leaves!=999, percent_herbivory!=999)


###ANOVA for height
summary(heightANOVA <- aov(height ~ diversity*warming, data=growth))
summary(heightANOVA <- aov(height ~ diversity, data=subset(growth, warming==0))) #only among unwarmed
summary(heightANOVA <- aov(height ~ diversity, data=subset(growth, warming==1))) #only among unwarmed

#drought*diversity interaction
ggplot(data=barGraphStats(data=growth, variable="height", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Treatment') +
  ylab('Height (cm)')

#drought only
heightDroughtPlot <- ggplot(data=barGraphStats(data=growth, variable="height", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Warming Treatment') +
  ylab('Height (cm)') +
  scale_fill_manual(values=c("blue", 'orange')) +
  theme(legend.position='none')

#diversity only -- unwarmed plots
heightDriversityUnwarmedPlot <- ggplot(data=barGraphStats(data=subset(growth, warming==0), variable="height", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Height (cm)') +
  theme(legend.position='none')

#diversity only -- warmed plots
heightDriversityWarmedPlot <- ggplot(data=barGraphStats(data=subset(growth, warming==1), variable="height", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Height (cm)') +
  scale_fill_gradient(low="orange", high="red") +
  theme(legend.position='none')



###ANOVA for leaf number
summary(leavesANOVA <- aov(leaf_num ~ diversity*warming, data=growth))

ggplot(data=barGraphStats(data=growth, variable="leaf_num", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaf Number')


###ANOVA for percent insect herbivory
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity*warming, data=growth))
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity, data=subset(growth, warming==0))) #only among unwarmed
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity, data=subset(growth, warming==1))) #only among unwarmed

ggplot(data=barGraphStats(data=growth, variable="percent_herbivory", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Percent Insect Herbivory')

#drought only
herbivoryDroughtPlot <- ggplot(data=barGraphStats(data=growth, variable="percent_herbivory", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Warming Treatment') +
  ylab('Insect Herbivory (%)') +
  scale_fill_manual(values=c("blue", 'orange')) +
  theme(legend.position='none')

#diversity only -- unwarmed plots
herbivoryDriversityUnwarmedPlot <- ggplot(data=barGraphStats(data=subset(growth, warming==0), variable="percent_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Insect Herbivory (%)') +
  theme(legend.position='none')

#diversity only -- warmed plots
herbivoryDriversityWarmedPlot <- ggplot(data=barGraphStats(data=subset(growth, warming==1), variable="percent_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Insect Herbivory (%)') +
  scale_fill_gradient(low="orange", high="red") +
  theme(legend.position='none')
