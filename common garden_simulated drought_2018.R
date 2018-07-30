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
pots <- read.csv('Outside_datafinal.csv')%>%
  mutate(bed=bed.., pot=plant..)%>%
  select(bed, pot)%>%
  unique()
trt <- read.csv('Outside_treatments.csv')%>%
  mutate(pot=plant..)%>%
  select(-plant..)%>%
  left_join(pots)
growth <- read.csv('Outside_datafinal.csv')%>%
  filter(height!='NA')%>%
  mutate(bed=bed.., pot=plant.., height=height, leaf_num=X..leaves, herbivorized_leaves=X..leaves.herb., percent_herbivory=avg..herb...leaf)%>%
  select(bed, pot, date, height, leaf_num, herbivorized_leaves, percent_herbivory)%>%
  left_join(trt)%>%
  #select only the last date (two weeks after rabbit herbivory)
  filter(date=='7/27/2018', height!=999, leaf_num!=999, herbivorized_leaves!=999, percent_herbivory!=999)


###repeated measures ANOVA for height
summary(heightANOVA <- aov(height ~ treatment, data=growth))

ggplot(data=barGraphStats(data=growth, variable="height", byFactorNames=c("treatment")), aes(x=as.factor(treatment), y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Treatment') +
  ylab('Height (cm)')


###repeated measures ANOVA for leaf number
summary(leavesANOVA <- aov(leaf_num ~ treatment, data=growth))

ggplot(data=barGraphStats(data=growth, variable="leaf_num", byFactorNames=c("treatment")), aes(x=as.factor(treatment), y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Treatment') +
  ylab('Leaf Number')


###repeated measures ANOVA for percent insect herbivory
summary(herbivoryANOVA <- aov(percent_herbivory ~ treatment, data=growth))

ggplot(data=barGraphStats(data=growth, variable="percent_herbivory", byFactorNames=c("treatment")), aes(x=as.factor(treatment), y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Treatment') +
  ylab('Percent Insect Herbivory')
