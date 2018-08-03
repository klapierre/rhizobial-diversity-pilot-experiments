library(tidyverse)
library(lme4)
library(lmerTest)
library(grid)

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
  select(-plant.., -bed)%>%
  filter(pot!='NA')%>%
  left_join(pots)%>%
  left_join(read.csv('garden_warming treatments.csv'))
growth <- read.csv('Outside_datafinal.csv')%>%
  mutate(bed=bed.., pot=plant.., height=height, leaf_num=X..leaves, herbivorized_leaves=X..leaves.herb., percent_herbivory=avg..herb...leaf, rabbit_removed_count=rabbit.removed, date_temp=date)%>%
  mutate(date=ifelse(date_temp=='6/28/2018', 'Jun28', ifelse(date_temp=='7/11/2018', 'Jul11', ifelse(date_temp=='7/19/2018', 'Jul19', 'Jul27'))))%>%
  select(bed, pot, date, height, leaf_num, herbivorized_leaves, percent_herbivory, rabbit_removed_count)%>%
  left_join(trt)
finalDate <- growth%>%
  filter(date=='Jul27', height!=999, leaf_num!=999, herbivorized_leaves!=999, percent_herbivory!=999)
leaves <- growth[-c(323,635),]%>% #note that catastrophic rabbit herbivory occurred July 13-15
  filter(leaf_num!=999)%>%
  select(bed, pot, date, leaf_num, strains, diversity, warming)%>%
    spread(key=date, value=leaf_num)%>%
  mutate(leaves_since_rabbit=Jul27-Jul19)


###ANOVA for height
summary(heightANOVA <- aov(height ~ diversity*warming, data=finalDate))
summary(heightANOVA <- aov(height ~ diversity, data=subset(finalDate, warming==0))) #only among unwarmed
summary(heightANOVA <- aov(height ~ diversity, data=subset(finalDate, warming==1))) #only among unwarmed

#drought*diversity interaction
ggplot(data=barGraphStats(data=finalDate, variable="height", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Treatment') +
  ylab('Height (cm)')

#drought only
heightDroughtPlot <- ggplot(data=barGraphStats(data=finalDate, variable="height", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Warming Treatment') +
  ylab('Height (cm)') +
  scale_fill_manual(values=c("blue", 'orange')) +
  theme(legend.position='none')
#export at 800x800

#diversity only -- unwarmed plots
heightDriversityUnwarmedPlot <- ggplot(data=barGraphStats(data=subset(finalDate, warming==0), variable="height", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Height (cm)') +
  theme(legend.position='none')
#export at 800x800

#diversity only -- warmed plots
heightDriversityWarmedPlot <- ggplot(data=barGraphStats(data=subset(finalDate, warming==1), variable="height", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Height (cm)') +
  scale_fill_gradient(low="orange", high="red") +
  theme(legend.position='none')
#export at 800x800



###ANOVA for leaf number - July 27th date
summary(leavesANOVA <- aov(Jul27 ~ diversity*warming, data=leaves))

ggplot(data=barGraphStats(data=subset(leaves, Jul27!='NA'), variable="Jul27", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaf Number')

###ANOVA for leaf number - July 11th date
summary(leavesANOVA <- aov(Jul11 ~ diversity*warming, data=leaves))

ggplot(data=barGraphStats(data=subset(leaves, Jul11!='NA'), variable="Jul11", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaf Number')


###ANOVA for percent insect herbivory
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity*warming, data=finalDate))
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity, data=subset(finalDate, warming==0))) #only among unwarmed
summary(herbivoryANOVA <- aov(percent_herbivory ~ diversity, data=subset(finalDate, warming==1))) #only among unwarmed

ggplot(data=barGraphStats(data=finalDate, variable="percent_herbivory", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Percent Insect Herbivory') +
  scale_fill_manual(values=c("blue", "orange"))

#drought only
herbivoryDroughtPlot <- ggplot(data=barGraphStats(data=finalDate, variable="percent_herbivory", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Warming Treatment') +
  ylab('Insect Herbivory (%)') +
  scale_fill_manual(values=c("blue", 'orange')) +
  theme(legend.position='none')
#export at 800x800

#diversity only -- averaged across warming treatment
herbivoryDriversityWarmedPlot <- ggplot(data=barGraphStats(data=finalDate, variable="percent_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Insect Herbivory (%)') +
  theme(legend.position='none')

#diversity only -- unwarmed plots
herbivoryDriversityUnwarmedPlot <- ggplot(data=barGraphStats(data=subset(finalDate, warming==0), variable="percent_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('') +
  ylab('Insect Herbivory (%)') +
  theme(legend.position='none') +
  annotate('text', x=0.7, y=7, label='(a)', size=8) +
  theme(axis.title.x=element_text(size=30, vjust=-0.35), axis.text.x=element_text(size=20), axis.title.y=element_text(size=30, angle=90, vjust=0.5), axis.text.y=element_text(size=20))
#export at 800x800

#diversity only -- warmed plots
herbivoryDriversityWarmedPlot <- ggplot(data=barGraphStats(data=subset(finalDate, warming==1), variable="percent_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean, fill=diversity)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Insect Herbivory (%)') +
  scale_fill_gradient(low="orange", high="red") +
  theme(legend.position='none') +
  annotate('text', x=0.7, y=8, label='(b)', size=8) +
  theme(axis.title.x=element_text(size=30, vjust=-0.35), axis.text.x=element_text(size=20), axis.title.y=element_text(size=30, angle=90, vjust=0.5), axis.text.y=element_text(size=20))
#export at 800x800

pushViewport(viewport(layout=grid.layout(2,1)))
print(herbivoryDriversityUnwarmedPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(herbivoryDriversityWarmedPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 800x1600


###recovery from rabbit
#rabbit leaves removed
rabbitRemoval <- growth%>%
  filter(date=='Jul19')%>%
  select(bed, pot, date, rabbit_removed_count, treatment_code, strains, diversity, warming)%>%
  left_join(leaves)%>%
  #calculate number of leaves removed by subtracting July 19th data from July 11th data (because in most cases, it seems the rabbit removed multiple leaves on one stem)
  mutate(rabbit_removed_calculated=Jul11-Jul19)%>%
  select(bed, pot, rabbit_removed_count, rabbit_removed_calculated, treatment_code, strains, diversity, warming, Jul19)%>%
  gather(key=rabbit_removed_metric, value=rabbit_removed, -bed, -pot, -strains, -diversity, -warming, -Jul19)%>%
  group_by(bed, pot,strains, diversity, warming, Jul19)%>%
  #find max across the counted rabbit removal stalks vs calculated number to get an estimated measure
  summarize(rabbit_removed_est=max(rabbit_removed))%>%
  ungroup()%>%
  filter(!is.na(rabbit_removed_est))%>%
  #calculate the proportion of leaves removed by dividing by total (Jul 19 leaves plus those removed)
  mutate(total_leaves=rabbit_removed_est+Jul19, rabbit_removed_proportion=rabbit_removed_est/total_leaves)

#by number of leaves removed
summary(rabbitRemovalANOVA <- aov(rabbit_removed_est ~ diversity*warming, data=rabbitRemoval))

ggplot(data=barGraphStats(data=rabbitRemoval, variable="rabbit_removed_est", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaves Removed by Rabbit') +
  scale_fill_manual(values=c("blue", "orange"))
#export at 800x800

ggplot(data=barGraphStats(data=rabbitRemoval, variable="rabbit_removed_est", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaves Removed by Rabbit')
#export at 800x800


#by proportion of leaves removed
summary(rabbitRemovalANOVA <- aov(rabbit_removed_proportion ~ diversity*warming, data=rabbitRemoval))

ggplot(data=barGraphStats(data=rabbitRemoval, variable="rabbit_removed_proportion", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Proportion Leaves Removed by Rabbit') +
  scale_fill_manual(values=c("blue", "orange"))
#export at 800x800

ggplot(data=barGraphStats(data=rabbitRemoval, variable="rabbit_removed_proportion", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Proportion Leaves Removed by Rabbit')
#export at 800x800

#estimated total leaf number - July 11th date
summary(leavesANOVA <- aov(total_leaves ~ diversity*warming, data=rabbitRemoval))

ggplot(data=barGraphStats(data=subset(leaves, Jul11!='NA'), variable="Jul11", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaf Number')



#leaves regrown since rabbit event
summary(rabbitRecoveryANOVA <- aov(leaves_since_rabbit ~ diversity*warming, data=leaves))

ggplot(data=barGraphStats(data=subset(leaves, !is.na(leaves_since_rabbit)), variable="leaves_since_rabbit", byFactorNames=c("diversity", "warming")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  xlab('Rhizobial Diversity Treatment') +
  ylab('Leaves Grown Since Rabbit') +
  scale_fill_manual(values=c("blue", "orange"))
#export at 800x800








