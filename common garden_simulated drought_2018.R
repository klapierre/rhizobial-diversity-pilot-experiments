################################################################################
##  common garden_simulated drought_2018.R: Examining the effects of rhizobial diversity and simulated drought on resistance to herbivory and growth rates.
##
##  Author: Kimberly Komatsu
##  Date created: July 29, 2018
################################################################################

library(nlme)
library(lsmeans)
library(ggeffects)
library(sjPlot)
library(grid)
library(performance)
library(tidyverse)

#laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch\\final data_komatsu approved')

#desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\SERC\\interns\\2018\\2018_REU_Esch\\final data_komatsu approved')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=20),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=20),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=16), legend.text=element_text(size=16))

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

#not in function
`%!in%` = Negate(`%in%`)

###########################################################################
###########################################################################

###read in data
trt <- read.csv('soy pilot_2018_common garden_treatments.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))
  
growthData <- read.csv('soy pilot_2018_common garden_growth.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  left_join(trt)%>%
  select(-num_flowers, -num_pods)%>%
  filter(height_cm<9000, num_leaves<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  rename(num_leaflets=num_leaves)%>% #leaflets were counted, divide by 3 to get num leaves
  mutate(num_leaves=num_leaflets/3)%>%
  select(-num_rabbit_herb)%>% #drop num_rabbit_herb, which is number of stems removed by rabbit (counted as clipped stems, which is not very accurate if the rabbit clipped a branch with many higher stems); will calculate num leaves removed by rabbit later on
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))%>% #gets day of year (doy) and day of experiment (doe)
  mutate(height_scaled=(height_cm-mean(height_cm))/sd(height_cm))

growthRate <- growthData%>%
  mutate(doy_cat=paste('doy', doy, sep='_'))%>%
  select(doy_cat, bed, warming, plant, treatment_code, strains, diversity, USDA_110, USDA_76, USDA_136, USDA_138, height_cm)%>%
  spread(key=doy_cat, value=height_cm)%>%
  mutate(rate_1=(doy_186-doy_179)/(186-179), rate_2=(doy_192-doy_186)/(192-186), rate_3=(doy_200-doy_192)/(200-192), rate_4=(doy_208-doy_200)/(208-200), rate_5=(doy_214-doy_208)/(214-208), rate_6=(doy_221-doy_214)/(221-214), rate_7=(doy_235-doy_221)/(235-221), rate_8=(doy_253-doy_235)/(253-235))%>%
  select(bed, warming, plant, treatment_code, strains, diversity, USDA_110, USDA_76, USDA_136, USDA_138, rate_1, rate_2, rate_3, rate_4, rate_5, rate_6, rate_7, rate_8)%>%
  gather(key=rate_period, value=growth_rate, rate_1:rate_8)

herbivoryData <- read.csv('soy pilot_2018_common garden_herbivory.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  filter(perc_herbivory<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  #calculate avg chewing insect herbivory per leaf
  group_by(bed, plant, date)%>%
  summarise(sum_perc_herbivory=sum(perc_herbivory), mean_perc_herbivory=mean(perc_herbivory))%>%
  ungroup()%>%
  left_join(growthData)%>%
  mutate(avg_perc_herbivory=sum_perc_herbivory/num_leaflets)%>%
  filter(!is.na(avg_perc_herbivory), avg_perc_herbivory<100000)%>%
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))

insectData <- read.csv('soy pilot_2018_common garden_insects.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  left_join(trt)%>%
  filter(Aphids<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))

flowerData <- read.csv('soy pilot_2018_common garden_growth.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  left_join(trt)%>%
  select(-height_cm, -num_leaves)%>%
  filter(num_flowers<9000, num_pods<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))

weightData <- read.csv('soy pilot_2018_common garden_bean weight.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))

fitnessData <- read.csv('soy pilot_2018_common garden_pods.csv')%>%
  select(-comments)%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  left_join(weightData)%>%
  select(-notes)%>%
  left_join(trt)%>%
  filter(total_pods<9000)%>%
  filter(aborted_pods<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  mutate(viable_pods=total_pods-aborted_pods)%>%
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))


###height--------- for all growth analyses, don't use doy=200, which is the measurement immediately following rabbit herbivory
##height absolute
#data check
ggplot(data=growthData, aes(x=doe, y=height_cm)) +
  geom_point() +
  facet_wrap(~plant)
hist(log10(growthData$height_cm))

#model -- repeated measures to examine the effects of treatments on growth rate (slope of height through time, fixed intercept)
summary(heightModel <- lme(log10(height_cm)~as.factor(diversity)*as.factor(warming)*doe,
                           data=subset(growthData, doy!=200&diversity!=0), 
                           random=~1|bed,
                           correlation=corAR1(form=~doe|bed/plant),
                           control=lmeControl(returnObject=T)))
anova(heightModel) #significant effects of warming and day of expt
check_model(heightModel)


#height figure
ggplot(data=barGraphStats(data=subset(growthData, doy!=200&diversity!=0), variable="height_cm", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(breaks=c(0,1), labels=c("Control", "Drought")) +
  xlab('') + ylab('Height (cm)') +
  annotate('text', x=1, y=33, label='a', size=8) +
  annotate('text', x=2, y=36, label='b', size=8)
#export at 800 x 800


# ##growth rate (cm/day)---------
# #data check
# ggplot(data=growthRate, aes(x=rate_period, y=growth_rate)) +
#   geom_point() +
#   facet_wrap(~plant)
# 
# summary(growthModel <- lme(growth_rate~as.factor(diversity)*as.factor(warming)*rate_period, data=subset(growthRate, rate_period!='rate_3'&diversity!=0), 
#                            random=~1|bed, 
#                            correlation=corAR1(form=~rate_period|bed/plant),
#                            control=lmeControl(returnObject=T), weights=varIdent(form=~1|diversity*warming)))
# anova(growthModel) #significant effects of warming and day of expt
# check_model(growthModel)
# 
#           
#           lmer(growth_rate~diversity*warming + (1|bed/plant), data=subset(growthRate, rate_period!='rate_3'))) #no effects
# check_model(growthModel)


###leaf number---------
#data check
ggplot(data=growthData, aes(x=doy, y=num_leaves)) +
  geom_point() +
  facet_wrap(~plant)

#model
summary(leafModel <- lme(log10(num_leaves)~as.factor(diversity)*as.factor(warming)*doe,
                         data=subset(growthData, doy!=200&diversity!=0&num_leaves>0), 
                         random=~1|bed, 
                         correlation=corAR1(form=~doe|bed/plant),
                         control=lmeControl(returnObject=T)))
anova(leafModel) #significant interaction between warming and day of expt
check_model(leafModel)

#leaf num figure
ggplot(barGraphStats(data=subset(growthData, doy!=200&diversity!=0&num_leaves>0), variable="num_leaves", byFactorNames=c("doe", "warming")), aes(x=doe, y=mean, color=as.factor(warming))) +
  geom_point(position=position_dodge(0.9), size=4) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_color_manual(values=c('#ffc425', '#f37735'), breaks=c(0,1), labels=c('Control', 'Drought'), name='Drought\nTreatment') +
  xlab('Day of Experiment') + ylab('Number of Leaves') +
  scale_y_continuous(trans='log10') + 
  theme(legend.position=c(0.1, 0.85))
#export at 800 x 600



###insect herbivory---------
#data check
ggplot(data=subset(herbivoryData, doy!=200), aes(x=doy, y=avg_perc_herbivory)) +
  geom_point() +
  facet_wrap(~plant)

#model
summary(insectherbModel <- lme(sqrt(avg_perc_herbivory)~as.factor(diversity)*as.factor(warming)*doe,
                               data=subset(herbivoryData, doy!=200&diversity!=0), 
                               random=~1|bed, 
                               correlation=corAR1(form=~doe|bed/plant),
                               control=lmeControl(returnObject=T)))
check_model(insectherbModel)
anova(insectherbModel) #significant effect of diversity and day of expt
lsmeans(insectherbModel, pairwise~as.factor(diversity), adjust="tukey")

#insect damage figure
insectFig <- ggplot(data=barGraphStats(data=subset(herbivoryData, doy!=200 & diversity!=0), variable="avg_perc_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  xlab('') + ylab('Insect Damage (%)') +
  annotate('text', x=1, y=3, label='a', size=8) +
  annotate('text', x=2, y=2.7, label='ab', size=8) +
  annotate('text', x=3, y=2.4, label='b', size=8) +
  annotate('text', x=0.5, y=3, label='(a)', size=8)
#export at 800 x 800


###rabbit herbivory---------
leavesPreRabbit <- growthData%>%
  filter(date=='7/11/2018')%>%
  rename(num_leaves_pre=num_leaves)%>%
  select(bed, warming, plant, num_leaves_pre)
rabbitData <- growthData%>%
  filter(date=='7/19/2018')%>%
  left_join(leavesPreRabbit)%>%
  mutate(percent_rabbit_herb=(num_leaves_pre-num_leaves)/num_leaves_pre*100)

summary(rabbitModel <- lme(percent_rabbit_herb~as.factor(diversity)*as.factor(warming),
            data=subset(rabbitData, diversity!=0),
            random=~1|bed))
anova(rabbitModel) #no effect

summary(lme(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==0), random=~1|bed)) #marginally sig diversity 
summary(lme(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==1), random=~1|bed)) #no effect

rabbitFig <- ggplot(data=barGraphStats(data=subset(rabbitData, date=='7/19/2018'&diversity!=0), variable="percent_rabbit_herb", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  scale_x_discrete(breaks=c(1,2,3), labels=c("1 strain", "2 strains", "3 strains")) +
  xlab('Rhizobial Diversity') + ylab('Rabbit Damage') +
  annotate('text', x=0.5, y=45, label='(b)', size=8)
#export at 800 x 800

#herbivory figure
pushViewport(viewport(layout=grid.layout(2,1)))
print(insectFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(rabbitFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 800x1200



###flower number---------
summary(flowerModel <- lme(num_flowers~as.factor(diversity)*as.factor(warming), data=subset(flowerData, num_flowers<9000&date=='8/23/2018'&diversity!=0), random=~1|bed))

anova(flowerModel) #warming effect

flowerFig <- ggplot(data=barGraphStats(data=subset(flowerData, date=='8/23/2018'&diversity!=0), variable="num_flowers", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  scale_x_discrete(breaks=c(0,1), labels=c('Control', 'Drought')) +
  xlab('') + ylab('Flower Number') +
  theme(axis.title.y=element_text(vjust=1)) +
  annotate('text', x=1, y=45, label='a', size=8) +
  annotate('text', x=2, y=25, label='b', size=8) +
  annotate('text', x=0.5, y=45, label='(a)', size=8)

flowerFig2 <- ggplot(data=barGraphStats(data=subset(flowerData, date=='8/23/2018'&diversity!=0), variable="num_flowers", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                    breaks=c(1,2,3),
                    labels=c('1 strain', '2 strains', '3 strains'),
                    name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Flower Number') +
  theme(axis.title.y=element_text(vjust=1), legend.position=c(0.98,0.98), legend.justification=c(1,1)) +
  annotate('text', x=0.7, y=44, label='a', size=8) +
  annotate('text', x=1, y=44, label='a', size=8) +
  annotate('text', x=1.3, y=50, label='a', size=8) +
  annotate('text', x=1.7, y=27, label='b', size=8) +
  annotate('text', x=2, y=24, label='b', size=8) +
  annotate('text', x=2.3, y=26, label='b', size=8) +
  annotate('text', x=0.5, y=50, label='(a)', size=8)


###pods---------
# #total pods
# summary(lme(total_pods~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, total_pods<9000), random=~1|bed)) #no effect
# summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==0))) #no effect
# summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==1))) #no effect
# 
# ggplot(data=barGraphStats(data=subset(fitnessData, total_pods<9000&diversity!=0), variable="total_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
#   geom_point(size=5, position=position_dodge(width=0.25)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
#   scale_color_manual(values=c('#0072B2', '#D55E00'),
#                      breaks=c(0,1),
#                      labels=c('Control', 'Drought'),
#                      name='Drought\nTreatment') +
#   xlab('Rhizobial Diversity') + ylab('Total Pod Number') +
#   theme(axis.title.y=element_text(vjust=1))
# #export at 665 x 610


#viable pods---------
summary(viableModel <- lme(viable_pods~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, viable_pods<9000&diversity!=0), random=~1|bed))
anova(viableModel) #no effect

summary(glm(viable_pods~diversity, data=subset(fitnessData, viable_pods<9000&warming==0))) #no effect
summary(glm(viable_pods~diversity, data=subset(fitnessData, viable_pods<9000&warming==1))) #no effect

viableFig <- ggplot(data=barGraphStats(data=subset(fitnessData, viable_pods>0&diversity!=0), variable="viable_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                    breaks=c(1,2,3),
                    labels=c('1 strain', '2 strains', '3 strains'),
                    name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Viable Pod Number') +
  theme(axis.title.y=element_text(vjust=1), legend.position='none') +
  annotate('text', x=0.5, y=125, label='(c)', size=8)
#export at 800x800


#aborted_pods
summary(abortedModel <- lme(aborted_pods~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, aborted_pods<9000&diversity!=0), random=~1|bed))
anova(abortedModel) #diversity and warming effects
lsmeans(abortedModel, pairwise~as.factor(diversity)*as.factor(warming), adjust="tukey")

summary(glm(aborted_pods~diversity, data=subset(fitnessData, aborted_pods<9000&warming==0))) #no effect
summary(glm(aborted_pods~diversity, data=subset(fitnessData, aborted_pods<9000&warming==1))) #no effect

abortedFig <- ggplot(data=barGraphStats(data=subset(fitnessData, aborted_pods<9000&diversity!=0), variable="aborted_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                     breaks=c(1,2,3),
                     labels=c('1 strain', '2 strains', '3 strains'),
                     name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Aborted Pod Number') +
  theme(axis.title.y=element_text(vjust=1), legend.position='none') +
  annotate('text', x=0.7, y=16, label='a', size=8) +
  annotate('text', x=1, y=13, label='ab', size=8) +
  annotate('text', x=1.3, y=11, label='b', size=8) +
  annotate('text', x=1.7, y=8, label='b', size=8) +
  annotate('text', x=2, y=8, label='b', size=8) +
  annotate('text', x=2.3, y=7, label='b', size=8) +
  annotate('text', x=0.5, y=16, label='(b)', size=8)
#export at 800x800


#total beans
summary(glm(total_beans~diversity*warming, data=subset(fitnessData, total_beans<9000&diversity!=0))) #no effect
summary(glm(total_beans~diversity, data=subset(fitnessData, total_beans<9000&warming==0))) #no effect
summary(glm(total_beans~diversity, data=subset(fitnessData, total_beans<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, total_beans<9000&diversity!=0), variable="total_beans", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Total Bean Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#healthy beans
summary(healthyModel <- lme(healthy_beans~diversity*warming, data=subset(fitnessData, healthy_beans<9000&diversity!=0), random=~1|bed)) 
anova(healthyModel) #marginal interaction

summary(glm(healthy_beans~diversity, data=subset(fitnessData, healthy_beans<9000&warming==0))) #no effect
summary(glm(healthy_beans~diversity, data=subset(fitnessData, healthy_beans<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, healthy_beans<9000&diversity!=0), variable="healthy_beans", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Healthy Bean Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#damaged beans
summary(damagedModel <- lme(damaged_beans~diversity*warming, data=subset(fitnessData, damaged_beans<9000&diversity!=0), random=~1|bed))
anova(damagedModel) #no effect

summary(glm(damaged_beans~diversity, data=subset(fitnessData, damaged_beans<9000&warming==0))) #no effect
summary(glm(damaged_beans~diversity, data=subset(fitnessData, damaged_beans<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, damaged_beans<9000&diversity!=0), variable="damaged_beans", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Healthy Bean Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#aborted beans
summary(abortedBeansModel <- lme(aborted_beans~diversity*warming, data=subset(fitnessData, aborted_beans<9000&diversity!=0), random=~1|bed)) 
anova(abortedBeansModel) #no effect
summary(glm(aborted_beans~diversity, data=subset(fitnessData, aborted_beans<9000&warming==0))) #no effect
summary(glm(aborted_beans~diversity, data=subset(fitnessData, aborted_beans<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, aborted_beans<9000&diversity!=0), variable="aborted_beans", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                    breaks=c(1,2,3),
                    labels=c('1 strain', '2 strains', '3 strains'),
                    name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Aborted Beans') +
  theme(axis.title.y=element_text(vjust=1), legend.position='none')
#export at 665 x 610


#bean weight
summary(weightModel <- lme(bean_weight_g~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, bean_weight_g<900&diversity!=0), random=~1|bed)) 
anova(weightModel) #no effect
summary(glm(bean_weight_g~diversity, data=subset(fitnessData, bean_weight_g<900&warming==0))) #no effect
summary(glm(bean_weight_g~diversity, data=subset(fitnessData, bean_weight_g<900&warming==1))) #diversity effect

weightFig <- ggplot(data=barGraphStats(data=subset(fitnessData, bean_weight_g<900&diversity!=0), variable="bean_weight_g", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                    breaks=c(1,2,3),
                    labels=c('1 strain', '2 strains', '3 strains'),
                    name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Bean Weight (g)') +
  theme(axis.title.y=element_text(vjust=1), legend.position='none') +
  annotate('text', x=0.5, y=53,label='(d)', size=8)
#export at 665 x 610


#flowers, pods, and weight
pushViewport(viewport(layout=grid.layout(2,2)))
print(flowerFig2, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(abortedFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(viableFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(weightFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 1800x1800



#aphid number
summary(aphidModel <- lme(Aphids~as.factor(diversity)*as.factor(warming), data=subset(insectData, Aphids<9000&diversity!=0), random=~1|bed))
anova(aphidModel) #warming and diversity interaction
lsmeans(aphidModel, pairwise~as.factor(diversity)*as.factor(warming), adjust="tukey")

summary(glm(Aphids~diversity, data=subset(insectData, Aphids<9000&warming==0))) #no effect
summary(glm(Aphids~diversity, data=subset(insectData, Aphids<9000&warming==1))) #diversity effect

ggplot(data=barGraphStats(data=subset(insectData, Aphids<9000&diversity!=0), variable="Aphids", byFactorNames=c("warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#f37735', '#00b159', '#00aedb'),
                    breaks=c(1,2,3),
                    labels=c('1 strain', '2 strains', '3 strains'),
                    name='Rhizobial\nDiversity') +
  scale_x_discrete(breaks=c(0,1),
                   labels=c('Control', 'Drought')) +
  xlab('') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1), legend.position=c(0.02,0.98), legend.justification=c(0,1)) +
  annotate('text', x=0.7, y=220,label='a', size=8) +
  annotate('text', x=1, y=150,label='a', size=8) +
  annotate('text', x=1.3, y=180,label='a', size=8) +
  annotate('text', x=1.7, y=1000,label='b', size=8) +
  annotate('text', x=2, y=580,label='a', size=8) +
  annotate('text', x=2.3, y=450,label='a', size=8)
#export at 800x800


#######
#monocultures-----------
#insect number
summary(glm(Aphids~strains*warming, data=subset(insectData, Aphids<9000&diversity==1))) #warming and diversity interaction

ggplot(data=barGraphStats(data=subset(insectData, Aphids<9000&diversity==1), variable="Aphids", byFactorNames=c("warming", "strains")), aes(x=as.factor(strains), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


###insect herbivory---------
summary(lmer(avg_perc_herbivory~strains*warming + (1|date), data=subset(herbivoryData, diversity==1))) #no effect
summary(glm(avg_perc_herbivory~strains*warming, data=subset(herbivoryData, date=='7/27/2018'&diversity==1))) #no effect
summary(glm(avg_perc_herbivory~strains, data=subset(herbivoryData, date=='7/27/2018'&warming==0&diversity==1))) #no effect 
summary(glm(avg_perc_herbivory~strains, data=subset(herbivoryData, date=='7/27/2018'&warming==1&diversity==1))) #no effect
summary(glm(avg_perc_herbivory~strains*warming, data=subset(herbivoryData, date=='8/23/2018'&diversity==1))) #no effect

ggplot(data=barGraphStats(data=subset(herbivoryData, date %in% c('7/27/2018', '8/23/2018') & diversity==1), variable="avg_perc_herbivory", byFactorNames=c("date", "warming", "strains")), aes(x=as.factor(strains), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))
#export at 665 x 610


###rabbit herbivory---------
summary(glm(percent_rabbit_herb~strains*warming, data=subset(rabbitData, diversity==1))) #no effect
summary(glm(percent_rabbit_herb~strains, data=subset(rabbitData, warming==0&diversity==1))) #no effect 
summary(glm(percent_rabbit_herb~strains, data=subset(rabbitData, warming==1&diversity==1))) #no effect

ggplot(data=barGraphStats(data=subset(rabbitData, date=='7/19/2018'&diversity==1), variable="percent_rabbit_herb", byFactorNames=c("date", "warming", "strains")), aes(x=as.factor(strains), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Rabbit Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))
#export at 665 x 610


#height absolute
summary(lmer(height_cm~strains*warming + (1|date), data=subset(growthData, diversity==1))) #warming effect
summary(glm(height_cm~strains*warming, data=subset(growthData, date=='7/27/2018'&diversity==1))) #early season, strain 4 is different
summary(glm(height_cm~strains*warming, data=subset(growthData, date=='8/23/2018'&diversity==1))) #end of season, strain 4 is different

ggplot(data=barGraphStats(data=subset(growthData, date %in% c('7/27/2018', '8/23/2018') & diversity==1), variable="height_cm", byFactorNames=c("date", "warming", "strains")), aes(x=as.factor(strains), y=mean, color=as.factor(warming))) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Height (cm)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))
#export at 1200 x 600


#total pods
summary(glm(total_pods~strains*warming, data=subset(fitnessData, total_pods<9000&diversity==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, total_pods<9000&diversity==1), variable="total_pods", byFactorNames=c("warming", "strains")), aes(x=as.factor(strains), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Total Pod Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#######
#selection vs complementarity-----------
#insect number
monoAphids <- insectData%>%
  filter(diversity==1&Aphids<9000)%>%
  group_by(warming, strains)%>%
  summarise(mono_aphid=mean(Aphids))%>%
  ungroup()%>%
  mutate(strains2=paste('strain', strains, sep='_'))%>%
  select(-strains)%>%
  spread(key=strains2, value=mono_aphid)
expectedAphids <- insectData%>%
  left_join(monoAphids)%>%
  mutate(expected_aphids=(strain_1*USDA_110 + strain_2*USDA_76 + strain_3*USDA_136 + strain_4*USDA_138)/diversity)%>%
  select(bed, plant, diversity, strains, warming, Aphids, expected_aphids)%>%
  rename(observed_aphids=Aphids)%>%
  gather(key='type', value='number', observed_aphids, expected_aphids)

ggplot(data=barGraphStats(data=subset(expectedAphids, number<9000&diversity!=0&warming==1), variable="number", byFactorNames=c("type", "diversity", "strains")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_aphids', 'observed_aphids'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610

ggplot(data=barGraphStats(data=subset(expectedAphids, number<9000&diversity!=0&warming==1), variable="number", byFactorNames=c("type", "diversity")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2, color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_aphids', 'observed_aphids'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


###insect herbivory---------
monoChew <- herbivoryData%>%
  filter(diversity==1&date=='7/27/2018')%>%
  group_by(warming, strains)%>%
  summarise(mono_chew=mean(avg_perc_herbivory))%>%
  ungroup()%>%
  mutate(strains2=paste('strain', strains, sep='_'))%>%
  select(-strains)%>%
  spread(key=strains2, value=mono_chew)
expectedChew <- herbivoryData%>%
  filter(date=='7/27/2018')%>%
  left_join(monoChew)%>%
  mutate(expected_chew=(strain_1*USDA_110 + strain_2*USDA_76 + strain_3*USDA_136 + strain_4*USDA_138)/diversity)%>%
  select(date, bed, plant, diversity, strains, warming, avg_perc_herbivory, expected_chew)%>%
  rename(observed_chew=avg_perc_herbivory)%>%
  gather(key='type', value='number', observed_chew, expected_chew)

ggplot(data=barGraphStats(data=subset(expectedChew, date %in% c('7/27/2018') & diversity!=0), variable="number", byFactorNames=c("type", "diversity", "strains")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_chew', 'observed_chew'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610

ggplot(data=barGraphStats(data=subset(expectedChew, date %in% c('7/27/2018') & diversity!=0), variable="number", byFactorNames=c("type", "diversity")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2, color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_chew', 'observed_chew'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610