################################################################################
##  common garden_simulated drought_2018.R: Examining the effects of rhizobial diversity and simulated drought on resistance to herbivory and growth rates.
##
##  Author: Kimberly Komatsu
##  Date created: July 29, 2018
################################################################################

library(nlme)
library(lsmeans)
library(pbkrtest)
library(car)
# library(ggeffects)
# library(sjPlot)
# library(grid)
library(performance)
# library(piecewiseSEM)
library(tidyverse)
library(ggpubr)

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


options(contrasts=c('contr.sum','contr.poly'))

#not in function
`%!in%` = Negate(`%in%`)

###########################################################################
###########################################################################

###read in data
trt <- read.csv('soy pilot_2018_common garden_treatments.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  mutate(strains_code=ifelse(strains=='1', 'A',
                             ifelse(strains=='2', 'B',
                                    ifelse(strains=='3', 'C',
                                           ifelse(strains=='4', 'D',
                                                  ifelse(strains=='1,2', 'AB',
                                                         ifelse(strains=='1,3', 'AC',
                                                                ifelse(strains=='1,4', 'AD',
                                                                       ifelse(strains=='2,3', 'BC',
                                                                              ifelse(strains=='2,4', 'BD',
                                                                                     ifelse(strains=='3,4', 'CD',
                                                                                            ifelse(strains=='1,2,3', 'ABC',
                                                                                                   ifelse(strains=='1,2,4', 'ABD',
                                                                                                          ifelse(strains=='1,3,4', 'ACD',
                                                                                                                 ifelse(strains=='2,3,4', 'BCD', 'control')))))))))))))))
  
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
  select(doy_cat, bed, warming, plant, treatment_code, strains, strains_code, diversity, USDA_110, USDA_76, USDA_136, USDA_138, height_cm)%>%
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
  summarise(sum_perc_herbivory=sum(perc_herbivory))%>%
  ungroup()%>%
  left_join(growthData)%>%
  mutate(avg_perc_herbivory=sum_perc_herbivory/num_leaflets)%>%
  filter(!is.na(avg_perc_herbivory), avg_perc_herbivory<100000)%>%
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))

growthDataInsectCovariate <- growthData%>%
  filter(doe==52)%>%
  select(plant, height_cm, num_leaflets)

insectData <- read.csv('soy pilot_2018_common garden_insects.csv')%>%
  mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
  left_join(trt)%>%
  filter(Aphids<9000)%>%
  filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
  mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))%>%
  left_join(growthDataInsectCovariate)

# flowerData <- read.csv('soy pilot_2018_common garden_growth.csv')%>%
#   mutate(bed=as.factor(bed), plant=as.factor(plant), warming=as.factor(warming))%>%
#   left_join(trt)%>%
#   select(-height_cm, -num_leaves)%>%
#   filter(num_flowers<9000, num_pods<9000)%>%
#   filter(plant %!in% c(635,900,1062,1073,2557))%>% #drop plants that died
#   mutate(date2=as.Date(date, format='%m/%d/%Y'), doy=as.numeric(strftime(date2, format = "%j")), doe=(doy-169))

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
# #data check
# ggplot(data=growthData, aes(x=doe, y=height_cm)) +
#   geom_point() +
#   facet_wrap(~plant)
# hist(log10(growthData$height_cm))

#model -- height at final date of measurement (pre-senescence)
summary(heightModel <- lme(log10(height_cm)~as.factor(diversity)*as.factor(warming),
                           data=subset(growthData, doe==66&diversity!=0), 
                           random=~1|bed))
check_model(heightModel)
anova.lme(heightModel, type='sequential') #significant effect of warming

#height figure
ggplot(data=barGraphStats(data=subset(growthData, doe==66&diversity!=0), variable="height_cm", byFactorNames=c("warming")), aes(x=as.factor(warming), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(breaks=c(0,1), labels=c("Control", "Drought")) +
  xlab('') + ylab('Height (cm)') +
  annotate('text', x=1, y=45, label='a', size=8) +
  annotate('text', x=2, y=51, label='b', size=8)
#export at 800 x 800



###leaf number---------
# #data check
# ggplot(data=growthData, aes(x=doy, y=num_leaves)) +
#   geom_point() +
#   facet_wrap(~plant)



###insect herbivory---------
# #data check
# ggplot(data=subset(herbivoryData, doy!=200), aes(x=doy, y=avg_perc_herbivory)) +
#   geom_point() +
#   facet_wrap(~plant)

#model
### start here and think about correlation structure
summary(insectherbModel <- lme(sqrt(avg_perc_herbivory)~as.factor(diversity)*as.factor(warming),
                               data=subset(herbivoryData, doy!=200&diversity!=0), 
                               random=~1|bed, 
                               correlation=corAR1(form=~doe|bed/plant),
                               control=lmeControl(returnObject=T)))
check_model(insectherbModel)
anova.lme(insectherbModel, type='sequential') #significant effect of diversity and marginally significant effect of warming
lsmeans(insectherbModel, pairwise~as.factor(diversity), adjust="tukey")


# #what you'd do with lmer
# summary(mod.pot<-lmer(sqrt(avg_perc_herbivory)~as.factor(diversity)*as.factor(warming) +
#                 (1|bed)+(1|doe), #this keeps random effects seperate
#               data=subset(herbivoryData, doy!=200&diversity!=0)))
# anova(mod.pot)
# means<-emmeans(mod.pot,~as.factor(diversity),mode = "satterthwaite") #gets the model estimates (can also backtransform somehow) to generate effect sizes (for example, could say that diversity level 3 has X% less damage than diversity level 1)
# plot(means, comparisons=T)


#insect damage figure - diversity
insectFig <- ggplot(data=barGraphStats(data=subset(herbivoryData, doy!=200 & diversity!=0), variable="avg_perc_herbivory", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  xlab('') + ylab('Insect Damage (%)') +
  annotate('text', x=1, y=3, label='a', size=8) +
  annotate('text', x=2, y=2.7, label='ab', size=8) +
  annotate('text', x=3, y=2.4, label='b', size=8) +
  annotate('text', x=0.5, y=3, label='(b)', size=8)
#export at 800 x 800

###aphid number---------
summary(aphidModel <- lme(Aphids~as.factor(diversity)*as.factor(warming), data=subset(insectData, diversity!=0), random=~1|bed))
check_model(aphidModel)
anova.lme(aphidModel, type='sequential') #warming and diversity interaction
lsmeans(aphidModel, pairwise~as.factor(diversity)*as.factor(warming), adjust="tukey")

summary(glm(Aphids~diversity, data=subset(insectData, warming==0))) #no effect
summary(glm(Aphids~diversity, data=subset(insectData, warming==1))) #diversity effect


aphidFig <- ggplot(data=barGraphStats(data=subset(insectData, diversity!=0), variable="Aphids", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, fill=as.factor(warming))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  scale_fill_manual(values=c('#666666', '#cc0000'),
                    breaks=c(0,1),
                    labels=c('control', 'drought'),
                    name='Drought\nTreatment') +
  xlab('') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1), legend.position=c(0.98,0.98), legend.justification=c(1,1), 
        axis.title.x=element_blank(), axis.text.x=element_blank()) +
  annotate('text', x=0.76, y=220,label='a', size=8) +
  annotate('text', x=1.24, y=1000,label='b', size=8) +
  annotate('text', x=1.76, y=150,label='a', size=8) +
  annotate('text', x=2.24, y=580,label='ab', size=8) +
  annotate('text', x=2.76, y=180,label='a', size=8) +
  annotate('text', x=3.24, y=450,label='a', size=8) +
  annotate('text', x=0.5, y=1000, label='(a)', size=8)
#export at 800x800


###rabbit herbivory---------
leavesPreRabbit <- growthData%>%
  filter(date=='7/11/2018')%>%
  rename(num_leaves_pre=num_leaves)%>%
  select(bed, warming, plant, num_leaves_pre)
rabbitData <- growthData%>%
  filter(date=='7/19/2018')%>%
  left_join(leavesPreRabbit)%>%
  mutate(leaves_rabbit_removed=(num_leaves_pre-num_leaves), 
         percent_rabbit_removed=(num_leaves_pre-num_leaves)/num_leaves_pre*100)%>% #calculates percent of leaves removed by rabbit
  mutate(leaves_rabbit_removed=ifelse(leaves_rabbit_removed<0, 0, leaves_rabbit_removed),
         percent_rabbit_herb=ifelse(percent_rabbit_removed<0, 0, percent_rabbit_removed)) #sets any plant that increased in leaf number after rabbit herbivory to 0 percent removed
  
#percent rabbit damage
summary(rabbitModel <- lme(percent_rabbit_herb~as.factor(diversity)*as.factor(warming),
            data=subset(rabbitData, diversity!=0),
            random=~1|bed))
check_model(rabbitModel)
anova.lme(rabbitModel, type='sequential')  #no effect

summary(lme(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==0), random=~1|bed)) #no effect 
summary(lme(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==1), random=~1|bed)) #no effect

rabbitFig <- ggplot(data=barGraphStats(data=subset(rabbitData, date=='7/19/2018'&diversity!=0), variable="percent_rabbit_herb", byFactorNames=c("diversity")), aes(x=as.factor(diversity), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) +
  scale_x_discrete(breaks=c(1,2,3), labels=c("1 strain", "2 strains", "3 strains")) +
  xlab('Rhizobial Diversity') + ylab('Rabbit Damage (%)') +
  annotate('text', x=0.5, y=45, label='(c)', size=8)
#export at 800 x 800

#herbivory figure
ggarrange(aphidFig, insectFig, rabbitFig,
          ncol = 1, nrow = 3)
#export at 800x1800



#viable pods---------
summary(viableModel <- lme(viable_pods~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, diversity!=0), random=~1|bed))
check_model(viableModel)
anova.lme(viableModel, type='sequential') #no effect

#aborted_pods
summary(abortedModel <- lme(aborted_pods~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, aborted_pods<9000&diversity!=0), random=~1|bed))
check_model(abortedModel)
anova.lme(abortedModel, type='sequential') #diversity and warming effects
lsmeans(abortedModel, pairwise~as.factor(diversity)*as.factor(warming), adjust="tukey")

#healthy beans
summary(healthyModel <- lme(healthy_beans~diversity*warming, data=subset(fitnessData, diversity!=0), random=~1|bed)) 
check_model(healthyModel)
anova.lme(healthyModel, type='sequential') #marginal interaction

#damaged beans
summary(damagedModel <- lme(damaged_beans~diversity*warming, data=subset(fitnessData, diversity!=0), random=~1|bed)) 
check_model(damagedModel)
anova.lme(damagedModel, type='sequential') #no effect

#aborted beans
summary(abortedBeansModel <- lme(aborted_beans~diversity*warming, data=subset(fitnessData, diversity!=0), random=~1|bed)) 
check_model(abortedBeansModel)
anova.lme(abortedBeansModel, type='sequential') #no effect

#bean weight
summary(weightModel <- lme(bean_weight_g~as.factor(diversity)*as.factor(warming), data=subset(fitnessData, bean_weight_g<900&diversity!=0), random=~1|bed)) 
check_model(weightModel)
anova.lme(weightModel, type='sequential') #no effect



###regressions between damage and fitness outcomes
#subset and merge data

#plot by doe to determine what date to select for regressions (data support averaging across all dates)
# ggplot(data=barGraphStats(data=subset(herbivoryData, diversity!=0&doy!=200), variable="avg_perc_herbivory", byFactorNames=c("warming", "diversity", "doe")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
#   facet_wrap(~doe)

herbivoryRegData <- herbivoryData%>%
  filter(doy!=200)%>%
  group_by(bed, plant, warming, diversity)%>%
  summarize(avg_insect_herb=mean(avg_perc_herbivory))%>%
  ungroup()
insectRegData <- insectData%>%select(bed, plant, warming, diversity, Aphids, height_cm, num_leaflets)
rabbitRegData <- rabbitData%>%select(bed, plant, warming, diversity, percent_rabbit_herb)
fitnessRegData <- fitnessData%>%select(bed, plant, warming, diversity, total_pods, aborted_pods, viable_pods, total_beans, healthy_beans, bean_weight_g)
regData <- herbivoryRegData%>%left_join(rabbitRegData)%>%left_join(fitnessRegData)%>%left_join(insectRegData)%>%filter(complete.cases(.))%>%filter(diversity!=0)

#are insect and rabbit herbivory correlated?
summary(lme(percent_rabbit_herb ~ sqrt(avg_insect_herb), data=regData, random=~1|bed)) #no correlation
summary(lme(percent_rabbit_herb ~ Aphids, data=regData, random=~1|bed)) #no correlation
summary(lm(sqrt(avg_insect_herb) ~ Aphids, data=regData)) #negative correlation


###residuals to account for plant size
healthyBeansResiduals <- as.data.frame(residuals(lm(healthy_beans~num_leaflets, data=regData)))%>%
  cbind(regData)
colnames(healthyBeansResiduals)[1] <-'residuals'

#healthy bean number
summary(lm(residuals ~ sqrt(avg_insect_herb), data=healthyBeansResiduals))
summary(lm(residuals ~ Aphids, data=healthyBeansResiduals)) #significant increase
summary(lm(residuals ~ percent_rabbit_herb, data=healthyBeansResiduals))


#figure
aphidBeansFig <- ggplot(data=regData, aes(x=Aphids, y=healthy_beans)) +
  geom_point(size=3, aes(shape=as.factor(diversity), color=as.factor(warming))) +
  xlab('Aphid Number') + ylab('Healthy Bean Number') +
  scale_color_manual(values=c('#666666', '#cc0000'),
                    breaks=c(0,1),
                    labels=c('control', 'drought'),
                    name='Drought\nTreatment') +
  scale_shape_manual(values=c(19,17,15),
                     breaks=c(1,2,3),
                     labels=c('1 strain', '2 strains', '3 strains'),
                     name='Diversity\nTreatment') +
  stat_smooth(method='lm', color='black', se=F) +
  scale_x_continuous(trans='log10')

insectBeansFig <- ggplot(data=regData, aes(x=avg_insect_herb, y=healthy_beans)) +
  geom_point(size=3, aes(shape=as.factor(diversity), color=as.factor(warming))) +
  xlab('Insect Damage (%)') + ylab('Healthy Bean Number') +
  scale_color_manual(values=c('#666666', '#cc0000'),
                     breaks=c(0,1),
                     labels=c('control', 'drought'),
                     name='Drought\nTreatment') +
  scale_shape_manual(values=c(19,17,15),
                     breaks=c(1,2,3),
                     labels=c('1 strain', '2 strains', '3 strains'),
                     name='Diversity\nTreatment') +
  # stat_smooth(method='lm', color='black', se=F) +
  scale_x_continuous(trans='sqrt')

rabbitBeansFig <- ggplot(data=regData, aes(x=percent_rabbit_herb, y=healthy_beans)) +
  geom_point(size=3, aes(shape=as.factor(diversity), color=as.factor(warming))) +
  xlab('Rabbit Damage (%)') + ylab('Healthy Bean Number') +
  scale_color_manual(values=c('#666666', '#cc0000'),
                     breaks=c(0,1),
                     labels=c('control', 'drought'),
                     name='Drought\nTreatment') +
  scale_shape_manual(values=c(19,17,15),
                     breaks=c(1,2,3),
                     labels=c('1 strain', '2 strains', '3 strains'),
                     name='Diversity\nTreatment')
  # stat_smooth(method='lm', color='black', se=F)

#aphid number, insect, and rabbit damage on healthy beans
ggarrange(aphidBeansFig, insectBeansFig, rabbitBeansFig,
          ncol = 3, nrow = 1,
          common.legend=T, legend='bottom')
#export at 1800x600





####### monocultures vs polycultures
#monocultures-----------
#insect number
summary(aphidsMonocultureModel <- lme(Aphids~strains_code*warming, data=subset(insectData, Aphids<9000&diversity==1), random=~1|bed)) 
check_model(aphidsMonocultureModel)
anova.lme(aphidsMonocultureModel, type='sequential') #warming 

aphidMonoculturePlot <- ggplot(data=barGraphStats(data=subset(insectData, Aphids<9000&diversity==1), variable="Aphids", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.9), width=0.2) +
  xlab('Rhizobial Strain') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1), axis.title.x=element_blank(), axis.text.x=element_blank()) +
  annotate('text', x=0.6, y=1000, label='(a)', size=8)
#export at 800x800


###insect herbivory---------
summary(insectherbModel <- lme(sqrt(avg_perc_herbivory)~as.factor(strains_code)*as.factor(warming),
                               data=subset(herbivoryData, doy!=200&diversity==1), 
                               random=~1|bed, 
                               correlation=corAR1(form=~doe|bed/plant),
                               control=lmeControl(returnObject=T)))
check_model(insectherbModel)
anova.lme(insectherbModel, type='sequential') #marginally significant effect of strain
lsmeans(insectherbModel, pairwise~as.factor(strains_code), adjust="tukey")

insectMonoculturePlot <- ggplot(data=barGraphStats(data=subset(herbivoryData, doy!=200&diversity==1), variable="avg_perc_herbivory", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(width=0.25), fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  # scale_color_manual(values=c('#0072B2', '#D55E00'),
  #                    breaks=c(0,1),
  #                    labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Invertebrate Herbivory (%)') +
  annotate('text', x=0.6, y=4.5, label='(b)', size=8)
#export at 800 x 800

# herbivoryData$strains_code <- factor(herbivoryData$strains_code, levels=c('A', 'B', 'C', 'D', 'AB', 'AC', 'AD', 'BC', 'BD', 'CD', 'ABC', 'ABD', 'BCD'))
#   
# ggplot(data=barGraphStats(data=subset(herbivoryData, doy!=200&diversity!=0), variable="avg_perc_herbivory", byFactorNames=c("strains_code", "doe")), aes(x=as.factor(strains_code), y=mean)) +
#   geom_bar(stat='identity', position=position_dodge(width=0.25), color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
#   xlab('Rhizobial Strain') + ylab('Healthy Bean Number') +
#   theme(axis.title.y=element_text(vjust=1)) +
#   facet_wrap(~doe)
# #export at 800x600


###rabbit herbivory---------
summary(rabbitMonocultureModel <- lme(percent_rabbit_herb~strains_code*warming, random=~1|bed,
                                      data=subset(rabbitData, diversity==1))) 
check_model(rabbitMonocultureModel)
anova.lme(rabbitMonocultureModel, type='sequential') #no effect 

rabbitMonoculturePlot <- ggplot(data=barGraphStats(data=subset(rabbitData, date=='7/19/2018'&diversity==1), variable="percent_rabbit_herb", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
  geom_bar(position=position_dodge(width=0.25), stat='identity', color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  # scale_color_manual(values=c('#0072B2', '#D55E00'),
  #                    breaks=c(0,1),
  #                    labels=c('Control', 'Drought')) +
  xlab('Rhizobial Strain') + ylab('Rabbit Herbivory (%)') +
  annotate('text', x=0.6, y=55, label='(c)', size=8)
#export at 800 x 800


#herbivory figure
ggarrange(aphidMonoculturePlot, insectMonoculturePlot, rabbitMonoculturePlot,
          ncol = 1, nrow = 3)
#export at 800x1800



###total pods-----------------------
summary(beansMonocultureModel <- lme(total_pods~strains_code*warming, random=~1|bed,
                                     data=subset(fitnessData, diversity==1))) 
check_model(beansMonocultureModel)
anova.lme(beansMonocultureModel, type='sequential') #strain
lsmeans(beansMonocultureModel, pairwise~as.factor(strains_code), adjust="tukey")


#healthy beans
summary(beansMonocultureModel <- lme(healthy_beans~strains_code*warming, random=~1|bed,
                                      data=subset(fitnessData, diversity==1))) 
check_model(beansMonocultureModel)
anova.lme(beansMonocultureModel, type='sequential') #strain
lsmeans(beansMonocultureModel, pairwise~as.factor(strains_code), adjust="tukey")

ggplot(data=barGraphStats(data=subset(fitnessData, diversity==1), variable="healthy_beans", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
  geom_bar(stat='identity', position=position_dodge(width=0.25), color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
    xlab('Rhizobial Strain') + ylab('Healthy Bean Number') +
  theme(axis.title.y=element_text(vjust=1)) +
  annotate('text', x=1, y=380,label='a', size=8) +
  annotate('text', x=2, y=320,label='ab', size=8) +
  annotate('text', x=3, y=335,label='ab', size=8) +
  annotate('text', x=4, y=270,label='b', size=8)
#export at 800x600

# fitnessData$strains_code <- factor(fitnessData$strains_code, levels=c('A', 'B', 'C', 'D', 'AB', 'AC', 'AD', 'BC', 'BD', 'CD', 'ABC', 'ABD', 'BCD'))
# 
# ggplot(data=barGraphStats(data=subset(fitnessData, diversity!=0), variable="healthy_beans", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
#   geom_bar(stat='identity', position=position_dodge(width=0.25), color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
#   xlab('Rhizobial Strain') + ylab('Healthy Bean Number') +
#   theme(axis.title.y=element_text(vjust=1))
# #export at 800x600


#total beans
summary(beansMonocultureModel <- lme(total_beans~strains_code*warming, random=~1|bed,
                                     data=subset(fitnessData, diversity==1))) 
check_model(beansMonocultureModel)
anova.lme(beansMonocultureModel) #strain
lsmeans(beansMonocultureModel, pairwise~as.factor(strains_code), adjust="tukey")

#bean weight
summary(beansMonocultureModel <- lme(bean_weight_g~strains_code*warming, random=~1|bed,
                                     data=subset(fitnessData, diversity==1))) 
check_model(beansMonocultureModel)
anova.lme(beansMonocultureModel) #strain
lsmeans(beansMonocultureModel, pairwise~as.factor(strains_code), adjust="tukey")

# ggplot(data=barGraphStats(data=subset(fitnessData, diversity==1), variable="bean_weight_g", byFactorNames=c("strains_code")), aes(x=as.factor(strains_code), y=mean)) +
#   geom_bar(stat='identity', position=position_dodge(width=0.25), color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
#   xlab('Rhizobial Strain') + ylab('Bean Weight (g)') +
#   theme(axis.title.y=element_text(vjust=1))
# #export at 800x600



#######
#selection vs complementarity-----------
# #aphid number
# monoAphids <- insectData%>%
#   filter(diversity==1&Aphids<9000)%>%
#   group_by(warming, strains)%>%
#   summarise(mono_aphid=mean(Aphids))%>%
#   ungroup()%>%
#   mutate(strains2=paste('strain', strains, sep='_'))%>%
#   select(-strains)%>%
#   spread(key=strains2, value=mono_aphid)
# expectedAphids <- insectData%>%
#   left_join(monoAphids)%>%
#   mutate(expected_aphids=(strain_1*USDA_110 + strain_2*USDA_76 + strain_3*USDA_136 + strain_4*USDA_138)/diversity)%>%
#   select(bed, plant, diversity, strains, warming, Aphids, expected_aphids)%>%
#   rename(observed_aphids=Aphids)%>%
#   gather(key='type', value='number', observed_aphids, expected_aphids)
# 
# ggplot(data=barGraphStats(data=subset(expectedAphids, number<9000&diversity!=0&warming==1), variable="number", byFactorNames=c("type", "diversity", "strains")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
#   geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
#   scale_shape_manual(values=c(17, 16),
#                      breaks=c('expected_aphids', 'observed_aphids'),
#                      labels=c('expected', 'observed')) +
#   xlab('Rhizobial Diversity') + ylab('Aphid Number') +
#   theme(axis.title.y=element_text(vjust=1))
# #export at 665 x 610
# 
# ggplot(data=barGraphStats(data=subset(expectedAphids, number<9000&diversity!=0&warming==1), variable="number", byFactorNames=c("type", "diversity")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
#   geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2, color='#D55E00') +
#   scale_shape_manual(values=c(17, 16),
#                      breaks=c('expected_aphids', 'observed_aphids'),
#                      labels=c('expected', 'observed')) +
#   xlab('Rhizobial Diversity') + ylab('Aphid Number') +
#   theme(axis.title.y=element_text(vjust=1))
# #export at 665 x 610


###insect herbivory---------
monoChew <- herbivoryData%>%
  filter(diversity==1&doy!=200)%>%
  group_by(warming, strains)%>%
  summarise(mono_chew=mean(avg_perc_herbivory))%>%
  ungroup()%>%
  mutate(strains2=paste('strain', strains, sep='_'))%>%
  select(-strains)%>%
  spread(key=strains2, value=mono_chew)
expectedChew <- herbivoryData%>%
  filter(doy!=200)%>%
  left_join(monoChew)%>%
  mutate(expected_chew=(strain_1*USDA_110 + strain_2*USDA_76 + strain_3*USDA_136 + strain_4*USDA_138)/diversity)%>%
  select(doy, bed, plant, diversity, strains, warming, avg_perc_herbivory, expected_chew)%>%
  rename(observed_chew=avg_perc_herbivory)%>%
  gather(key='type', value='number', observed_chew, expected_chew)

ggplot(data=barGraphStats(data=subset(expectedChew, doy!=200 & diversity!=0), variable="number", byFactorNames=c("type", "diversity", "strains")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_chew', 'observed_chew'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610

ggplot(data=barGraphStats(data=subset(expectedChew, doy!=200 & diversity!=0), variable="number", byFactorNames=c("type", "diversity")), aes(x=as.factor(diversity), y=mean, shape=as.factor(type))) +
  geom_point(size=5, position=position_dodge(width=0.25), color='#D55E00') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2, color='#D55E00') +
  scale_shape_manual(values=c(17, 16),
                     breaks=c('expected_chew', 'observed_chew'),
                     labels=c('expected', 'observed')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610
