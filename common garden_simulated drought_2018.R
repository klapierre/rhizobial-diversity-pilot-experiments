################################################################################
##  common garden_simulated drought_2018.R: Examining the effects of rhizobial diversity and simulated drought on resistance to herbivory and growth rates.
##
##  Author: Kimberly Komatsu
##  Date created: July 29, 2018
################################################################################

library(lme4)
library(lmerTest)
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
summary(heightModel <- lmer(log10(height_cm)~diversity*warming + doe + (doe|bed), data=subset(growthData, doy!=200)))
anova(heightModel)
check_model(heightModel)

# #not repeated measures - test
# summary(heightModel2 <- lmer(log10(height_cm)~diversity*warming +(1|bed), data=subset(growthData, doe==84)))
# anova(heightModel2)
# 
# plot_model(heightModel2, type='pred', terms = c("diversity", "warming"))

#diversity effect
heightModelDiversity <- ggpredict(heightModel, terms = c("doe", "bed", "diversity"))%>% # this gives overall predictions for the model
  rename(doe=x, bed=group, diversity=facet)

#why don't the lines have non-parallel slopes!?
heightDiversityPlot <- ggplot(data=subset(growthData, doy!=200), aes(x=doe, y=height_cm, color=as.factor(diversity))) +
  geom_point() +
  # geom_line(data = cbind(growthData, pred = predict(heightModel)), aes(y = pred), size = 1) + 
  geom_smooth(data=heightModelDiversity, aes(x=doe, y=10^(predicted), group=as.factor(bed), color=as.factor(diversity)), size=3, method='lm', se=F) +
  scale_color_manual(values=c('#ffc425', '#f37735', '#00b159', '#00aedb'), breaks=c(0,1,2,3), labels=c('uninoc.', '1 strain', '2 strains', '3 strains')) +
  xlab('Day of Experiment') + ylab('Height (cm)') +
  scale_y_continuous(trans='log10')

ggplot(data=heightModelDiversity, aes(x=doe, y=10^predicted)) +
  geom_point()
  geom_smooth(method='lm', se=F)


#warming effect
heightModelWarming <- ggpredict(heightModel, terms = c("doe", "warming"), type='re')%>% # this gives overall predictions for the model
  rename(doe=x, warming=group)

heightWarmingPlot <- ggplot(data=subset(growthData, doy!=200), aes(x=doe, y=height_cm, color=as.factor(warming))) +
  geom_point() +
  geom_line(data=heightModelWarming, aes(x=doe, y=10^(predicted), color=as.factor(warming)), size=3) +
  scale_color_manual(values=c('#666666', '#cc0000'), breaks=c(0,1), labels=c('ambient', 'warmed')) +
  xlab('Day of Experiment') + ylab('Height (cm)') +
  scale_y_continuous(trans='log10')
#export at 1200 x 600


##growth rate (cm/day)---------
#data check
ggplot(data=growthRate, aes(x=rate_period, y=growth_rate)) +
  geom_point() +
  facet_wrap(~plant)

summary(growthModel <- lmer(growth_rate~diversity*warming + (1|bed/plant), data=subset(growthRate, rate_period!='rate_3'))) #no effects
check_model(growthModel)


###leaf number---------
#data check
ggplot(data=growthData, aes(x=doy, y=num_leaves)) +
  geom_point() +
  facet_wrap(~plant)

#model
summary(leafModel <- lmer(log(num_leaves)~diversity*warming + (1|bed/plant), data=subset(growthData, doy!=200 & num_leaves>0))) #no effect
check_model(leafModel)



###insect herbivory---------
#data check
ggplot(data=subset(herbivoryData, doy!=200), aes(x=doy, y=avg_perc_herbivory)) +
  geom_point() +
  facet_wrap(~plant)

#model
summary(insectherbModel <- lmer(sqrt(avg_perc_herbivory)~diversity*warming +doe + (doe|bed), data=subset(herbivoryData, doy!=200)))
anova(insectherbModel) #no effect (neither for pre rabbit, post rabbit, or any certain date)
check_model(insectherbModel)

# #nlme option
# library(nlme)
# library(car)
# correlation = Structure(form  = ~ time | subjvar)
# 
# model = lme(sqrt(avg_perc_herbivory) ~ diversity*warming*doy,
#             random = ~1|bed/plant,
#             correlation = corAR1(),
#             data=subset(herbivoryData, doy!=200 & diversity>0),
#             method="REML")
# Anova(model)



ggplot(data=barGraphStats(data=subset(herbivoryData, doy!=200 & diversity>0), variable="avg_perc_herbivory", byFactorNames=c("doy", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~doy, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))
#export at 665 x 610



###rabbit herbivory---------
leavesPreRabbit <- growthData%>%
  filter(date=='7/11/2018')%>%
  rename(num_leaves_pre=num_leaves)%>%
  select(bed, warming, plant, num_leaves_pre)
rabbitData <- growthData%>%
  filter(date=='7/19/2018')%>%
  left_join(leavesPreRabbit)%>%
  mutate(percent_rabbit_herb=(num_leaves_pre-num_leaves)/num_leaves_pre*100)

summary(glm(percent_rabbit_herb~diversity*warming, data=subset(rabbitData, diversity!=0))) #no effect
summary(glm(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==0))) #marginally sig diversity effect 
summary(glm(percent_rabbit_herb~diversity, data=subset(rabbitData, warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(rabbitData, date=='7/19/2018'&diversity!=0), variable="percent_rabbit_herb", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Rabbit Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))
#export at 665 x 610


###flower number---------
summary(glm(num_flowers~diversity*warming, data=subset(growthData, num_flowers<9000&date=='8/23/2018'&diversity!=0))) #no effect
summary(glm(num_flowers~diversity, data=subset(growthData, num_flowers<9000&date=='8/23/2018'&warming==0))) #no effect
summary(glm(num_flowers~diversity, data=subset(growthData, num_flowers<9000&date=='8/23/2018'&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(growthData, date=='8/23/2018'&diversity!=0), variable="num_flowers", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Flower Number') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
  theme(strip.text.x = element_text(size=20),
        strip.background = element_rect(colour="black", fill="white"),
        panel.spacing = unit(2, "lines"))


###pods---------
#total pods
summary(glm(total_pods~diversity*warming, data=subset(fitnessData, total_pods<9000))) #no effect
summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==0))) #no effect
summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, total_pods<9000), variable="total_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Total Pod Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#viable pods---------
summary(glm(viable_pods~diversity*warming, data=subset(fitnessData, viable_pods<9000&diversity!=0))) #no effect
summary(glm(viable_pods~diversity, data=subset(fitnessData, viable_pods<9000&warming==0))) #no effect
summary(glm(viable_pods~diversity, data=subset(fitnessData, viable_pods<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, viable_pods>0&diversity!=0), variable="viable_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Viable Pod Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#aborted_pods
summary(glm(aborted_pods~diversity*warming, data=subset(fitnessData, aborted_pods<9000&diversity!=0))) #diversity and warming effects
summary(glm(aborted_pods~diversity, data=subset(fitnessData, aborted_pods<9000&warming==0))) #no effect
summary(glm(aborted_pods~diversity, data=subset(fitnessData, aborted_pods<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, aborted_pods<9000&diversity!=0), variable="aborted_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Aborted Pod Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


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
summary(glm(healthy_beans~diversity*warming, data=subset(fitnessData, healthy_beans<9000&diversity!=0))) #marginal warming effect
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
summary(glm(damaged_beans~diversity*warming, data=subset(fitnessData, damaged_beans<9000&diversity!=0))) #no effect
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
summary(glm(aborted_beans~diversity*warming, data=subset(fitnessData, aborted_beans<9000&diversity!=0))) #no effect
summary(glm(aborted_beans~diversity, data=subset(fitnessData, aborted_beans<9000&warming==0))) #no effect
summary(glm(aborted_beans~diversity, data=subset(fitnessData, aborted_beans<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, aborted_beans<9000&diversity!=0), variable="aborted_beans", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Healthy Bean Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#bean weight
summary(glm(bean_weight_g~diversity*warming, data=subset(fitnessData, bean_weight_g<9000&diversity!=0))) #warming and diversity effect
summary(glm(bean_weight_g~diversity, data=subset(fitnessData, bean_weight_g<9000&warming==0))) #no effect
summary(glm(bean_weight_g~diversity, data=subset(fitnessData, bean_weight_g<9000&warming==1))) #diversity effect

ggplot(data=barGraphStats(data=subset(fitnessData, bean_weight_g<9000&diversity!=0), variable="bean_weight_g", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Bean Weight (g)') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


#insect number
summary(glm(Aphids~diversity*warming, data=subset(insectData, Aphids<9000&diversity!=0))) #warming and diversity interaction
summary(glm(Aphids~diversity, data=subset(insectData, Aphids<9000&warming==0))) #no effect
summary(glm(Aphids~diversity, data=subset(insectData, Aphids<9000&warming==1))) #diversity effect

ggplot(data=barGraphStats(data=subset(insectData, Aphids<9000&diversity!=0), variable="Aphids", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Aphid Number') +
  theme(axis.title.y=element_text(vjust=1))
#export at 665 x 610


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