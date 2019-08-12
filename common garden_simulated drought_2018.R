################################################################################
##  common garden_simulated drought_2018.R: Examining the effects of rhizobial diversity and simulated drought on resistance to herbivory and growth rates.
##
##  Author: Kimberly Komatsu
##  Date created: July 29, 2018
################################################################################

library(lme4)
library(lmerTest)
library(grid)
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


###########################################################################
###########################################################################

###read in data
trt <- read.csv('soy pilot_2018_common garden_treatments.csv')
growthData <- read.csv('soy pilot_2018_common garden_growth.csv')%>%
  left_join(trt)%>%
  filter(height_cm<9000, num_leaves<9000)
herbivoryData <- read.csv('soy pilot_2018_common garden_herbivory.csv')%>%
  filter(perc_herbivory<9000)%>%
  #calculate avg herbivory per leaf
  group_by(bed, plant, date)%>%
  summarise(sum_perc_herbivory=sum(perc_herbivory))%>%
  ungroup()%>%
  left_join(growthData)%>%
  mutate(avg_perc_herbivory=sum_perc_herbivory/num_leaves)%>%
  filter(!is.na(avg_perc_herbivory), avg_perc_herbivory<100000)
insectData <- read.csv('soy pilot_2018_common garden_insects.csv')%>%
  left_join(trt)
fitnessData <- read.csv('soy pilot_2018_common garden_pods.csv')%>%
  select(-comments)%>%
  left_join(read.csv('soy pilot_2018_common garden_bean weight.csv'))%>%
  select(-notes)%>%
  left_join(trt)%>%
  mutate(viable_pods=total_pods-aborted_pods)

#relative to rhizobial controls
growthRelative <- growthData%>%
  filter(diversity==0)%>%
  rename(height_ctl=height_cm, num_flowers_ctl=num_flowers, num_pods_ctl=num_pods)%>%
  select(bed, date, warming, height_ctl, num_flowers_ctl, num_pods_ctl)%>%
  right_join(growthData)%>%
  filter(diversity!=0)%>%
  mutate(height_rel=(height_cm-height_ctl)/height_ctl, flw_rel=(num_flowers-num_flowers_ctl)/num_flowers_ctl, pods_rel=(num_pods-num_pods_ctl)/num_pods_ctl)


###height---------
#height absolute
summary(lmer(height_cm~diversity*warming + (1|date), data=growthData)) #warming effect
summary(glm(height_cm~diversity*warming, data=subset(growthData, date=='7/27/2018'&diversity!=0))) #interactive effect
      summary(glm(height_cm~diversity, data=subset(growthData, date=='7/27/2018'&warming==0))) #diversity increases height
      summary(glm(height_cm~diversity, data=subset(growthData, date=='7/27/2018'&warming==1))) #no significant effect
summary(glm(height_cm~diversity*warming, data=subset(growthData, date=='8/23/2018'&diversity!=0))) #end of season, no differences anymore

# #height relative to rhizobial control
# summary(lmer(height_rel~diversity*warming + (1|date), data=growthRelative)) #warming effect
# summary(glm(height_rel~diversity*warming, data=subset(growthRelative, date=='7/27/2018'))) #interative effects (div increases height in unwarmed, but no effect on height in warmed)
#       summary(glm(height_rel~diversity, data=subset(growthRelative, date=='7/27/2018'&warming==0))) #diversity increases height
#       summary(glm(height_rel~diversity, data=subset(growthRelative, date=='7/27/2018'&warming==1))) #no significant effect
# summary(glm(height_rel~diversity*warming, data=subset(growthRelative, date=='8/23/2018'))) #end of season, no differences anymore

# ggplot(data=barGraphStats(data=subset(growthRelative, !is.na(height_rel)), variable="height_rel", byFactorNames=c("date", "warming", "diversity")), aes(x=diversity, y=mean, color=as.factor(warming))) +
#   geom_point() +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
#   scale_color_manual(values=c('#0072B2', '#D55E00')) +
#   facet_wrap(~date, scales='free')

# ggplot(data=barGraphStats(data=subset(growthRelative, date %in% c('7/27/2018', '8/23/2018') & !is.na(height_rel)), variable="height_rel", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) + 
#   geom_point(size=5, position=position_dodge(width=0.25)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
#   scale_color_manual(values=c('#0072B2', '#D55E00')) +
#   xlab('Rhizobial Diversity') + ylab('Height Relative to Control') +
#   facet_wrap(~date, scales='free')

ggplot(data=barGraphStats(data=subset(growthData, date %in% c('7/27/2018', '8/23/2018') & diversity!=0), variable="height_cm", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
    geom_point(size=5) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
    scale_color_manual(values=c('#0072B2', '#D55E00'),
                       breaks=c(0,1),
                       labels=c('Control', 'Drought')) +
    xlab('Rhizobial Diversity') + ylab('Height (cm)') +
    theme(axis.title.y=element_text(vjust=1)) +
    facet_wrap(~date, scales='free') +
    theme(strip.text.x = element_text(size=20),
          strip.background = element_rect(colour="black", fill="white"),
          panel.spacing = unit(2, "lines"))
#export at 1200 x 600


# ggplot(data=subset(growthData, date %in% c('7/27/2018', '8/23/2018') & diversity!=0), aes(x=as.factor(diversity), y=height_cm, color=as.factor(warming))) +
#   geom_boxplot() +
#   geom_point() +
#   scale_color_manual(values=c('#0072B2', '#D55E00'),
#                      breaks=c(0,1),
#                      labels=c('Control', 'Drought')) +
#   xlab('Rhizobial Diversity') + ylab('Height (cm)') +
#   theme(axis.title.y=element_text(vjust=1)) +
#   facet_wrap(~date, scales='free') +
#   theme(strip.text.x = element_text(size=20),
#         strip.background = element_rect(colour="black", fill="white"),
#         panel.spacing = unit(2, "lines"))
# #export at 1200 x 600


# ggplot(data=barGraphStats(data=subset(growthRelative, date %in% c('7/27/2018', '8/23/2018') & !is.na(height_rel)), variable="height_rel", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(warming), y=mean, fill=as.factor(diversity))) + 
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
#   facet_wrap(~date, scales='free')


###leaf number---------
summary(lmer(num_leaves~diversity*warming + (1|date), data=growthData)) #no effect
summary(glm(num_leaves~diversity*warming, data=subset(growthData, date=='7/27/2018'&diversity!=0))) #no effect
summary(glm(num_leaves~diversity, data=subset(growthData, date=='7/27/2018'&warming==0))) #no effect
summary(glm(num_leaves~diversity, data=subset(growthData, date=='7/27/2018'&warming==1))) #no effect
summary(glm(num_leaves~diversity*warming, data=subset(growthData, date=='8/23/2018'&diversity!=0))) #no effect


###insect herbivory---------
summary(lmer(avg_perc_herbivory~diversity*warming + (1|date), data=herbivoryData)) #marginally significant interaction
summary(glm(avg_perc_herbivory~diversity*warming, data=subset(herbivoryData, date=='7/27/2018'&diversity!=0))) #no effect
summary(glm(avg_perc_herbivory~diversity, data=subset(herbivoryData, date=='7/27/2018'&warming==0))) #no effect 
summary(glm(avg_perc_herbivory~diversity, data=subset(herbivoryData, date=='7/27/2018'&warming==1))) #no effect
summary(glm(avg_perc_herbivory~diversity*warming, data=subset(herbivoryData, date=='8/23/2018'&diversity!=0))) #no effect

summary(lmer(avg_perc_herbivory~diversity*warming + (1|date), data=herbivoryData)) #marginally significant

ggplot(data=barGraphStats(data=subset(herbivoryData, date %in% c('7/27/2018', '8/23/2018') & diversity!=0), variable="avg_perc_herbivory", byFactorNames=c("date", "warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
  geom_point(size=5, position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.25), width=0.2) +
  scale_color_manual(values=c('#0072B2', '#D55E00'),
                     breaks=c(0,1),
                     labels=c('Control', 'Drought')) +
  xlab('Rhizobial Diversity') + ylab('Invertebrate Herbivory (%)') +
  theme(axis.title.y=element_text(vjust=1)) +
  facet_wrap(~date, scales='free') +
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
  mutate(percent_rabbit_herb=(num_leaves_pre-num_leaves)/num_leaves_pre)

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
summary(glm(total_pods~diversity*warming, data=subset(fitnessData, total_pods<9000&diversity!=0))) #no effect
summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==0))) #no effect
summary(glm(total_pods~diversity, data=subset(fitnessData, total_pods<9000&warming==1))) #no effect

ggplot(data=barGraphStats(data=subset(fitnessData, total_pods<9000&diversity!=0), variable="total_pods", byFactorNames=c("warming", "diversity")), aes(x=as.factor(diversity), y=mean, color=as.factor(warming))) +
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