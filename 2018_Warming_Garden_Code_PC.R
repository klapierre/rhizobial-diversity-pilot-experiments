#################################################################################
# Name: Soy Bean Warming 2018 First Look
# Coder: Jamie Pullen, Kathryn Bloodworth
# Date: 22 August 2018. 
# Purpose: This code is for exploring the data collected at the warming garden during the summer of 2018. Warming x bacterial diversity experiment. Kim L, John P, Kathryn B, Jamie P, Nicole E. 
####################################################################################
#8/22/2018  Initial code setup for ingesting the plant level data, merging the data sets(growth, herbivory, and treatment), graphing, and quick stats.  



#####################################################################################
#Working Directory *setting to the main 2018 folder so that we can pull in temperature data and plant level data.
setwd("C:/Users/pullenj/Dropbox (Smithsonian)/Common_Garden/2018")

#packages
library(tidyverse)
library(readxl)
library(summarytools)
library(lme4)

#### Data ####
Growth<- read_csv("Data/Height_Herbivory_Measurements/Garden_Height_Leaves_Herbivory_Data_2018.csv", 
                  col_types = cols(date = col_date(format = "%m/%d/%y")))#, #might need to adjust date format (Y or y)
#na = "9999")

Herbivory <- read_csv("Data/Height_Herbivory_Measurements/Warming_Garden_Percent_Herbivory_2018.csv", 
                      col_types = cols(date = col_date(format = "%m/%d/%Y")))#,  #might need to adjust date format (Y or y)
#na = "9999")

Treatment <- read_csv("Data/Height_Herbivory_Measurements/Diversity_Treatment_2018.csv")

#### Data wrangling ####
Herbivory[Herbivory==9999]<-0


Herbivory_sum<-Herbivory%>%
  group_by(date, plant_num, bed_num)%>%
  summarise(dmg_sum=sum(perc_herbivory), rust_sum=sum(percent_rust))

merge<-merge(Growth, Herbivory_sum, by=c("plant_num", "date", "bed_num"), all.x=T)

individuals<-merge%>%
  mutate(dmg_sum=replace_na(dmg_sum, 0), rust_sum=replace_na(rust_sum, 0))

individuals[individuals==9999]<-NA

view(dfSummary(individuals))   

#adding avg dmg. find a different way to change 9999 to na and combine these three steps.
individuals2<-individuals%>%
  mutate(dmg_avg=dmg_sum/num_leaves, rust_avg=rust_sum/num_leaves)%>%
  select(-dmg_sum, -rust_sum)%>%
  left_join(Treatment)


#### graphing ####

ggplot(individuals2, aes(date, height_cm, color=as.factor(warming)))+
  geom_jitter()+
  geom_smooth()


ggplot(individuals2, aes(date, num_leaves, color=as.factor(warming)))+
  geom_jitter()+
  geom_smooth()

ggplot(subset(individuals2, dmg_avg<=75), aes(date, dmg_avg, color=as.factor(warming)))+ #subsetted outliers ****need to be checked
  geom_jitter()+
  geom_smooth()

ggplot(individuals2, aes(as.factor(diversity), num_leaves))+
  geom_boxplot()+
  facet_wrap(~warming)


########

lmer(height_cm~diversity*warming+(1|plant_num),data=individuals2)

lmer(num_leaves~diversity*warming+(1|plant_num),data=individuals2)
lmer(dmg_avg~diversity*warming*date+(1|plant_num),data=subset(individuals2, dmg_avg!="NA"& dmg_avg!="NaN"& dmg_avg!="Inf"))

