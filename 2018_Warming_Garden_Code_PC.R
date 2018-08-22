#################################################################################
# Name: Soy Bean Warming 2018 First Look
# Coder: Jamie Pullen, Kathryn Bloodworth
# Date: 22 August 2018. 
# Purpose: This code is for exploring the data collected at the warming garden during the summer of 2018. Warming x bacterial diversity experiment. Kim L, John P, Kathryn B, Jamie P, Nicole E. 
####################################################################################
#8/22/2018  Initial code setup for ingesting the plant level data, merging the data sets(growth, herbivory, and treatment), graphing, and quick stats.  



#####################################################################################
#Working Directory *setting to the main 2018 folder so that we can pull in temperature data and plant level data.

#Jamie - WD
setwd("C:/Users/pullenj/Dropbox (Smithsonian)/Common_Garden/2018")

#Kathryn -WD
setwd("/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/")

#packages
#install.packages("tidyverse")
library(tidyverse)
#install.packages("readxl")
library(readxl)
#install.packages("summarytools")
library(summarytools)
#install.packages("lme4")
library(lme4)

#### Data ####
#Read in csv into new data frame and change bed_num, and warming to factors, and change date to proper date format 
Growth<- read_csv("Data/Height_Herbivory_Measurements/Garden_Height_Leaves_Herbivory_Data_2018.csv", 
                  #might need to adjust date format (Y or y, depending on whether year is read in as 4 numbers (Y) or 2 (y))
                  col_types = cols(bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14", "15", "16")), date = col_date(format = "%m/%d/%Y"), warming = col_factor(levels = c("0","1"))))

#Read in csv into new data frame, changing date to proper format, and bed number and warming to factors
Herbivory <- read_csv("Data/Height_Herbivory_Measurements/Warming_Garden_Percent_Herbivory_2018.csv", 
                      col_types = cols(bed_num = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")), date = col_date(format = "%m/%d/%Y"), warming = col_factor(levels = c("1", "0")))) #might need to adjust date format (Y or y)

#Read in csv into new data frame, changing bed_num and diversity to factors
Treatment <- read_csv("Data/Height_Herbivory_Measurements/Diversity_Treatment_2018.csv",
                      col_types = cols(bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16")), diversity = col_factor(levels = c("1","2", "3", "0"))))

#### Data wrangling ####

#In the herbivory data frame, change any 9999 (meaning dead plant) to 0s
Herbivory[Herbivory==9999]<-0


#Make a new data frame called Herbivory_sum
Herbivory_sum<-Herbivory%>%
  #group by date, plant_num, and bed_num
  group_by(date, plant_num, bed_num)%>%
  #find the sum of percent herbivory and percent rust by above groupings
  summarise(dmg_sum=sum(perc_herbivory), rust_sum=sum(percent_rust))

#Merge Growth data frame and Herbivory_sum data frame by plant_num, date, and bed_num
#merge<-merge(Growth, Herbivory_sum, by=c("plant_num", "date", "bed_num"), all.x=T)

#Make a new data frame called individuals
#individuals<-merge%>%
  #replace any 'na' values with zero
  #mutate(dmg_sum=replace_na(dmg_sum, 0), rust_sum=replace_na(rust_sum, 0))
  
##Joined two above steps into one##
#Make a new data frame called individuals from Growth data frame.  Left_join data frame Herbivory_sum and replace all "na"s in dmg_sum and rust_sum with 0s
individuals<-Growth%>%
  left_join(Herbivory_sum)%>%
  mutate(dmg_sum=replace_na(dmg_sum, 0), rust_sum=replace_na(rust_sum, 0))


#individuals[individuals==9999]<-NA

#View summary of data frame individuals 
view(dfSummary(individuals))   

#adding avg dmg. find a different way to change 9999 to na and combine these three steps.
individuals2<-individuals%>%
  mutate(dmg_avg=dmg_sum/num_leaves, rust_avg=rust_sum/num_leaves)%>%
  select(-dmg_sum, -rust_sum)%>%
  left_join(Treatment)



#### graphing ####

#Make plot from data frame individuals2 in order to graph the height of each plant by date with warming/not in two different colors
ggplot(individuals2, aes(date, height_cm, color=(warming)))+
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()

#Make plot from data frame individuals2 in order to graph the number of leaves of each plant by date with warming/not in two different colors
ggplot(individuals2, aes(date, num_leaves, color=as.factor(warming)))+
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()

#Make plot from data frame individuals2 in order to graph the average damage per leaf by date with warming/not in two different colors
ggplot(subset(individuals2, dmg_avg<=75), aes(date, dmg_avg, color=as.factor(warming)))+ #subsetted outliers ****need to be checked
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()

#Make a plot from data frame individuals2 in order to graph number of leaves per plant by rhizobial diversity treatment
ggplot(individuals2, aes(as.factor(diversity), num_leaves))+
  geom_boxplot()+
  #separate out warmed and not warmed treatments
  facet_wrap(~warming)



########

lmer(height_cm~diversity*warming+(1|plant_num),data=individuals2)

lmer(num_leaves~diversity*warming+(1|plant_num),data=individuals2)
lmer(dmg_avg~diversity*warming*date+(1|plant_num),data=subset(individuals2, dmg_avg!="NA"& dmg_avg!="NaN"& dmg_avg!="Inf"))

