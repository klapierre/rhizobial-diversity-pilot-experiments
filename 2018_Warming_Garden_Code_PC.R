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
setwd("/Users/kathrynbloodworth/Dropbox (Smithsonian)/SERC Ecosystem Conservation/projects/Common_Garden/2018/Data")


#packages
#install.packages("tidyverse")
library(tidyverse)
#install.packages("readxl")
library(readxl)
#install.packages("summarytools")
library(summarytools)
#install.packages("lme4")
library(lme4)
#install.packages("lmerTest")
library(lmerTest)


#### Data ####
#Read in csv into new data frame and change bed_num, and warming to factors, and change date to proper date format 
Growth<- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Height_Herbivory_Measurements/Garden_Height_Leaves_Herbivory_Data_2018.csv",col_types = cols(bed_num = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16")), date = col_date(format = "%m/%d/%Y"), height_cm = col_number(), num_flowers = col_number(),num_leaves = col_number(), num_pods = col_number(),num_rabbit_herb = col_number(), warming = col_factor(levels = c("0","1")))) #might need to adjust date format (Y or y, depending on whether year is read in as 4 numbers (Y) or 2 (y))
        

#Read in csv into new data frame, changing date to proper format, and bed number and warming to factors
Herbivory <- read_csv("Data/Height_Herbivory_Measurements/Warming_Garden_Percent_Herbivory_2018.csv", 
                      col_types = cols(bed_num = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")), date = col_date(format = "%m/%d/%Y"), warming = col_factor(levels = c("1", "0")))) #might need to adjust date format (Y or y)

#Read in csv into new data frame, changing bed_num and diversity to factors
Treatment <- read_csv("Data/Height_Herbivory_Measurements/Diversity_Treatment_2018.csv",
                      col_types = cols(bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16")), diversity = col_factor(levels = c("1","2", "3", "0"))))

#Read in csv into new data frame, changing bed_num and diversity to factors
Aphids<-read_csv("Data/Aphids/2018_Insect_Counts.csv", 
                 col_types = cols(Aphids = col_number(), Large_Ants = col_number(), Predator_Ladybug = col_number(),Predator_Spider = col_number(), Small_Ants = col_number(),bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13","14", "15", "16")), date = col_date(format = "%m/%d/%Y"), warming = col_factor(levels = c("1", "0"))))%>%
  left_join(Treatment)
        
#### Data wrangling ####

#Make a new data frame called Herbivory_sum
Herbivory_sum<-Herbivory%>%
  #Filter out 9999
  filter(perc_herbivory!=9999)%>%
  #group by date, plant_num, and bed_num
  group_by(date, plant_num, bed_num, warming)%>%
  #find the sum of percent herbivory and percent rust by above groupings
  summarise(dmg_sum=sum(perc_herbivory), rust_sum=sum(percent_rust))

#Make a new data frame called Individuals from Growth data frame.  Left_join data frame Herbivory_sum
Individuals <-Growth%>%
  left_join(Herbivory_sum)%>%
  mutate(dmg_sum=ifelse(num_leaves==9999,9999,dmg_sum))%>%
  mutate(rust_sum=ifelse(num_leaves==9999,9999,rust_sum))

##Unnecessary now? Keeping incase we need to go back to it
#Make a new data frame called individuals from Growth data frame.  Left_join data frame Herbivory_sum and replace all "na"s in dmg_sum and rust_sum with 0s
#Individuals<-Growth%>%
  #left_join(Herbivory_sum)%>%
  #mutate(dmg_sum=replace_na(dmg_sum, 0), rust_sum=replace_na(rust_sum, 0))
#individuals[individuals==9999]<-NA

#View summary of data frame individuals 
view(dfSummary(Individuals))   

#adding avg dmg. find a different way to change 9999 to na and combine these three steps.
Individuals2<-Individuals%>%
  mutate(dmg_avg=dmg_sum/num_leaves, rust_avg=rust_sum/num_leaves)%>%
  select(-dmg_sum, -rust_sum,-Comments)%>%
  left_join(Treatment)%>%
  select(-treatment_code)

#Averages per bed number and treatment type
Average_Height_Leaves<-Individuals2%>%
  filter(height_cm!=9999)%>%
  filter(num_leaves!=9999)%>%
  group_by(bed_num,warming,diversity,date)%>%
  summarise(Height_avg=mean(height_cm),num_leaves_avg=mean(num_leaves),num_flowers_avg=mean(num_flowers))

#make new data frame for flower number
Flower_Number<-Individuals2%>%
  filter(num_leaves!=9999)%>%
  group_by(bed_num,date)%>%
  mutate(num_plants=length(plant_num))%>%
  ungroup()%>%
  select(bed_num,date,warming,num_flowers,diversity,num_plants)%>%
  group_by(bed_num,date,num_plants)%>%
  summarise(flower_sum=sum(num_flowers))%>%
  ungroup()%>%
  mutate(flower_avg=flower_sum/num_plants)

Pod_Number_Avg<-Individuals2%>%
  select(bed_num,date,warming,num_pods,diversity,plant_num)%>%
  filter(num_pods!=9999)%>%
  group_by(bed_num,date)%>%
  mutate(num_plants=length(plant_num))%>%
  ungroup()%>%
  group_by(bed_num,date,num_plants)%>%
  summarise(pod_sum=sum(num_pods))%>%
  ungroup()%>%
  mutate(Pod_avg=pod_sum/num_plants)

Pod_Number<-Individuals2%>%
  left_join(Pod_Number_Avg)

#Make a new data frame to join Aphid data with treatment and warming data
Aphid_Avg<-Aphids%>%
  left_join(Treatment)%>%
  filter(Aphids!=9999)%>%
  group_by(bed_num,warming,diversity,date)%>%
  summarise(aphid_avg=mean(Aphids),large_ants_avg=mean(Large_Ants),small_ants_avg=mean(Small_Ants),predator_spider_avg=mean(Predator_Spider),predator_ladybug_avg=mean(Predator_Ladybug))


#Make new data frame for error bars
Growth_data_Summary<-Individuals2%>%
  #Group data by the columns "Watershed" and "exclosure"
  group_by(bed_num,warming,diversity)%>%
  #In this data frame, summarize the data.  Make a new column named "Richness_Std" and calculate the standard deviation from the column "Richness".  Also calculate the mean and length from the column "Richness" and place them into their own columns.
  summarize(num_leaves_Std=sd(num_leaves),num_leaves_Mean=mean(num_leaves),num_leaves_n=length(num_leaves))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(num_leaves_St_Error=num_leaves_Std/sqrt(num_leaves_n))

#### graphing ####

theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=20), plot.title =
               element_text(size=20, vjust=2), panel.grid.major= element_line("grey"
               ),
             panel.grid.minor=element_blank(), legend.title=element_text(face="bold"),
             legend.text=element_text(size=30), panel.background=element_rect((fill="white")))

#Make plot from data frame individuals2 in order to graph the height of each plant by date with warming/not in two different colors
ggplot(subset(Individuals2,height_cm!=9999), aes(date, height_cm, color=(warming)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Height (cm)"
  ylab("Height(cm)")

#FACET WRAP Make plot from data frame individuals2 in order to graph the height of each plant by date with warming/not in two different colors
ggplot(subset(Individuals2,height_cm!=9999), aes(date, height_cm, color=(diversity)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Height (cm)"
  ylab("Height(cm)")+
  facet_wrap(~warming)


#Make a plot from data frame Average_Height_Leaves in order to graph average height by rhizobial diversity treatment
ggplot(Average_Height_Leaves, aes(as.factor(diversity), Height_avg,fill=warming))+
  geom_bar(stat = "identity", position="dodge")+
  #Label the x-axis "Rhizobial Diversity"
  xlab("Rhizobial Diversity")+
  #Label the y-axis "Height(cm)"
  ylab("Height(cm)")+
  expand_limits(y=80)

#Make plot from data frame individuals2 in order to graph the number of leaves of each plant by date with warming/not in two different colors
ggplot(subset(Individuals2,num_leaves!=9999), aes(date, num_leaves, color=as.factor(warming)))+
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Number of Leaflets"
  ylab("Number of Leaflets")

#Make plot from data frame individuals2 in order to graph the number of leaves of each plant by date with warming/not in two different colors
ggplot(subset(Individuals2,num_leaves!=9999), aes(date, num_leaves, color=as.factor(diversity)))+
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Number of Leaflets"
  ylab("Number of Leaflets")+
  facet_wrap(~warming)

#Make a plot from data frame individuals2 in order to graph avg numberof leaves per plant by rhizobial diversity treatment
ggplot(Average_Height_Leaves, aes(as.factor(diversity), num_leaves_avg,fill=warming))+
  geom_bar(stat = "identity", position="dodge")+
#geom_errorbar((aes(ymin=num_leaves_Mean-num_leaves_St_Error,ymax=num_leaves_Mean+num_leaves_St_Error),position=position_dodge(0.9),width=0.2))+
  #Label the x-axis "Rhizobial Diversity"
xlab("Rhizobial Diversity")+
  #Label the y-axis "Number of Leaflets"
  ylab("Numberof Leaflets")+
  expand_limits(y=120)

#Make plot from data frame individuals2 in order to graph the average damage per leaf by date with warming/not in two different colors
ggplot(subset(Individuals2, dmg_avg<=75 & date!="2018-08-23"), aes(date, dmg_avg, color=as.factor(warming)))+ #subsetted outliers ****need to be checked - % above 100
  geom_jitter()+
  #Add smooth line across data
  geom_smooth()+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Average % Herbivory per Leaf")+
  expand_limits(y=50)

#Make plot from data frame individuals2 in order to graph the average damage per leaf by date with warming/not in two different colors
ggplot(subset(Individuals2, dmg_avg<=75& date!="2018-08-23"), aes(date, dmg_avg, color=as.factor(diversity)))+ #subsetted outliers ****need to be checked - % above 100
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Average % Herbivory per Leaf")+
  expand_limits(y=50)+
  facet_wrap(~warming)
  

#Make plot from data frame individuals2 in order to graph the number of flowers of each plant by date with warming/not in two different colors
ggplot(subset(Individuals2, num_flowers<=100), aes(x=warming, y=num_flowers))+#subsetted outliers ****need to remove 9999s
  geom_bar(stat = "identity", position="dodge")+
#Label the x-axis "Date"
  xlab("Warming")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Number of Flowers")+
  expand_limits(y=100)

#Number of Leaves by date
ggplot(subset(Individuals2, num_flowers<=200), aes(date, num_flowers, color=as.factor(diversity)))+ 
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Number of Flowers")+
  expand_limits(y=50)+
  facet_wrap(~warming)

#Make a plot from data frame individuals2 in order to graph avg number of flowers per plant by rhizobial diversity treatment
ggplot(Average_Height_Leaves, aes(as.factor(diversity), num_flowers_avg,fill=warming))+
  geom_bar(stat = "identity", position="dodge")+
  #Label the x-axis "Rhizobial Diversity"
  xlab("Rhizobial Diversity")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Number of Flowers")+
  expand_limits(y=100)

#Make a plot from data frame individuals2 in order to graph avg number of flowers per plant by rhizobial diversity treatment
ggplot(Aphid_Avg, aes(as.factor(diversity), aphid_avg,fill=warming))+##not aphid average
  geom_boxplot()+
  #Label the x-axis "Rhizobial Diversity"
  xlab("Rhizobial Diversity")+
  #Label the y-axis "Average Number of Aphids"
  ylab("Average Number of Aphids")

#Make a plot from data frame individuals2 in order to graph avg number of flowers per plant by rhizobial diversity treatment
ggplot(subset(Pod_Number, num_pods!=9999), aes(as.factor(diversity),num_pods,fill=warming))+
  geom_boxplot()+
  #Label the x-axis "Rhizobial Diversity"
  xlab("Rhizobial Diversity")+
  #Label the y-axis "Average Number of Aphids"
  ylab("Pods")+
  facet_wrap(~date,drop = TRUE)

#Number of Leaves by date
ggplot(subset(Individuals2, num_pods!=9999), aes(date, num_pods, color=as.factor(diversity)))+ 
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Number of Pods")+
  facet_wrap(~warming)

#Number of Leaves by date
ggplot(subset(Individuals2, num_pods<=200), aes(date, num_pods, color=as.factor(diversity)))+ 
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(se=F)+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Average % Herbivory per Leaf"
  ylab("Number of pods")+
  expand_limits(y=50)+
  facet_wrap(~warming, drop = TRUE)

########


summary(aov(Aphids~diversity*warming,data = Aphids))
lmer(Aphids~diversity*warming+(1|bed_num),data = Aphids)

summary (Mixed_Model_Aphids <- lmer(Aphids~diversity*warming+(1|bed_num),data = Aphids))
anova(Mixed_Model_Aphids)

summary(Mixed_Model_Height<- lmer(height_cm~diversity*warming+(1|plant_num),data=Individuals2))
anova(Mixed_Model_Height)

summary(Mixed_Model_Leaves<-lmer(num_leaves~diversity*warming*date+(1|plant_num),data=Individuals2))
anova(Mixed_Model_Leaves)

summary(Mixed_Model_Dmg<-lmer(dmg_avg~diversity*warming*date+(1|plant_num),data=subset(Individuals2, dmg_avg!="NA"& dmg_avg!="NaN"& dmg_avg!="Inf"& date!="2018-08-23")))
anova(Mixed_Model_Dmg)

summary (Mixed_Model_Pods_09_07 <- lmer(num_pods~diversity*warming+(1|bed_num),data = subset(Pod_Number, date=="2018-09-07"&num_pods!=9999)))
anova(Mixed_Model_Pods)


summary (Mixed_Model_Pods_09_20 <- lmer(num_pods~diversity*warming+(1|bed_num),data = subset(Pod_Number, date=="2018-09-20"&num_pods!=9999)))
anova(Mixed_Model_Pods)
