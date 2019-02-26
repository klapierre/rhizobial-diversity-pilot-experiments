#### Soy Div Working Group Analysis ####
#Coders: Kathryn Bloodworth, Kristina Borst, Ben Taylor
#Purspose: Learning Stuff about Soybeans

#### Working Directories ####

#Ben WD
setwd("C:/Users/Benton/Desktop/Work Files/Courses Taught/R Study Group_SERC/Soy Div_R_Coding_Group_SERC")

#Kathryn WD
#Mac
setwd("~/Dropbox (Smithsonian)/SERC Ecosystem Conservation/R_Coding_Working_Group/Soy_Div_2018_R_Coding_Group")
#PC
setwd("/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/")

#Kristina WD

#### Install Packages ####
#install the package ggplot2
#install.packages("ggplot2")
library(ggplot2)#load the package ggplot2
#install.packages("summarytools")
library(summarytools)
#install the package tidyverse
#install.packages("tidyverse")
library(tidyverse)#load the package tidyverse

#### Read in the dataframes ####
#read in "Garden_Height_Leaves_Herbivory_Data_2018.csv" and call it plnts - change the column "bed number" to factors, change the column "date" to date structure, change the column "height_cm", "num_flowers","num_leaves","num_pods" and "num_rabbit_herb" to numbers instead of integers and change the column "warming" to factors
plnts <- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Height_Herbivory_Measurements/Garden_Height_Leaves_Herbivory_Data_2018.csv",col_types = cols(bed_num = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16")), date = col_date(format = "%m/%d/%Y"),height_cm = col_number(), num_flowers = col_number(), num_leaves = col_number(), num_pods = col_number(),num_rabbit_herb = col_number(), warming = col_factor(levels = c("0","1"))))

#read in "Diversity_Treatment_2018.csv" to the data frame div - change "bed_num" to factors
div <- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Height_Herbivory_Measurements/Diversity_Treatment_2018.csv",col_types = cols(bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16"))))

#read in "Garden_Height_Leaves_Herbivory_Data_2018.csv" to the data frame herb - change "bed_bum" to factors, "date" to date format, and "warming" to factors
herb <- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Height_Herbivory_Measurements/Garden_Height_Leaves_Herbivory_Data_2018.csv", col_types = cols(bed_num = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8","9", "10", "11", "12", "13", "14","15", "16")), date = col_date(format = "%m/%d/%Y"),warming = col_factor(levels = c("1","0"))))

#read in "2018_Insect_Counts.csv" to the data frame ins - change "bed_bum" to factors, "date" to date format, and "warming" to factors
ins <- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Aphids/2018_Insect_Counts.csv",col_types = cols(bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")), date = col_date(format = "%m/%d/%Y"), warming = col_factor(levels = c("1","0"))))

#read in "Soybean_Harvest_Pods_Beans_Data_2018.csv" to the data frame soybean_harvest
soybean_harvest <- read_csv("C:/Users/bloodworthk/Dropbox (Smithsonian)/SERC Ecosystem Conservation/Projects/Common_Garden/2018/Data/Soybean_Harvest_Data/Soybean_Harvest_Pods_Beans_Data_2018.csv",col_types = cols(aborted_pods = col_number(),abrt_beans = col_number(), bed_num = col_factor(levels = c("1","2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13","14", "15", "16")), brown_pods = col_number(),dmg_beans = col_number(), green_pods = col_number(),harvest_date = col_date(format = "%m/%d/%Y"),hlthy_beans = col_number(), total_beans = col_number(),total_pods = col_number(), warming = col_factor(levels = c("1","0"))))



####BNT version - Base R
### Starting to put the dataframes together for analyses ###
#Merging the diversity treatments into the plant height dataframe
div.mrg<-div[,c(2,5)]#makes a smaller diversity dataframe with just plant number and diversity
dat<-merge(plnts,div.mrg,by="plant_num", all.x=T,all.y=T)#makes a new dataframe with plnts and div combined


#Merging the insect count data into the plant height dataframe
ins.mrg<-ins[,c(3:8)]#gets just the columns in the ins data we need to merge
dat<-merge(dat,ins.mrg,by="plant_num",all.x=T,all.y=T)#merges the insect counts into the big "dat" dataframe
setdiff(plnts$plant_num,ins.mrg$plant_num)

#Merging the final assessment of herbivory into the dataframe
last.herb<-herb[herb$date=="8/24/2018",]
herb.mrg<-data.frame("plant_num"=sort(unique(last.herb$plant_num)),
                     "perc_herbivory"=with(last.herb, tapply(perc_herbivory,plant_num,mean, na.rm=T)))
dat<-merge(dat,herb.mrg,by="plant_num",all.x=T,all.y=T)


dat.fin<-dat[dat$date=="2018-09-20",]
dat.fin<-dat.fin[!is.na(dat.fin$plant_num),]

dat.fin$nasties<-(dat.fin$Aphids+dat.fin$Large_Ants+dat.fin$Small_Ants)

### KB Version - How to do above steps with tidyverse - main difference is that you don't need to make multiple data frames - you can merge two directly together and then remove columns - we could even combine step two and add in data frame "ins" in this step 
#make a new data frame called "plnts_div" from data frame plnts
plnts_div<-plnts%>%
  #join together data frame "div" using whatever columns the two data frames have in common
  left_join(div)%>%
  #remove the columns listed below
  select(-strains,-Strain_1,-Strain_2,-Strain_3,-Strain_4)%>%
  #make a new column called "div_trtmnt" and paste the columns "diversity" and "treatment", separating them with a period
  mutate(div_trtmnt=paste(diversity,treatment_code,sep = "."))%>%
  #make a new column called "trtmnt_mono" and if the column diversity is 0, place a 0 in the new column, if the diversity is 2, place a 2 in the new column, if the diversity is 3, place a 3 in the new column - for all other empty cells, copy what is in the column "div_trtmnt" -- this allows us to tease apart the mono culture data 
  mutate(trtmnt_mono=ifelse(diversity==0,0, ifelse(diversity==2,2,ifelse(diversity==3,3,div_trtmnt))))%>%
  select(-div_trtmnt)

#Make a new data frame called "avg_height" using the data frame "plnts_div"
avg_height<-plnts_div%>%
  #filter out any 9999 in the column "height_cm"
  filter(height_cm!=9999)%>%
  #group by "bed_num", "warming", "div_trtmnt", and "date" 
  group_by(warming,trtmnt_mono,date)%>%
  #take the mean of "height_cm" according to above groupings and put in new column called "height_avg", take length, standard deviation and standard error
  summarise(height_length=length(height_cm),height_avg=mean(height_cm),height_std=sd(height_cm),height_sderror=height_std/sqrt(height_length))%>%
  #ungroup the data that was previously grouped
  ungroup()
#NOT TIDYVERSE - this is from the package Summary Tools
view(dfSummary(avg_height))

#Make plot from data frame avg_height in order to graph the height of each plant by diversity treatment with warming/not in two different colors#
ggplot(avg_height, aes(x=trtmnt_mono, y=height_avg, fill=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Average Height(cm)")#+
#add error bars
#geom_errorbar(aes(ymin=Richness_Mean-Richness_St_Error,ymax=Richness_Mean+Richness_St_Error),position=position_dodge(0.9),width=0.2)+

#Make plot from data frame avg_height in order to graph the height of each plant by date with warming/not in two different colors#
ggplot(avg_height, aes(x=date, y=height_avg, color=(warming)))+#filter out 9999
  geom_jitter()+
  #make a smooth line across the data with the linear model model
  geom_smooth(method = "lm")+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Height (cm)"
  ylab("Average Height(cm)")


#Make a new data frame called "avg_flowers" using the data frame "plnts_div"
avg_flowers<-plnts_div%>%
  #filter out any 9999 in the column "num_flowers"
  filter(num_flowers!=9999)%>%
  #group by "bed_num", "warming", "div_trtmnt", and "date" 
  group_by(bed_num,warming,trtmnt_mono,date)%>%
  #take the mean of "num_flowers" according to above groupings and put in new column called "num_flowers_avg"
  summarise(num_flowers_avg=mean(num_flowers))
#NOT TIDYVERSE - this is from the package Summary Tools
view(dfSummary(avg_flowers)) 

###Make a quick ggplot graph to visualize our data
#Make plot from data frame avg_flowers in order to graph the height of each plant by date with warming/not in two different colors#
ggplot(avg_flowers, aes(x=date, y=num_flowers_avg, color=(warming)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(method="loess")+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Height (cm)"
  ylab("Average Height(cm)")

#Make plot from data frame avg_flowers in order to graph the height of each plant by diversity treatment with warming/not in two different colors#
ggplot(avg_flowers, aes(x=trtmnt_mono, y=num_flowers_avg, color=(warming)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(method = "lm")+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Average Flowers")

#Make a new data frame called "avg_leaves" using the data frame "plnts_div"
avg_leaves<-plnts_div%>%
  #filter out any 9999 in the column "num_leaves"
  filter(num_leaves!=9999)%>%
  #group by "bed_num", "warming", "div_trtmnt", and "date" 
  group_by(bed_num,warming,trtmnt_mono,date)%>%
  #take the mean of "num_leaves" according to above groupings and put in new column called "num_leaves_avg"
  summarise(num_leaves_avg=mean(num_leaves))
#NOT TIDYVERSE - this is from the package Summary Tools
view(dfSummary(avg_leaves)) 

###Make a quick ggplot graph to visualize our data
#Make plot from data frame avg_leaves in order to graph the height of each plant by date with warming/not in two different colors#
ggplot(avg_leaves, aes(x=date, y=num_leaves_avg, color=(warming)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(method="lm")+
  #Label the x-axis "Date"
  xlab("Date")+
  #Label the y-axis "Height (cm)"
  ylab("Average Height(cm)")

#Make plot from data frame avg_flowers in order to graph the height of each plant by diversity treatment with warming/not in two different colors#
ggplot(avg_leaves, aes(x=trtmnt_mono, y=num_leaves_avg, color=(warming)))+#filter out 9999
  geom_jitter()+
  #Add smooth line across data
  geom_smooth(method = "lm")+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Number of Leaves")

#### Harvest Data ####

#Make data frame merging harvest data with treatment data
harvest_trtmnt<-soybean_harvest%>%
  left_join(div)%>%
  mutate(trtmnt=paste(diversity,strains,sep="."))%>%
  mutate(trtmnt_mono=ifelse(diversity==0,0, ifelse(diversity==2,2,ifelse(diversity==3,3,trtmnt))))%>%
  select(-treatment_code,-strains,-Strain_1,-Strain_2,-Strain_3,-Strain_4,-trtmnt)%>%
  filter(total_pods!=9999)

  
#Total Pods  
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=total_pods, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Total Number Pods")

summary(Mixed_Model_Total_Pods<-lmer(total_pods~diversity*warming+(1|plant_num),data=harvest_trtmnt))
anova(Mixed_Model_Total_Pods)

#Green Pods
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=green_pods, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Green Pods")

#Brown Pods
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=brown_pods, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Brown Pods")

#Aborted Pods
ggplot(subset(harvest_trtmnt,aborted_pods!=9999), aes(x=trtmnt_mono, y=aborted_pods, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Aborted Pods")


#Total Beans
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=total_beans, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Total Beans")

#Healthy Beans
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=hlthy_beans, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Healthy Beans")

#Damaged Beans
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=dmg_beans, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Damaged Beans")

#Aborted Beans
ggplot(harvest_trtmnt, aes(x=trtmnt_mono, y=abrt_beans, color=(warming)))+#filter out 9999
  geom_boxplot()+
  #Label the x-axis "Date"
  xlab("Bacterial Diversity")+
  #Label the y-axis "Height (cm)"
  ylab("Aborted Beans")



