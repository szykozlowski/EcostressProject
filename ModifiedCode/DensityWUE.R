library(readr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(dplyr)

# read in ET data
ET_June_19 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_June_2019/ET-June-2019-ECO3ETPTJPL-001-results.csv")
ET_June_20 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_June_2020/ET-June-2020-ECO3ETPTJPL-001-results.csv")
ET_July_18 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_July_2018/ET-July-2018-ECO3ETPTJPL-001-results.csv")
ET_July_19 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_July_2019/ET-July-2019-ECO3ETPTJPL-001-results.csv")
ET_July_20 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_July_2020/ET-July-2020-ECO3ETPTJPL-001-results.csv")
ET_Aug_18 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_Aug_2018/ET-Aug-2018-ECO3ETPTJPL-001-results.csv")
ET_Aug_19 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_Aug_2019/ET-Aug-2019-ECO3ETPTJPL-001-results.csv")
ET_Aug_20 = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_ET_Aug_2020/ET-Aug-2020-ECO3ETPTJPL-001-results.csv")

# rbind
all_Swiz_ET_data = rbind(ET_June_19, ET_June_20, ET_July_18, ET_July_19, ET_July_20,
                         ET_Aug_18, ET_Aug_19, ET_Aug_20)

colnames(all_Swiz_ET_data)[5] = "Date_total"
colnames(all_Swiz_ET_data)[9] = "ETcanopy"
colnames(all_Swiz_ET_data)[10] = "ETdaily"

# Adjust time zone + seperate date and time
all_Swiz_ET_data$date = as_datetime(all_Swiz_ET_data$Date_total)
all_Swiz_ET_data$date.adj = all_Swiz_ET_data$date + hours(1)
all_Swiz_ET_data = separate(data=all_Swiz_ET_data, col=date.adj, into=c("date","time"), sep =" ")

# extract "m-d" from "Y-m-d"
n_last = 5
substr(all_Swiz_ET_data$date, nchar(all_Swiz_ET_data$date) - n_last + 1, nchar(all_Swiz_ET_data$date))
all_Swiz_ET_data$m_d = substr(all_Swiz_ET_data$date, nchar(all_Swiz_ET_data$date) - n_last + 1, nchar(all_Swiz_ET_data$date))

# Make $ID character
all_Swiz_ET_data$ID = as.character(all_Swiz_ET_data$ID)

# seperate into seasonal
Swiz_ET_Data_Summer_18 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2018-01-01" & all_Swiz_ET_data$date <= "2018-12-31",]
Swiz_ET_Data_Summer_19 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2019-01-01" & all_Swiz_ET_data$date <= "2019-12-31",]
Swiz_ET_Data_Summer_20 = all_Swiz_ET_data[all_Swiz_ET_data$date >= "2020-01-01" & all_Swiz_ET_data$date <= "2020-12-31",]


# Calculate GPP 

# make coulumns integers to match ET column
all_Swiz_WUE_data$`Orbit Number` = as.integer(all_Swiz_WUE_data$`Orbit Number`)
all_Swiz_WUE_data$`Build ID` = as.integer(all_Swiz_WUE_data$`Build ID`)
all_Swiz_WUE_data$`Scene ID` = as.integer(all_Swiz_WUE_data$`Scene ID`)

# make unique ID by pasting all 8 columns together
all_Swiz_WUE_data$Unique = paste(all_Swiz_WUE_data$Category, all_Swiz_WUE_data$ID, all_Swiz_WUE_data$Latitude,
                                 all_Swiz_WUE_data$Longitude, all_Swiz_WUE_data$Date_total, all_Swiz_WUE_data$`Orbit Number`,
                                 all_Swiz_WUE_data$`Build ID`, all_Swiz_WUE_data$`Scene ID`)

all_Swiz_ET_data$Unique = paste(all_Swiz_ET_data$Category, all_Swiz_ET_data$ID, all_Swiz_ET_data$Latitude,
                                all_Swiz_ET_data$Longitude, all_Swiz_ET_data$Date_total, all_Swiz_ET_data$Orbit.Number,
                                all_Swiz_ET_data$Build.ID, all_Swiz_ET_data$Scene.ID)

# only need these three columns from ET to make 1 big data fram
Essential_Swiz_ET_data = all_Swiz_ET_data[, c("ETcanopy", "ETdaily", "Unique")]

# merge, final data fram w all data and sites matched up via $unique 
All_Swiz_WUE_ET_Data = merge(all_Swiz_WUE_data, Essential_Swiz_ET_data)

# WUE (GPP/ET) - no calcs necessary 
head(All_Swiz_WUE_ET_Data$WUEavg)

# WUE (GPP/T) = (WUE x ET Daily) / (ET Daily x ET canopy/100)
All_Swiz_WUE_ET_Data$WUE_GPP_by_T = (All_Swiz_WUE_ET_Data$WUEavg * All_Swiz_WUE_ET_Data$ETdaily) / (All_Swiz_WUE_ET_Data$ETdaily * All_Swiz_WUE_ET_Data$ETcanopy/100)

# Remove inf values
All_Swiz_WUE_ET_Data = All_Swiz_WUE_ET_Data %>% 
  filter_if(~is.numeric(.), all_vars(!is.infinite(.)))

# remove values greater than 10
All_Swiz_WUE_ET_Data = All_Swiz_WUE_ET_Data %>% 
  filter(WUE_GPP_by_T <= 5)

# ggplot geom_density histogram of WUE (GPP/ET) compared to WUE (GPP/T)

# WUE (GPP/T)
Plot_Swiz_WUE_GPP_T = ggplot(All_Swiz_WUE_ET_Data, aes(x=WUE_GPP_by_T)) + geom_histogram(aes(y=..density..), binwidth=.05,color="#87CEFA", fill="white") +
  geom_density(alpha=.2, fill="#87CEFA") + xlim(0,4) + ylim(0,1.25) + 
  geom_vline(aes(xintercept=mean(WUE_GPP_by_T)), color="blue", linetype="dashed", size=1) +
  xlab("WUE (GPP/T)")

# WUE (GPP/ET)
Plot_Swiz_WUE_GPP_ET = ggplot(All_Swiz_WUE_ET_Data, aes(x=WUEavg)) + geom_histogram(aes(y=..density..), binwidth=.05,color="#F08080", fill="white") +
  geom_density(alpha=.2, fill="#F08080") + xlim(0,4) + ylim(0,1.25) + 
  geom_vline(aes(xintercept=mean(WUEavg)), color="red", linetype="dashed", size=1) +
  xlab("WUE (GPP/ET)")

# 2 Panel Image
# Plot_Swiz_Both_WUE = ggarrange(Plot_Swiz_WUE_GPP_ET, Plot_Swiz_WUE_GPP_T,
                               # labels = c("A", "B"),
                               # ncol = 1, nrow = 2)
plot(Plot_Swiz_WUE_GPP_T)
plot(Plot_Swiz_WUE_GPP_ET)