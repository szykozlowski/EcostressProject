library(readr)
library(tidyr)
library(lubridate)
library(ggplot2)

all_Swiz_WUE_data = read_csv("C:/Users/Szymon/Desktop/EcostressProject/Data/SwitzerlandPointData-2020-2022-ECO4WUE-001-results.csv")
colnames(all_Swiz_WUE_data)[5] = "Date_total"
colnames(all_Swiz_WUE_data)[9] = "WUEavg"


# Adjust time zone + seperate date and time
all_Swiz_WUE_data$date = as_datetime(all_Swiz_WUE_data$Date)
all_Swiz_WUE_data$date.adj = all_Swiz_WUE_data$date + hours(1)
all_Swiz_WUE_data = separate(data=all_Swiz_WUE_data, col=date.adj, into=c("date","time"), sep =" ")

# extract "m-d" from "Y-m-d"
n_last = 5
substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))

all_Swiz_WUE_data$m_d = substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))
# Make $ID character
all_Swiz_WUE_data$ID = as.character(all_Swiz_WUE_data$ID)

# seperate into seasonal
Swiz_WUE_Data_Summer_20 = all_Swiz_WUE_data[all_Swiz_WUE_data$Date >= "2020-01-01" & all_Swiz_WUE_data$Date <= "2020-12-31",]
Swiz_WUE_Data_Summer_21 = all_Swiz_WUE_data[all_Swiz_WUE_data$Date >= "2021-01-01" & all_Swiz_WUE_data$Date <= "2021-12-31",]
Swiz_WUE_Data_Summer_22 = all_Swiz_WUE_data[all_Swiz_WUE_data$Date >= "2022-01-01" & all_Swiz_WUE_data$Date <= "2022-12-31",]

Swiz_WUE_Data_Summer_20 <- Swiz_WUE_Data_Summer_20[order(Swiz_WUE_Data_Summer_20$Date_total),]
Swiz_WUE_Data_Summer_21 <- Swiz_WUE_Data_Summer_21[order(Swiz_WUE_Data_Summer_21$Date_total),]
Swiz_WUE_Data_Summer_22 <- Swiz_WUE_Data_Summer_22[order(Swiz_WUE_Data_Summer_22$Date_total),]
all_Swiz_WUE_data <- all_Swiz_WUE_data[order(all_Swiz_WUE_data$Date_total),]

# Read in RF_Temp_RH csv
Swiz_RF_Temp_RH_data = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Swiz_RF_Temp_RH/switz-sitelocations.csv", na.strings = "#N/A")

# Rename columns to "Latitude" and "Longitude"
colnames(Swiz_RF_Temp_RH_data)[3] = "Latitude"
colnames(Swiz_RF_Temp_RH_data)[4] = "Longitude"

# merge by lat and long
All_Swiz_WUE_RF_Temp_RH_data.2 = merge(all_Swiz_WUE_data,Swiz_RF_Temp_RH_data,by=c("Latitude","Longitude"))

All_Swiz_WUE_RF_Temp_RH_data = All_Swiz_WUE_RF_Temp_RH_data.2[All_Swiz_WUE_RF_Temp_RH_data.2$WUEavg<6,]

test = do.call(data.frame, aggregate(All_Swiz_WUE_RF_Temp_RH_data$WUEavg, list(All_Swiz_WUE_RF_Temp_RH_data$ID,All_Swiz_WUE_RF_Temp_RH_data$mean.annual.P, All_Swiz_WUE_RF_Temp_RH_data$mean.annual.T, All_Swiz_WUE_RF_Temp_RH_data$mean.annual.RH), 
                                     FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))))

names(test) = c("ID", "mean.annual.P", "mean.annual.T","mean.annual.RH", "WUEavg", "WUEavg.sd")
par(mfrow = c(3,1))

# All Summers - UNSURE WHAT GOES INTO Y AND COLOR
Swiz_loess_WUE_All = ggplot(All_Swiz_WUE_RF_Temp_RH_data, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("All Summers") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))  +
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2020 WUE
Swiz_loess_WUE_Summer_20 = ggplot(Swiz_WUE_Data_Summer_20, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2018") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2021 WUE
Swiz_loess_WUE_Summer_21 = ggplot(Swiz_WUE_Data_Summer_21, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2019") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2022 WUE
Swiz_loess_WUE_Summer_22 = ggplot(Swiz_WUE_Data_Summer_22, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2020") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))


# Plot all of the WUE Data for Every Summer
plot(Swiz_loess_WUE_All)

# Plot WUE Data for Indivudual Summers
plot(Swiz_loess_WUE_Summer_20)
plot(Swiz_loess_WUE_Summer_21)
plot(Swiz_loess_WUE_Summer_22)


plot(data = test, WUEavg ~ mean.annual.P)
plot(data = test, WUEavg ~ mean.annual.T)
plot(data = test, WUEavg ~ mean.annual.RH)