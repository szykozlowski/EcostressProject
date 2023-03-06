library(readr)
library(tidyr)
library(lubridate)
library(ggplot2)


# Adjust time zone + seperate date and time
Swiz_WUE_Data_Summer_19_csv = read_csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Ecostress2018-ECO4WUE-001-results.csv")
colnames(Swiz_WUE_Data_Summer_19_csv)[5] = "Date_total"
colnames(Swiz_WUE_Data_Summer_19_csv)[9] = "WUEavg"

# extract "m-d" from "Y-m-d" 
n_last = 5
substr(Swiz_WUE_Data_Summer_19_csv$date, nchar(Swiz_WUE_Data_Summer_19_csv$date) - n_last + 1, nchar(Swiz_WUE_Data_Summer_19_csv$date))

# AUGHHHH
Swiz_WUE_Data_Summer_19_csv$m_d = substr(Swiz_WUE_Data_Summer_19_csv$date, nchar(Swiz_WUE_Data_Summer_19_csv$date) - n_last + 1, nchar(Swiz_WUE_Data_Summer_19_csv$date))
all_Swiz_WUE_data$m_d = substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))

# Sort and Date
Swiz_WUE_Data_Summer_19_csv$ID = as.character(Swiz_WUE_Data_Summer_19_csv$ID)
Swiz_WUE_Data_Summer_19 = Swiz_WUE_Data_Summer_19_csv[Swiz_WUE_Data_Summer_19_csv$Date >= "2018-06-01" & Swiz_WUE_Data_Summer_19_csv$Date <= "2018-09-31",]
Swiz_WUE_Data_Summer_19 <- Swiz_WUE_Data_Summer_19[order(Swiz_WUE_Data_Summer_19$Date_total),]


Swiz_RF_Temp_RH_data = read.csv("C:/Users/Szymon/Desktop/EcostressProject/Data/Swiz_RF_Temp_RH/switz-sitelocations.csv", na.strings = "#N/A")

# Rename columns to "Latitude" and "Longitude"
colnames(Swiz_RF_Temp_RH_data)[3] = "Latitude"
colnames(Swiz_RF_Temp_RH_data)[4] = "Longitude"

# Remove values greater than 6
Swiz_WUE_Data_Summer_19_data.2 = merge(Swiz_WUE_Data_Summer_19,Swiz_RF_Temp_RH_data,by=c("Latitude","Longitude"))
All_Swiz_WUE_RF_Temp_RH_data = All_Swiz_WUE_RF_Temp_RH_data.2[All_Swiz_WUE_RF_Temp_RH_data.2$WUEavg<6,]
All_Swiz_WUE_RF_Temp_RH_data_Summer = All_Swiz_WUE_RF_Temp_RH_data_Summer.2[All_Swiz_WUE_RF_Temp_RH_data_Summer.2$WUEavg<6,]

#Single Year
Swiz_WUE_Data_Summer_19_data = Swiz_WUE_Data_Summer_19_data.2[Swiz_WUE_Data_Summer_19_data.2$WUEavg<6,]

# Single Year
test3 = do.call(data.frame, aggregate(Swiz_WUE_Data_Summer_19_data$WUEavg, list(Swiz_WUE_Data_Summer_19_data$ID,Swiz_WUE_Data_Summer_19_data$mean.annual.P, Swiz_WUE_Data_Summer_19_data$mean.annual.T, Swiz_WUE_Data_Summer_19_data$mean.annual.RH), 
                                      FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))))

# Summer 2019 WUE
Swiz_loess_WUE_Summer_19 = ggplot(Swiz_WUE_Data_Summer_19_data, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2018") + ylim(0,12) +
  
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

plot(Swiz_loess_WUE_Summer_19)


