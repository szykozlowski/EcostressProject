install.packages("reshape2")

library(raster)
library(reshape2)
library(rgdal)

######
# Add country outline: 
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Swiz_Shapefile/Switzerland_shapefile")
switz = readOGR("ch_1km.shp")
# reproject spatial outline to match raster data
switz_RP = spTransform(switz, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

######
# Add mask:
forestmask = readOGR("C:/Users/Szymon/Desktop/EcostressProject/Data/Swiz_Forest_Cover_Shapefile/Vector_Landuse_CH/VEC200_LandCover.shp")

forestmask_RP = spTransform(forestmask,crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


#######################################################################

#https://strimas.com/post/processing-large-rasters-in-r/

#/Users/goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF

# JULY 2018 TIF

# Set working directory to local source
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
library(rstudioapi)




# setwd("~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF")
# Swiz_rastlist_JULY_18 <- list.files(path = setwd("~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
# allrasters_Swiz_JULY_18 <- lapply(Swiz_rastlist_JULY_18, raster)
# allrasters_Swiz_JULY_18.n <- lapply(allrasters_Swiz_JULY_18,FUN=crop,switz_RP)

# Felton's code

# Set working directory to local source
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))
library(rstudioapi)
library(raster)

test<-"./../Switz_WUE_July_2018_TIF/"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop
	
	}


# end Felton's code
	


###allrasters_Swiz_JULY_18.n2 <- lapply(allrasters_Swiz_JULY_18,FUN=mask,forestmask_RP)


# AUG 2018 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2018_TIF")
Swiz_rastlist_AUG_18 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2018_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_AUG_18 <- lapply(Swiz_rastlist_AUG_18, raster)

# JUNE 2019 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_June_2019_TIF")
Swiz_rastlist_JUNE_19 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_June_2019_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_JUNE_19 <- lapply(Swiz_rastlist_JUNE_19, raster)

# JULY 2019 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_July_2019_TIF")
Swiz_rastlist_JULY_19 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_July_2019_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_JULY_19 <- lapply(Swiz_rastlist_JULY_19, raster)

# AUG 2019 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2019_TIF")
Swiz_rastlist_AUG_19 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2019_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_AUG_19 <- lapply(Swiz_rastlist_AUG_19, raster)

# JUNE 2020 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_June_2020_TIF")
Swiz_rastlist_JUNE_20 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_June_2020_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_JUNE_20 <- lapply(Swiz_rastlist_JUNE_20, raster)

# JULY 2020 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_July_2020_TIF")
Swiz_rastlist_JULY_20 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_July_2020_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_JULY_20 <- lapply(Swiz_rastlist_JULY_20, raster)

# AUG 2020 TIF
setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2020_TIF")
Swiz_rastlist_AUG_20 <- list.files(path = setwd("C:/Users/Szymon/Desktop/EcostressProject/Data/Switz_WUE_Aug_2020_TIF"),pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters_Swiz_AUG_20 <- lapply(Swiz_rastlist_AUG_20, raster)

# Summer 2018 #######################

# Extend July 2018
# # large_dims = extent(5, 11, 45.5, 48)
# # extend_whole = extend(allrasters_Swiz_JULY_18[[3]], large_dims)
# # extended_allrasters_Swiz_JULY_18 = lapply(allrasters_Swiz_JULY_18, resample, extend_whole, method = "bilinear")
# # stack_Swiz_JULY_18 = stack(extended_allrasters_Swiz_JULY_18)
# # stack_Swiz_JULY_18.agg = aggregate(stack_Swiz_JULY_18,fact=20)

# Extend July 2018
large_dims = extent(5, 11, 45.5, 48)
dummy.raster=raster(ext=large_dims,res=c(.05,.05))
extended_allrasters_Swiz_JULY_18 = lapply(allrasters_Swiz_JULY_18, resample, dummy.raster, method = "bilinear")
stack_Swiz_JULY_18 = stack(extended_allrasters_Swiz_JULY_18)
stack_Swiz_JULY_18.agg = aggregate(stack_Swiz_JULY_18,fact=20)

# Extend AUG 2018
extended_allrasters_Swiz_AUG_18 = lapply(allrasters_Swiz_AUG_18, resample, dummy.raster, method = "bilinear")
stack_Swiz_AUG_18 = stack(extended_allrasters_Swiz_AUG_18,quick=TRUE)
stack_Swiz_AUG_18.agg = aggregate(stack_Swiz_AUG_18,fact=20)


# All Summer 18 Stack
stack_Swiz_Summer_18 = stack(stack_Swiz_JULY_18.agg,stack_Swiz_AUG_18.agg)
stack_Swiz_Summer_18_sd = calc(stack_Swiz_Summer_18, fun = sd)

sd(stack_Swiz_Summer_18, na.rm = TRUE)

# Plot Summer 18 WUE + SD
#plot(mean(stack_Swiz_Summer_18, na.rm = TRUE), col= brewer.pal(9,"RdYlBu"))
#plot(mean(stack_Swiz_Summer_18_sd, na.rm = TRUE), col= brewer.pal(9,"RdYlBu"))

mean_stack_Swiz_Summer_18 = mean(stack_Swiz_Summer_18, na.rm = TRUE)
plot(mean_stack_Swiz_Summer_18, col= brewer.pal(9,"RdYlBu"))


plot(stack_Swiz_Summer_18_sd, na.rm = TRUE, col= brewer.pal(9,"RdYlBu"))


# Summer 2019 #######################
# Extend June 2019
extended_allrasters_Swiz_JUNE_19 = lapply(allrasters_Swiz_JUNE_19, resample, extend_whole, method = "bilinear")
stack_Swiz_JUNE_19 = stack(extended_allrasters_Swiz_JUNE_19)

# Extend JULY 2019
extended_allrasters_Swiz_JULY_19 = lapply(allrasters_Swiz_JULY_19, resample, extend_whole, method = "bilinear")
stack_Swiz_JULY_19 = stack(extended_allrasters_Swiz_JULY_19)

# Extend AUG 2019
extended_allrasters_Swiz_AUG_19 = lapply(allrasters_Swiz_AUG_19, resample, extend_whole, method = "bilinear")
stack_Swiz_AUG_19 = stack(extended_allrasters_Swiz_AUG_19)

# All Summer 19 Stack + sd
stack_Swiz_Summer_19 = stack(stack_Swiz_JUNE_19, stack_Swiz_JULY_19, stack_Swiz_AUG_19)
stack_Swiz_Summer_19_sd = calc(stack_Swiz_Summer_19, sd)

# Plot Summer 19 WUE
#plot(mean(stack_Swiz_Summer_19, na.rm = TRUE), col= brewer.pal(9,"RdYlBu"))

mean_stack_Swiz_Summer_19 = mean(stack_Swiz_Summer_19, na.rm = TRUE)
plot(mean_stack_Swiz_Summer_19, col= brewer.pal(9,"RdYlBu"))


# Summer 2020 #######################
# Extend June 2020
extended_allrasters_Swiz_JUNE_20 = lapply(allrasters_Swiz_JUNE_20, resample, extend_whole, method = "bilinear")
stack_Swiz_JUNE_20 = stack(extended_allrasters_Swiz_JUNE_20)

# Extend JULY 2020
extended_allrasters_Swiz_JULY_20 = lapply(allrasters_Swiz_JULY_20, resample, extend_whole, method = "bilinear")
stack_Swiz_JULY_20 = stack(extended_allrasters_Swiz_JULY_20)

# Extend AUG 2020
extended_allrasters_Swiz_AUG_20 = lapply(allrasters_Swiz_AUG_20, resample, extend_whole, method = "bilinear")
stack_Swiz_AUG_20 = stack(extended_allrasters_Swiz_AUG_20)

# All Summer 20 Stack + sd
stack_Swiz_Summer_20 = stack(stack_Swiz_JUNE_20, stack_Swiz_JULY_20, stack_Swiz_AUG_20)
stack_Swiz_Summer_20_sd = calc(stack_Swiz_Summer_20, sd)

# Plot All Summer 20 WUE + SD? - how to plot SD?
plot(mean(stack_Swiz_Summer_20, na.rm = TRUE), col= brewer.pal(9,"RdYlBu"))

mean_stack_Swiz_Summer_20 = mean(stack_Swiz_Summer_20, na.rm = TRUE)
plot(mean_stack_Swiz_Summer_20, col= brewer.pal(9,"RdYlBu"))


# mask images - on hold
setwd("~/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Swiz_Forest_Cover_Shapefile/Vector_Landuse_CH")
shpfile = readOGR("VEC200_LandCover.shp",layer="VEC200_LandCover")
plot(shpfile)


#mask(stack_Swiz_Summer_20, )

#We need to remake the loess plot with lines for each site colored by rainfall. 
#And then remake it by sites based on their species. So, can you plot lines by site 
#but color by mean annual rainfall? I sent that to you in a separate file. 

# WUE ####################

# read in csv files
WUE_June_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_June_2019/June-2019-ECO4WUE-001-results.csv")
WUE_June_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_June_2020/June-2020-ECO4WUE-001-results.csv")
WUE_July_18 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2018/July-2018-ECO4WUE-001-results.csv")
WUE_July_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2019/July-2019-ECO4WUE-001-results.csv")
WUE_July_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_July_2020/July-2020-ECO4WUE-001-results.csv")
WUE_Aug_18 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2018/August-2018-ECO4WUE-001-results.csv")
WUE_Aug_19 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2019/August-2019-ECO4WUE-001-results.csv")
WUE_Aug_20 = read_csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_Aug_2020/August-2020-ECO4WUE-001-results.csv")

# rbind
all_Swiz_WUE_data = rbind(WUE_June_19, WUE_June_20, WUE_July_18, WUE_July_19, WUE_July_20,
                          WUE_Aug_18, WUE_Aug_19,WUE_Aug_20)

colnames(all_Swiz_WUE_data)[5] = "Date_total"
colnames(all_Swiz_WUE_data)[9] = "WUEavg"

# Adjust time zone + seperate date and time
all_Swiz_WUE_data$date = as_datetime(all_Swiz_WUE_data$Date_total)
all_Swiz_WUE_data$date.adj = all_Swiz_WUE_data$date + hours(1)
all_Swiz_WUE_data = separate(data=all_Swiz_WUE_data, col=date.adj, into=c("date","time"), sep =" ")

# extract "m-d" from "Y-m-d"
n_last = 5
substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))
all_Swiz_WUE_data$m_d = substr(all_Swiz_WUE_data$date, nchar(all_Swiz_WUE_data$date) - n_last + 1, nchar(all_Swiz_WUE_data$date))

# Make $ID character
all_Swiz_WUE_data$ID = as.character(all_Swiz_WUE_data$ID)

# seperate into seasonal
Swiz_WUE_Data_Summer_18 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2018-01-01" & all_Swiz_WUE_data$date <= "2018-12-31",]
Swiz_WUE_Data_Summer_19 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2019-01-01" & all_Swiz_WUE_data$date <= "2019-12-31",]
Swiz_WUE_Data_Summer_20 = all_Swiz_WUE_data[all_Swiz_WUE_data$date >= "2020-01-01" & all_Swiz_WUE_data$date <= "2020-12-31",]

# LOESS graphs

# Read in RF_Temp_RH csv
Swiz_RF_Temp_RH_data = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Swiz_RF_Temp_RH/switz-sitelocations.csv", na.strings = "#N/A")

# Rename columns to "Latitude" and "Longitude"
colnames(Swiz_RF_Temp_RH_data)[3] = "Latitude"
colnames(Swiz_RF_Temp_RH_data)[4] = "Longitude"

# merge by lat and long
All_Swiz_WUE_RF_Temp_RH_data.2 = merge(all_Swiz_WUE_data,Swiz_RF_Temp_RH_data,by=c("Latitude","Longitude"))
head(All_Swiz_WUE_RF_Temp_RH_data.2)

All_Swiz_WUE_RF_Temp_RH_data = All_Swiz_WUE_RF_Temp_RH_data.2[All_Swiz_WUE_RF_Temp_RH_data.2$WUEavg<6,]
str(All_Swiz_WUE_RF_Temp_RH_data)
max(All_Swiz_WUE_RF_Temp_RH_data$WUEavg)

nrow(All_Swiz_WUE_RF_Temp_RH_data)

test = do.call(data.frame, aggregate(All_Swiz_WUE_RF_Temp_RH_data$WUEavg, list(All_Swiz_WUE_RF_Temp_RH_data$ID,All_Swiz_WUE_RF_Temp_RH_data$mean.annual.P, All_Swiz_WUE_RF_Temp_RH_data$mean.annual.T, All_Swiz_WUE_RF_Temp_RH_data$mean.annual.RH), 
                 FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))))
                 
names(test) = c("ID", "mean.annual.P", "mean.annual.T","mean.annual.RH", "WUEavg", "WUEavg.sd")
str(test)
par(mfrow = c(3,1))
plot(data = test, WUEavg ~ mean.annual.P)
plot(data = test, WUEavg ~ mean.annual.T)
plot(data = test, WUEavg ~ mean.annual.RH)

# All Summers - UNSURE WHAT GOES INTO Y AND COLOR
Swiz_loess_WUE_All = ggplot(All_Swiz_WUE_RF_Temp_RH_data, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("All Summers") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))  +
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2018 WUE
Swiz_loess_WUE_Summer_18 = ggplot(Swiz_WUE_Data_Summer_18, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2018") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2019 WUE
Swiz_loess_WUE_Summer_19 = ggplot(Swiz_WUE_Data_Summer_19, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2019") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Summer 2020 WUE
Swiz_loess_WUE_Summer_20 = ggplot(Swiz_WUE_Data_Summer_20, aes(x = as.POSIXct(m_d,format = "%m-%d"), y = WUEavg, color = ID)) + 
  geom_point() + geom_smooth(se = FALSE) + theme_bw() + ggtitle("Summer 2020") +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) + 
  xlab("Date") + ylab(expression(paste("WUE (g C kg"^"-1 ",H[2],"O)")))

# Adjust WUE

# read in ET data
ET_June_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_June_2019/ET-June-2019-ECO3ETPTJPL-001-results.csv")
ET_June_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_June_2020/ET-June-2020-ECO3ETPTJPL-001-results.csv")
ET_July_18 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2018/ET-July-2018-ECO3ETPTJPL-001-results.csv")
ET_July_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2019/ET-July-2019-ECO3ETPTJPL-001-results.csv")
ET_July_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_July_2020/ET-July-2020-ECO3ETPTJPL-001-results.csv")
ET_Aug_18 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2018/ET-Aug-2018-ECO3ETPTJPL-001-results.csv")
ET_Aug_19 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2019/ET-Aug-2019-ECO3ETPTJPL-001-results.csv")
ET_Aug_20 = read.csv("/Users/brandonbernardo/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_ET_Aug_2020/ET-Aug-2020-ECO3ETPTJPL-001-results.csv")

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
Plot_Swiz_Both_WUE = ggarrange(Plot_Swiz_WUE_GPP_ET, Plot_Swiz_WUE_GPP_T,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)





