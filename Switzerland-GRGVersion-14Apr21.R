install.packages("reshape2")

library(scico)
library(ggplot2)
library(sf)
library(gridExtra)
library(raster)
library(reshape2)
library(rgdal)

#--------------------------------------------------------------#
# Set a WD (Goldsmith Approach)
#--------------------------------------------------------------#

setwd("~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/")

#--------------------------------------------------------------#
# Upload
#--------------------------------------------------------------#

########JULY 2018

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2018_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop
	}
	
large_dims = extent(5, 11, 45.5, 48)
dummy.raster=raster(ext=large_dims)
extended_allrasters_Swiz_July_18 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_July_18 = stack(extended_allrasters_Swiz_July_18)
writeRaster(stack_Swiz_July_18,"stack_Swiz_July_18.tif")

########AUGUST 2018

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2018_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_Aug_18 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_Aug_18 = stack(extended_allrasters_Swiz_Aug_18)
writeRaster(stack_Swiz_Aug_18,"stack_Swiz_Aug_18.tif")

########JUNE 2019

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2019_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_June_19 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_June_19 = stack(extended_allrasters_Swiz_June_19)
writeRaster(stack_Swiz_June_19,"stack_Swiz_June_19.tif")

########JULY 2019

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2019_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_July_19 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_July_19 = stack(extended_allrasters_Swiz_July_19)
writeRaster(stack_Swiz_July_19,"stack_Swiz_July_19.tif")

########AUG 2019

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2019_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_Aug_19 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_Aug_19 = stack(extended_allrasters_Swiz_Aug_19)
writeRaster(stack_Swiz_Aug_19,"stack_Swiz_Aug_19.tif")

########JUNE 2020

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_June_2020_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_June_20 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_June_20 = stack(extended_allrasters_Swiz_June_20)
writeRaster(stack_Swiz_June_20,"stack_Swiz_June_20.tif")

########JULY 2020

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_July_2020_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_July_20 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_July_20 = stack(extended_allrasters_Swiz_July_20)
writeRaster(stack_Swiz_July_20,"stack_Swiz_July_20.tif")

########AUG 2020

test<-"~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Switz_WUE_Aug_2020_TIF"
test.dir<-dir(test,full.names = T)

loop.list <- list()

for(i in test.dir){
  
  test<-raster(i)
	test.loop = aggregate(test,fact=200)
	loop.list[[i]] = test.loop	
	}
	
extended_allrasters_Swiz_Aug_20 = lapply(loop.list, resample, dummy.raster, method = "bilinear")
stack_Swiz_Aug_20 = stack(extended_allrasters_Swiz_Aug_20)
writeRaster(stack_Swiz_Aug_20,"stack_Swiz_Aug_20.tif")

#--------------------------------------------------------------#
# Add country outline and DEM 
#--------------------------------------------------------------#

setwd("~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Swiz_Shapefile/")
switz = readOGR("G1L12.shp",layer="G1L12")
# reproject spatial outline to match raster data
switz_RP = spTransform(switz, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

dem = raster("DHM200.asc")
proj21781 = "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs"
projection(dem) = proj21781
dem_RP = spTransform(dem, crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

dem=readOGR("dhm25_l.shp")

#--------------------------------------------------------------#
# Add forest mask:
#--------------------------------------------------------------#

forestmask = readOGR("~goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/Bernardo_All_Swiz_Data/Switzerland copy/Swiz_Forest_Cover_Shapefile/Vector_Landuse_CH/VEC200_LandCover.shp")
forestmask_RP = spTransform(forestmask,crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#--------------------------------------------------------------#
#All summer plots: 
#--------------------------------------------------------------#

stack_Swiz1=stack(stack_Swiz_July_18,stack_Swiz_Aug_18,stack_Swiz_June_19,stack_Swiz_July_19,stack_Swiz_Aug_19,stack_Swiz_June_20,stack_Swiz_July_20,stack_Swiz_Aug_20)
stack_Swiz=mask(stack_Swiz1,forestmask_RP,inverse=TRUE)
stack_Swiz.mean = calc(stack_Swiz,fun=median,na.rm=TRUE)
stack_Swiz.sd = calc(stack_Swiz,fun=sd,na.rm=TRUE)

stack_Swiz_mean_pts = rasterToPoints(stack_Swiz.mean,spatial=TRUE)
stack_Swiz_df = data.frame(stack_Swiz_mean_pts)

Swiz_mean = ggplot()+
	geom_raster(data = stack_Swiz_df,aes(x = x, y = y, fill = layer))+
	scale_fill_scico(palette = 'batlow',direction=-1)+
	labs(fill="WUE")+
		ggtitle("Median summer WUE")+
	geom_polygon(data=switz_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
	theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+ 
	coord_quickmap()

#####Standard deviation 

stack_Swiz_sd_pts = rasterToPoints(stack_Swiz.sd,spatial=TRUE)
stack_Swiz_sd_df = data.frame(stack_Swiz_sd_pts)

Swiz_sd = ggplot()+
	geom_raster(data = stack_Swiz_sd_df,aes(x = x, y = y, fill = layer))+
	scale_fill_scico(palette = 'batlow',direction=-1)+
	labs(fill="WUE")+
	ggtitle("Standard deviation of summer WUE")+
	geom_polygon(data=switz_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
	theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+ 
	coord_quickmap()

grid.arrange(Swiz_mean,Swiz_sd,ncol=2)

#--------------------------------------------------------------#
#Each summer plots: 
#--------------------------------------------------------------#

summer2018.all=stack(stack_Swiz_July_18,stack_Swiz_Aug_18)
summer2018=mask(summer2018.all,forestmask_RP,inverse=TRUE)
summer2019.all=stack(stack_Swiz_June_19,stack_Swiz_July_19,stack_Swiz_Aug_19)
summer2019=mask(summer2019.all,forestmask_RP,inverse=TRUE)
summer2020.all=stack(stack_Swiz_June_20,stack_Swiz_July_20,stack_Swiz_Aug_20)
summer2020=mask(summer2020.all,forestmask_RP,inverse=TRUE)

summer2018.mean = calc(summer2018,fun=median,na.rm=TRUE)
summer2019.mean = calc(summer2019,fun=median,na.rm=TRUE)
summer2020.mean = calc(summer2020,fun=median,na.rm=TRUE)

summer2018_pts = rasterToPoints(summer2018.mean,spatial=TRUE)
summer2018_df = data.frame(summer2018_pts)

summer2019_pts = rasterToPoints(summer2019.mean,spatial=TRUE)
summer2019_df = data.frame(summer2019_pts)

summer2020_pts = rasterToPoints(summer2020.mean,spatial=TRUE)
summer2020_df = data.frame(summer2020_pts)

summer2018.plot = ggplot()+
	geom_raster(data = summer2018_df,aes(x = x, y = y, fill = layer))+
	scale_fill_scico(palette = 'batlow',limits=c(0,3.5),direction=-1)+
	labs(fill="WUE")+
	geom_polygon(data=switz_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
	ggtitle("Summer 2018")+
	theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+ 
	coord_quickmap()
	
summer2019.plot = ggplot()+
	geom_raster(data = summer2019_df,aes(x = x, y = y, fill = layer))+
	scale_fill_scico(palette = 'batlow',limits=c(0,3.5),direction=-1)+
	labs(fill="WUE")+
	geom_polygon(data=switz_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
	ggtitle("Summer 2019")+
	theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+ 
	coord_quickmap()
	
summer2020.plot = ggplot()+
	geom_raster(data = summer2020_df,aes(x = x, y = y, fill = layer))+
	scale_fill_scico(palette = 'batlow',limits=c(0,3.5),direction=-1,alpha=0.5)+
	labs(fill="WUE")+
	ggtitle("Summer 2020")+
	geom_polygon(data=switz_RP,aes(x=long, y=lat, group=group),alpha=0,color="black")+
	theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+ 
	coord_quickmap()
	
grid.arrange(summer2018.plot,summer2019.plot,summer2020.plot,ncol=3)

        	####annotation_scale(data=summer2020_df,plot_unit="m",location="br")+


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





