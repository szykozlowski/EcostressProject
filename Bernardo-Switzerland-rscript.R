setwd("/Users/goldsmit/Dropbox/NASA-ECOSTRESS/BernardoProject/Switzerland/")


#------------------------------------------------------#
#13C Data 
#------------------------------------------------------#

switz13C = read.csv("IAP_Summer2019_CN.csv",header=TRUE)
str(switz13C)

agg.switz13C = do.call(data.frame,aggregate(switz13C$Delta_13C,list(switz13C$siteID),FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),n=length(x))))
names(agg.switz13C) = c("siteID","mean13C","sd13C","n13C")
agg.switz13C

###Convert 13C to ci/ca: 

agg.switz13C$big.delta = ((-8.6-agg.switz13C$mean13C)/(1000+agg.switz13C$mean13C)*1000)
agg.switz13C$ci.ca = (((-8.6-agg.switz13C$mean13C)/(1000+agg.switz13C$mean13C)*1000)-4.4)/(27-4.4)
agg.switz13C$a.gs = (400*(1-agg.switz13C$ci.ca))/1.6

write.csv(agg.switz13C,"carbon13validationdata.csv")

#------------------------------------------------------#
#WUE ECOSTRESS DATA
#------------------------------------------------------#

####No ECOSTRESS data for site 1066? 

wue = read.csv("Ecostress.switzerland.csv",header=TRUE,na.strings=c("#DIV/0!","#N/A"))
str(wue)
wue.Aug19 = wue[which(wue$Year == "2019"),]
str(wue.Aug19)

agg.wue = do.call(data.frame,aggregate(wue.Aug19$GPP.T,list(wue.Aug19$ID),FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),n=length(x))))
names(agg.wue) = c("siteID","meanGPP.T","sdGPP.T","nGPP.T")

###Cribbed together data on June/July VPD and scaled to generate iWUE: 
agg.wue$iWUE=agg.wue$meanGPP.T*agg.wue$VPD

merged = merge(agg.wue,agg.switz13C)
str(merged)

plot(data=merged,meanGPP.T~ci.ca,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)))
summary(lm(data=merged,meanGPP.T~ci.ca))

wue.2 = wue[wue$GPP.T <8 & wue$GPP.T >0.2,]

agg.wue = do.call(data.frame,aggregate(wue.2$GPP.T,list(wue.2$ID),FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),min=min(x,na.rm=TRUE),max=max(x,na.rm=TRUE),n=length(x))))
names(agg.wue) = c("siteID","meanGPP.T","sdGPP.T","minGPP.T","maxGPP.T","nGPP.T")
write.csv(agg.wue,"agg.wue.csv")


#------------------------------------------------------#
#Basic Relationships
#------------------------------------------------------#

data.raw = read.csv("ECOSTRESSdata.switzerland-sitelocations-17Feb21.csv",header=TRUE,na.strings="#N/A")
data = data.raw[which(data.raw$meanGPP.T < 10),]

par(mfrow=c(3,3),mar=c(5,4.5,1.5,1),cex=1.02)
plot(meanGPP.T~elevation,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Elevation (m asl)")
mtext("(A)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~slope.degree,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Slope (\u00B0)")
mtext("(B)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~aspect.degree,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Aspect (\u00B0)")
mtext("(C)",side=3,line=-1.1,adj=.95)

plot(meanGPP.T~mean.annual.T,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Mean annual temperature (\u00B0C)")
mtext("(D)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~mean.annual.RH,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Mean annual relative humidity (%)")
mtext("(E)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~mean.annual.P,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Mean annual precipitation (mm)")
mtext("(F)",side=3,line=-1.1,adj=.95)

plot(meanGPP.T~total.n.deposition,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab=expression(paste("Total nitrogen deposition (kg N ",ha^-1, yr^-1,")")))
mtext("(G)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~oxidized.nitrogen,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab=expression(paste(NO[3],""^"-"," (kg N ",ha^-1, yr^-1,")")))
mtext("(H)",side=3,line=-1.1,adj=.95)
plot(meanGPP.T~reduced.nitrogen,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab=expression(paste(NH[4],""^"+"," (kg N ",ha^-1, yr^-1,")")))
mtext("(I)",side=3,line=-1.1,adj=.95)

####Second Figure

par(mar=c(5,4.5,1.5,1),cex=1.02)
boxplot(meanGPP.T~species.composition,data=data,las=1,ylab=expression(paste("WUE (GPP/T) g C ",kg^-1," ",H[2],O)),xlab="Dominant tree species")

###########

data.sd = data[data$sdGPP.T <30,]

par(mfrow=c(3,3))
plot(sdGPP.T~elevation,data=data.sd,las=1)
plot(sdGPP.T~slope.degree,data=data.sd,las=1)
plot(sdGPP.T~aspect.degree,data=data.sd,las=1)

plot(sdGPP.T~mean.annual.T,data=data.sd,las=1)
plot(sdGPP.T~mean.annual.RH,data=data.sd,las=1)
plot(sdGPP.T~mean.annual.P,data=data.sd,las=1)

plot(sdGPP.T~total.n.deposition,data=data.sd,las=1)
plot(sdGPP.T~oxidized.nitrogen,data=data.sd,las=1)
plot(sdGPP.T~reduced.nitrogen,data=data.sd,las=1)

tapply(data$meanGPP.T,data$species.composition,FUN=mean,na.rm=T)
tapply(data$meanGPP.T,data$species.composition,FUN=sd,na.rm=T)

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

        	annotation_scale(data=summer2020_df,plot_unit="m",location="br")+





















################Deprecated####################


setwd("/Users/goldsmit/Desktop/BernardoProject/Switzerland/WUEFiles/")

wue = do.call(rbind,lapply(list.files(path = "/Users/goldsmit/Desktop/BernardoProject/Switzerland/WUEFiles/"),read.csv))
write.csv(wue,"wue.csv")

setwd("/Users/goldsmit/Desktop/BernardoProject/Switzerland/ETFiles/")

et = do.call(rbind,lapply(list.files(path = "/Users/goldsmit/Desktop/BernardoProject/Switzerland/ETFiles/"),read.csv))
write.csv(et,"et.csv")
wue.et = merge(wue,et,by=wue$Date)

agg.wue = do.call(data.frame,aggregate(wue$ECO4WUE_001_Water_Use_Efficiency_WUEavg,list(wue$ID),FUN=function(x) c(mean=mean(x,na.rm=T),sd=sd(x,na.rm=T),n=length(x))))
names(agg.wue) = c("siteID","meanWUE","sdWUE","nWUE")

