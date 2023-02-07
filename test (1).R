library(raster)
library(reshape2)

# get SD

#create directory to draw all files from
wue_dir <- dir('~/Documents/Chapman Research/ECOSTRESS Project/Switzerland/', full.names = T) #I assume these are .tif files
#wue_dir <-wue_dir[-c(1,32,33,34)] #If you want to remove certain files (or just choose a couple to test code)
wue_stack <-raster::stack(wue_dir) #stack all files

#make raster stack into a dataframe
wue_stack_df<-rasterToPoints(wue_stack)
wue_stack_df <- as.data.frame(wue_stack_df)

#this should produce a large 'wide' format dataframe where each 'stack' is a column
# if this is the case, then try this

wue_stack_df  <- melt(wue_stack_df, 
                      id.vars = c("x", "y"),
                      variable.name = "month") #melt to long format

#clean up month column
wue_stack_df$month<-gsub('~/Documents/Chapman Research/ECOSTRESS Project/Switzerland/Switz_WUE_',
                        '', npp_stack_df$month)

# an so on until you have a dataframe with columns: x,y,month,wue (or whatever you name them)

#get sd
wue_sd<-aggregate(wue~x+y,sd,data=wue_stack_df)
wue_sd_raster<-rasterFromXYZ(wue_sd)
#some code for color gradients:

# color gradient/plot
et= c("red", "orange", "yellow",'green','cyan3','purple')
#et <- c('red', 'white','blue')
bks_et<- quantile(wue_stack_df$wue, probs=seq(0.0, 1, by=0.10), na.rm = TRUE)
bkcols.et <- colorRampPalette(et)(length(bks_et)-1)
#proj4string(x)<-CRS(aea.proj)
r.range <- round(c(minValue(wue_stack_df), maxValue(wue_stack_df)))

plot(wue_sd_raster,breaks = bks_et,axes=F,box=F,legend = F, col = bkcols.et,
     legend.width=0.80,legend.shrink=1,main='Time Series',cex.main=1,
     axis.args=list(at=seq(r.range[1], r.range[2], 10),
                    labels=seq(r.range[1], r.range[2], 10),
                    cex.axis=0.75),
     legend.args=list(expression(paste('')), side=4, font=2,adj=0.5, line=2.5, cex=0.9))

# maybe turn the shapfile into raster first?
#https://stackoverflow.com/questions/35096133/converting-shapefile-to-raster

masked = crop(orginal, forests.raster) 
npp_50_threshold_raster <-mask(masked, forests.raster)

