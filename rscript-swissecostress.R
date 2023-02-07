setwd("/Users/goldsmit/Dropbox/Switzerland_IAPDatabase/Proposals_Reports/InterimProgress-2020/ECOSTRESSSwitzerland")

##------------------------------------------------##
##ET Data as a Predictor of SOI
##------------------------------------------------##

ecos = read.csv("2019-Sampling-Sites-Switzerland-ECO3ETPTJPL-001-results.csv",header=TRUE)
str(ecos)

soi.2019 = read.csv("2015.2019merge.csv",header=TRUE)
str(soi.2019)

##Aggregate the data:

ecos$Trans = ecos$ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily*(ecos$ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETcanopy/100)

ecos.agg = aggregate(ecos$ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily,list(ecos$ID,ecos$Year,ecos$Month),FUN=median)
names(ecos.agg) = c("siteID","Year","Month","ETDaily")
str(ecos.agg)

##Subset to the summer: 
ecos.summer = ecos.agg[which(ecos.agg$Month == c(6:8)),]
str(ecos.summer)

###Look at the data: 
boxplot(ETDaily~Year,data=ecos.summer,las=1)
boxplot(ETDaily~siteID,data=ecos.summer,las=1)
boxplot(ETDaily~Month,data=ecos.summer,las=1)

##Merge the two data streams: 
	##I tried to plot the all summer all years vs. soi - no pattern
	
ecos.site = aggregate(ecos.summer$ETDaily,list(ecos.summer$siteID),FUN=median)
names(ecos.site) = c("siteID","ETDaily")

ecos.soi = merge(soi.2019,wue.agg)
str(ecos.soi)

#The plot: 

plot(SOI.2019~WUE,data=ecos.soi,las=1) ###No pattern across all summer all years
abline(lm(SOI.2019~WUE,data=ecos.soi))

##------------------------------------------------##
##WUE 
##------------------------------------------------##

wue = read.csv("2019-Sampling-Sites-Switzerland-ECO4WUE-001-results.csv",header=TRUE)
str(wue)

wue.agg = aggregate(wue$ECO4WUE_001_Water_Use_Efficiency_WUEavg,list(wue$ID),FUN=mean)
names(wue.agg) = c("siteID","WUE")
str(wue.agg)

boxplot(wue.agg$WUE~wue.agg$siteID)

