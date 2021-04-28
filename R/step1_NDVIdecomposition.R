#############################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step 1: decomposition of NDVI time-series of random MODIS   
#forested pixels over Spain downloaded and pre-processed throught Google earth engine
#################################
#rm(list = ls())#empty environment
library(DBEST)
library(zoo)
library(xts)
library(graphics)
library(stats)
library(Hmisc)
library(stringr)
library(proj4)
library(sp)
library(gdata)
##########function to turn Null into NA
nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}
######################################
setwd("./Khoury&Coomes_GCB/data")#set working directory
######################################organize MODIS dat for analysis 
modispts1=read.csv("modisNDVI1.csv", header=T,stringsAsFactors=FALSE)
modispts2=read.csv("modisNDVI2.csv", header=T,stringsAsFactors=FALSE)
modispts3=read.csv("modisNDVI3.csv", header=T,stringsAsFactors=FALSE)
modispts=rbind(modispts1,modispts2,modispts3)
modispts<-modispts[2:9]#filter not needed columns
modispts=unique(modispts)#remove dubplicates
modispts$date<-as.Date(as.character(modispts$date), format='%Y%m%d')#format date
modispts$NDVI<-modispts$NDVI*0.0001 #scale NDVI values
coords<-unique(modispts[,6:7]) #get unique coordinates
coords$points<-c(1:5026)
##########create initial dataframe with pixel Id, NUmber of msing observation per time-series, latitude, longitude,
# slope at the location, domiminant species at the location,density of all species at location, elevation and location,
#and type of DBEST detected 4 biggest changes: gradual vs abrupt (abrupt assumed being caused by fires,clearcut,sattelite sensor problems)
Dbestdata<-c()
############## decompose series using DBEST
for (p in coords$points) { tryCatch({
     coordp<-coords[p,1:2]
     ts=subset(modispts,latitude==coordp$latitude & longitude==coordp$longitude)
     dom=unique(ts$dominantspp)
     slop<-unique(ts$slope)
     dens=unique(ts$density)
     elev=unique(ts$elevation)
     ind<-zoo(order.by=seq.Date(from = as.Date('2000-02-01'), to=as.Date('2017-12-01'), by='month'))
     tszoo<-zoo(ts$NDVI,ts$date)
     numb<-length(coredata(tszoo))
     tszoo<-merge(tszoo,ind)
     tszoo<-na.approx(tszoo, rule=2, maxgap=3)
     #detrendzoo<-lm(coredata(tszoo)~index(tszoo))$residuals
     tf=anyNA(tszoo) #make sure no nas left in datasets
     pts=ts(tszoo,start=c(2000,02,01), end=c(2017,12,01), frequency=12)
     dbestmod=paste("xdbestmod", p, sep="_")
     dbestmod=assign(dbestmod, DBEST(data=pts, data.type="cyclical",
                                     seasonality=12, algorithm="change detection",
                                     breakpoints.no = 4, first.level.shift=0.1,
                                     second.level.shift=0.2, duration=24,
                                     distance.threshold="default", alpha=0.05, plot="off"))
     dbestmod=nullToNA(dbestmod)
    
     Dbestdata= rbind(Dbestdata, data.frame(ID=paste("xdbestmod", p, sep="_"),misobs=tf,Long=unique(ts$longitude),Lat=unique(ts$latitude),
                                            slope=slop,Dominant=dom,Density=dens,elevation=elev,changeType=paste(dbestmod$ChangeType, collapse = ', '),stringsAsFactors = F))
    },error=function(e){cat("ERROR:",conditionMessage(e),"\n")})}

##transform dataframe into spatial dataframe and set crs.
coordinates(Dbestdata)<-c('Long','Lat')
wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
proj4string(Dbestdata)<-wgs

################################removing pixels with any abrupt change 
Dbestdata2<-Dbestdata[Dbestdata$ChangeType!='1, 0, 0, 1',]
Dbestdata2<-Dbestdata2[Dbestdata2$ChangeType!='1, 0, 0, 0',]
Dbestdata2<-Dbestdata2[Dbestdata2$ChangeType!='1, 0, 1, 0',]
Dbestdata2<-Dbestdata2[Dbestdata2$ChangeType!='0, 1, 0, 0',]
Dbestdata2<-Dbestdata2[Dbestdata2$ChangeType!='0, 0, 0, 1',]
Dbestdata2<-Dbestdata2[Dbestdata2$ChangeType!='0, 0, 0',]

#44 extreme changes
#647 missing data
#from 3855 to 3307 missing values
#from 3307 to 3271 extreme changes
############################# clear environment and save environemnt with dataframe and individual time-series.
rm(modispts,ts,tszoo,pts,coordp,coords,dbestmod,tf,modispts1,modispts2,modispts3,Dbestdata)
save.image('./Khoury&Coomes_GCB/processed/step1_NDVIdata.RData')
