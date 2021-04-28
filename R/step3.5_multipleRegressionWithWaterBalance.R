################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step3.5: Determining average water balance (WB) of Spain, computed from GLDAS 2.1 environmental data (WB=precipitation-potential evapotranpiration). 
#Repetition of step 3 using WB instead of SPEI (sensitivity analysis) so coefficient alphaC, betaC and epsilonC all relate to WB and not SPEI here.
################################ libraries ############################################################################
library(tseries)
library(rgdal)
library(DBEST)
library(zoo)
library(SPEI)
library(lubridate)
library(Hmisc)
library(gdata)
################################################## creating time-series of water balance ########################################################################
####load precipitation,temprature,windspeed,mean dialy incoming solar radiation to compute WB timeseries and exctract aveage for pixels
setwd('./Khoury&Coomes_GCB/data/MODIS/water_balance')
precipitation<-read.csv('precip_table.csv',stringsAsFactors = F)
precipitation<-precipitation[,c(3:6,8)]
precipitation$precip<-precipitation$precip*86400 #convert precipitation from precipitation rate	in kg/m^2/s to mm/day
precipitation$date<-as.Date(as.character(precipitation$date), format='%Y%m%d')
precipitation$days<-monthDays(precipitation$date)
precipitation$precipmon<-precipitation$precip*precipitation$days#mm/month
precipitation$precip<-NULL
precipitation$days<-NULL
tmin<-read.csv('tmin_table.csv',stringsAsFactors = F)
tmin$date<-as.Date(as.character(tmin$date), format='%Y%m%d')
tmin<-tmin[,c(3:6,8)]
tmin$tmin<-tmin$tmin-273.15 #convert temperature from kelvin to celcius
tmax<-read.csv('tmax_table.csv',stringsAsFactors = F)
tmax$date<-as.Date(as.character(tmax$date), format='%Y%m%d')
tmax<-tmax[,c(3:6,8)]
tmax$tmax<-tmax$tmax-273.15 #convert temperature from kelvin to celcius
tmean<-read.csv('tmean_table.csv',stringsAsFactors = F)
tmean$date<-as.Date(as.character(tmean$date), format='%Y%m%d')
tmean<-tmean[,c(3:6,8)]
tmean$tmean<-tmean$tmean-273.15 #convert temperature from kelvin to celcius
wind<-read.csv('wind_table.csv',stringsAsFactors = F)
wind$date<-as.Date(as.character(wind$date), format='%Y%m%d')
wind<-wind[,c(3:6,8)]# wind m/s
rs<-read.csv('rs_table.csv',stringsAsFactors = F)
rs$date<-as.Date(as.character(rs$date), format='%Y%m%d')
rs<-rs[,c(3:8)]
rs$Rs<-rs$Rs*0.0864#convert from W/m^2 to MJ m-2 d-1

#####combine tables
gldas<-merge(precipitation,tmax, by=c('ID','date','latitude','longitude'))
gldas<-merge(gldas,tmin, by=c('ID','date','latitude','longitude'))
gldas<-merge(gldas,wind, by=c('ID','date','latitude','longitude'))
gldas<-merge(gldas,rs, by=c('ID','date','latitude','longitude'))
rm(rs,precipitation,tmean,tmax,tmin,wind)#delete raw data

######Create the time-series of WB and populate it
waterbalance1<-data.frame(unique(gldas$date))
colnames(waterbalance1)<-c('Date')
####in waterbalance1 the column name indicate the pixel ID and the row represent the time-series for these pixels IDs at scale 1 (i.e. water accumulation not acounted for) 
for (l in gldas$ID){tryCatch({ 
  ts=subset(gldas,ID==l)
  ts=unique(ts)
  pet<-penman(Tmin=ts$tmin,Tmax=ts$tmax,U2=ts$wind,lat=lonlat$latitude, Rs=ts$Rs,z=unique(ts$elevation), na.rm=T)#calculate penman montheith PET for each plot
  ts$pet<-pet[,1]
  ts$wb<-ts$precip-ts$pet#calculate WB
  tszoo<-zoo(ts[4:15],ts$date)
  ind<-zoo(order.by=seq.Date(from = as.Date('2000-02-01'), to=as.Date('2017-12-01'), by='month'))
  tszoo<-merge(tszoo,ind)
  tszoo<-na.approx(tszoo, rule=2, maxgap=4)
  tf=anyNA(tszoo)
  waterbalance1<-cbind(waterbalance1,coredata(tszoo$wb))
  colnames(waterbalance1)[length(waterbalance1)]<-l
},error=function(e){cat("ERROR:",conditionMessage(e),"\n")})}

########################### Get average yearly water balance at each pixel to be used in resilience analysis ############################# 
WB<-rowsum(waterbalance1[,2:3900],rep(1,215))/18
WB<-t(WB)
env<-data.frame(rownames(WB), WB[,1])
colnames(env)<-c('ID','WB')

########################################################## working with waterbalance1 ##################################################################################
years2<-seq(1,215, by=1)
resmultiregtrend<-c()
for (y in Dbestdata2$ID){
  tryCatch({
    x<-get(y)
    xcoord<-Dbestdata2@coords[Dbestdata2$ID==y]
    green<-lm(x$Trend~index(x$Trend))
    s1<-waterbalance1[,y]
    s1<-lowess(coredata(s1),f=0.1)$y
    s2<-lm(s1~years2)
    betaC<-coef(summary(s2))[1,1]#intercept
    alphaC<-coef(summary(s2))[2,1]#trend in wb
    epsilonC<-s2$residuals#detrended wb 
    #s<-lag(s,k=ccflag,na.pad = T)#could be used to evaluate lag of greenness response, however does not include drought accumulation, result obtain were not very usefull  
    ndvilm0<-lm(x$Trend~years2)#model ndvi and time
    alphaG<-coef(summary(ndvilm0))[2,1]#observed greening
    betaG<-coef(summary(ndvilm0))[1,1]#observed greening intercept
    ndvilm1<-lm(x$Trend~years2+s1) #model with spei and time
    alphaP<-coef(summary(ndvilm1))[2,1]#potential greening if not affected by drought
    gammaP<-coef(summary(ndvilm1))[3,1]#short-term responsiveness
    betaP<-coef(summary(ndvilm1))[1,1]#intercept
    epsilonP<-ndvilm1$residuals
    #alphaG<-alphaP+gammaG*alphaC  #actual growth
    teffect<-alphaG-alphaP#long term responsiveness
    pteffect<-(alphaG-alphaP)/sqrt((coef(summary(ndvilm1))[2,2])^2+(coef(summary(ndvilm0))[2,2])^2)
    minwb<-min(s1)
    #summerspei<-mean(spei12[jun,y])
    resmultiregtrend<-rbind(resmultiregtrend,data.frame(ID=y,wbR2=summary(ndvilm1)$r.squared,wbpval= lmp(ndvilm1), 
                                               wbtime=alphaG, wbtime_pval=coef(summary(ndvilm0))[2,4],wbestimgrowth=alphaP,
                                               wbteffect=teffect,wbteffect_pval=pteffect,wb=gammaG,wb_pval=coef(summary(ndvilm1))[3,4],
                                              wbtrend=alphaC,wbtrend_pval=coef(summary(s2))[2,4],minwb=minspei,long=xcoord[1],lat=xcoord[2]))
  },error=function(e){cat("ERROR:",conditionMessage(e),"/n")})}
####### see Suplmentary Information Table S1 for reference to models and methods for information about equations
#wbR2 adjusted r2 from model M1 with WB instead of SPEI 
#wbpval pvalue of previous model
#wbtime coefficient of greenness change as determined from M4 with WB instead of SPEI
#wbtime_pval pvalue of coefficient of greenness change 
#wbestimgrowth potential growth in absence of climate change as determined with WB instead of SPEI
#wbteffect climate change effect on greeness as determined from M4 using WB instead of SPEI
#wbteffect_pval pvalue of climate change effect on greenness
#wb covariance between WB and NDVI instead of SPEI and NDVI as determined from model M4
#wb_pval pvalue of WB coefficient in M4 using WB instead of SPEI
#wbtrend trend of WB taken from M1 when using WB instead of SPEI
#wbtrend_pval pvalue of M1 when using WB instead of SPEI
#minwb determining minimum WB from time-series of WB
#long longitude of evaluate pixel
#lat latitude of evaluated pixel
coordinates(resmultiregtrend)<-c('long','lat')#transform into sptial points
proj4string(resmultiregtrend)<-Dbestdata2@proj4string

#create metrics table to be compared with metrics detemined with SPEI in figure 4 and suplementary information
ndvimetrics_wb<-resmultiregtrend
ndvimetrics_wb<-merge(ndvimetrics_wb,env,by='ID')

################################################################################################
keep(ndvimetrics_wb,sure=T)
save.image("./Khoury&Coomes/processed/step3.5_regressionDataWithWB.RData")

