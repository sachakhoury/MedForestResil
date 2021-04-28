#############################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step4: Decomposition of NDVI time-series running DBEST package to extract short-term resilience 
#metrics: sens=ndvi,lai loss; recov=ndvi,lai gain; resil=ndvi,lai gain/loss.
################################# libraries ########################################################################
library(DBEST)
library(zoo)
library(xts)
library(graphics)
library(stats)
library(Hmisc)
library(stringr)
library(raster)
###################################################################################################################
load("./Khoury&Coomes/processed/step3_timeseriesRegressionData.RData")
dates<-c(1:215)
DbestdataNDVI<-c()
############## decompose series using DBEST# #######################################################################
for (p in ndvimetrics@data$ID) { tryCatch({
    tszoo<-get(p)$Data
    pts=ts(tszoo,start=c(2000,02,01), end=c(2017,12,01), frequency=12)
    dbestmod=p
    dbestmod=assign(dbestmod, DBEST(data=pts, data.type="cyclical",
                                    seasonality=12, algorithm="change detection",breakpoints.no = 30,
                                    first.level.shift=0.1,
                                    second.level.shift=0.2, duration=24,
                                    distance.threshold=0.1, alpha=0.05, plot="fig1"))#Dbest decomposition
    scale<- ndvimetrics@data[ndvimetrics@data$ID==p,13]#makesure it's the correct column for spei-scale
    spi<-get(paste('spei',scale, sep=''))[,p]
    spi<-lowess(coredata(spi),f=0.1)$y
    # extracting first resilience and drought metrics 
    sens1<-min(dbestmod$Change[dbestmod$Change<0&dbestmod$Start<170])
    decrate1<-sens1/(dbestmod$Duration[dbestmod$Change==sens1])
    start1<-dates[dbestmod$Start[dbestmod$Change==sens1]+(dbestmod$Duration[dbestmod$Change==sens1]/2)]
    minspi1<-min(spi[dbestmod$Start[dbestmod$Change==sens1]:dbestmod$End[dbestmod$Change==sens1]],na.rm=T)
    
    num=min(dbestmod$Start[dbestmod$Start>=dbestmod$End[dbestmod$Change==sens1]])
    recov1<-dbestmod$Change[dbestmod$Start==num]
    elast1<-recov1/(dbestmod$End[dbestmod$Change==recov1]-dbestmod$End[dbestmod$Change==sens1])
    maxspi1<-max(spi[dbestmod$Start[dbestmod$Change==recov1]:dbestmod$End[dbestmod$Change==recov1]],na.rm=T)
    if(recov1<0){
      recov1=dbestmod$Fit[num]-dbestmod$Fit[dbestmod$End[dbestmod$Change==sens1]]
      elast1=recov1/(num-dbestmod$End[dbestmod$Change==sens1])
      maxspi1<-max(spi[dbestmod$End[dbestmod$Change==sens1]:num],na.rm=T)
      }
    # extracting second resilience and drought metrics
    sens2<-min(dbestmod$Change[dbestmod$Change<0 & dbestmod$Change>sens1&dbestmod$Start<170])
    decrate2<-sens2/(dbestmod$Duration[dbestmod$Change==sens2])
    start2<-dates[dbestmod$Start[dbestmod$Change==sens2]+(dbestmod$Duration[dbestmod$Change==sens2]/2)]
    minspi2<-min(spi[dbestmod$Start[dbestmod$Change==sens2]:dbestmod$End[dbestmod$Change==sens2]],na.rm=T)
   
    num=min(dbestmod$Start[dbestmod$Start>=dbestmod$End[dbestmod$Change==sens2]])
    recov2<-dbestmod$Change[dbestmod$Start==num]
    elast2<-recov2/(dbestmod$End[dbestmod$Change==recov2]-dbestmod$End[dbestmod$Change==sens2])
    maxspi2<-max(spi[dbestmod$Start[dbestmod$Change==recov2]:dbestmod$End[dbestmod$Change==recov2]],na.rm=T)
    if(recov2<0){
      recov2=dbestmod$Fit[num]-dbestmod$Fit[dbestmod$End[dbestmod$Change==sens2]]
      elast2=recov2/(num-dbestmod$End[dbestmod$Change==sens2])
      maxspi2<-max(spi[dbestmod$End[dbestmod$Change==sens2]:num],na.rm=T) 
    }
    
  sens<-mean(sens1,sens2,na.rm=T)*(-1)
  recov<-mean(recov1,recov2,na.rm=T) 
  elast<-mean(elast1,elast2,na.rm=T)
  decrate<-mean(decrate1,decrate2,na.rm=T)*(-1)
  minspei<-mean(minspi1,minspi2,na.rm=T)
  maxspei<-mean(maxspi1,maxspi2,na.rm=T)
  resil<-recov/sens
  t<-abs((dbestmod$Start[dbestmod$Change==sens1])-(dbestmod$Start[dbestmod$Change==sens2]))
  ndvimean=mean(dbestmod$Trend)
  
  DbestdataNDVI= rbind(DbestdataNDVI, data.frame(ID=p,ndvimean=ndvimean,
                     Duration=t,senstv=sens,derate=decrate,recov=recov,elast=elast,resil=resil,minspei=minspei,maxspei=maxspei,drought1=start1,drought2=start2,stringsAsFactors = F))
},error=function(e){cat("ERROR:",conditionMessage(e),"/n")})}

DbestdataNDVI<-DbestdataNDVI[DbestdataNDVI$recov>0,]
DbestdataNDVI$middrought<-(DbestdataNDVI$drought1+DbestdataNDVI$drought2)/2
#ndvi metrics table with: ID-pixel id, ndvimean-average NDVI over pixel, Duration-duration of the drought,
#senstv-NDVI loss, decrate-NDVI loss rate, recov-NDVI gain after drought, elast-NDVI gain rate
#resil-NDVI gain over loss, minspei-minimum SPEI during drought, maxspei-maximum SPEI after during recovery period,
#drought1-month of first strongest drought drought2-month of second strongest drought

#############################################Merging dataset and creating resilience metrics######################
ndvimetrics<-merge(ndvimetrics,DbestdataNDVI,by=c('ID'),all.x=FALSE)
##################################################################

DbestdataLAI<-c()
############## decompose series using DBEST
for (p in Dbestdata2@data$ID) { tryCatch({
    tszoo<-get(p)$Data
    pts=ts(tszoo,start=c(2000,02,01), end=c(2017,12,01), frequency=12)
    dbestmod=p
    dbestmod=assign(dbestmod, DBEST(data=pts, data.type="cyclical",
                                    seasonality=12, algorithm="change detection",breakpoints.no = 30,
                                    first.level.shift=0.1,
                                    second.level.shift=0.2, duration=24,
                                    distance.threshold=0.1, alpha=0.05, plot="fig1"))
    scale<- laimetrics@data[laimetrics@data$ID==p,13]#makesure it's the correct column for spei-scale
    spi<-get(paste('spei',scale, sep=''))[,p]
    spi<-lowess(coredata(spi),f=0.1)$y
    # extracting first resilience and drought metrics 
    sens1<-min(dbestmod$Change[dbestmod$Change<0&dbestmod$Start<170])
    decrate1<-sens1/(dbestmod$Duration[dbestmod$Change==sens1])
    start1<-dates[dbestmod$Start[dbestmod$Change==sens1]+(dbestmod$Duration[dbestmod$Change==sens1]/2)]
    minspi1<-min(spi[dbestmod$Start[dbestmod$Change==sens1]:dbestmod$End[dbestmod$Change==sens1]],na.rm=T)
  
    num=min(dbestmod$Start[dbestmod$Start>=dbestmod$End[dbestmod$Change==sens1]])
    recov1<-dbestmod$Change[dbestmod$Start==num]
    elast1<-recov1/(dbestmod$End[dbestmod$Change==recov1]-dbestmod$End[dbestmod$Change==sens1])
    maxspi1<-max(spi[dbestmod$Start[dbestmod$Change==recov1]:dbestmod$End[dbestmod$Change==recov1]],na.rm=T)
    if(recov1<0){
      recov1=dbestmod$Fit[num]-dbestmod$Fit[dbestmod$End[dbestmod$Change==sens1]]
      elast1=recov1/(num-dbestmod$End[dbestmod$Change==sens1])
      maxspi1<-max(spi[dbestmod$End[dbestmod$Change==sens1]:num],na.rm=T)
    }
    # extracting second resilience and drought metrics
    sens2<-min(dbestmod$Change[dbestmod$Change<0 & dbestmod$Change>sens1&dbestmod$Start<170])
    decrate2<-sens2/(dbestmod$Duration[dbestmod$Change==sens2])
    start2<-dates[dbestmod$Start[dbestmod$Change==sens2]+(dbestmod$Duration[dbestmod$Change==sens2]/2)]
    minspi2<-min(spi[dbestmod$Start[dbestmod$Change==sens2]:dbestmod$End[dbestmod$Change==sens2]],na.rm=T)
    
    
    num=min(dbestmod$Start[dbestmod$Start>=dbestmod$End[dbestmod$Change==sens2]])
    recov2<-dbestmod$Change[dbestmod$Start==num]
    elast2<-recov2/(dbestmod$End[dbestmod$Change==recov2]-dbestmod$End[dbestmod$Change==sens2])
    maxspi2<-max(spi[dbestmod$Start[dbestmod$Change==recov2]:dbestmod$End[dbestmod$Change==recov2]],na.rm=T)
    if(recov2<0){
      recov2=dbestmod$Fit[num]-dbestmod$Fit[dbestmod$End[dbestmod$Change==sens2]]
      elast2=recov2/(num-dbestmod$End[dbestmod$Change==sens2])
      maxspi2<-max(spi[dbestmod$End[dbestmod$Change==sens2]:num],na.rm=T) 
    }
  sens<-mean(sens1,sens2,na.rm=T)*(-1)
  recov<-mean(recov1,recov2,na.rm=T) 
  elast<-mean(elast1,elast2,na.rm=T)
  decrate<-mean(decrate1,decrate2,na.rm=T)*(-1)
  minspei<-mean(minspi1,minspi2,na.rm=T)
  maxspei<-mean(maxspi1,maxspi2,na.rm=T)
  resil<-recov/sens
  t<-abs((dbestmod$Start[dbestmod$Change==sens1])-(dbestmod$Start[dbestmod$Change==sens2]))
  LAImean=mean(dbestmod$Trend)
  
  DbestdataLAI= rbind(DbestdataLAI, data.frame(ID=p,LAImean=LAImean, Duration=t,senstv=sens,decrate=decrate,recov=recov,elast=elast,resil=resil,
                 minspei=minspei,maxspei=maxspei,drought1=start1,dought1=start2,stringsAsFactors = F))
},error=function(e){cat("ERROR:",conditionMessage(e),"/n")})}

DbestdataLAI<-DbestdataLAI[DbestdataLAI$recov>0,]
DbestdataLAI$middrought<-(DbestdataLAI$drought1+DbestdataLAI$drought2)/2
#lai metrics table with: ID-pixel id, LAImean-average Lai over pixel, Duration-duration of the drought,
#senstv-lai loss, decrate-lai loss rate, recov-lai gain after drought, elast-lai gain rate
#resil-lai gain over loss, minspei-minimum SPEI during drought, maxspei-maximum SPEI after during recovery period,
#drought1-month of first strongest drought drought2-month of second strongest drought
###################Merging dataset and creating resilience metrics######################
laimetrics<-merge(laimetrics,DbestdataLAI,by=c('ID'),all.x=FALSE)

##clean-up
keep(laimetrics,ndvimetrics,wgs,sure=T)
#######################################################################################
save.image('./Khoury&Coomes/processed/step4_segmentation&continuousData.RData')
#######################################################################################



