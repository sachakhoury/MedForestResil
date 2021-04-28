#################################### 
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step3: Multiple linear regression with NDVI, SPEI and TIME computingshort-term (CGW) and
#long-term climate change effect (CCI) in the continuous analysis of forest canopy resilience.
###############################################
library(tseries)
library(rgdal)
library(DBEST)
library(zoo)
library(gdata)
####################functions###################################################
####get the highest correlation from the  cross correlation function between NDVI and the different SPEI scales
Find_Max_CCF<- function(a,b){
  d <- ccf(coredata(a), b, plot = F)
  cor = d$acf[,,1]
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  res=subset(res, lag ==0)
  res_max = res[which.max(res$cor),]
  return(res_max)
} 

##### get overall p-value from lm object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
######################################determining relationship between LAI and NDVI for year 2013 before the biggest drought happens#######################
modislai=read.csv("./Khoury&Coomes_GCB/data/MODIS/modislai2013.csv", header=T,stringsAsFactors=FALSE)
modislai<-modislai[2:6]
modislai$Lai<-modislai$Lai/10
modislai$date<-as.Date(as.character(modislai$date), format='%Y%m%d')
modislai1<-modislai[modislai$date=='2013-08-01',]
#modislaimean<-aggregate(modislai,FUN=mean, by=list(modislai$ID))
#modislaimean$ID<-NULL
#colnames(modislaimean)[1]<-'ID'
modisNDVI=read.csv("./Khoury&Coomes_GCB/data/MODIS/modisndvi2013.csv", header=T,stringsAsFactors=FALSE)
modisNDVI<-modisNDVI[2:6]
modisNDVI$NDVI<-modisNDVI$NDVI*0.0001
modisNDVI$date<-as.Date(as.character(modisNDVI$date), format='%Y%m%d')
modisNDVI1<-modisNDVI[modisNDVI$date=='2013-08-01',]
#modisNDVImean<-aggregate(modisNDVI,FUN=mean, by=list(modisNDVI$ID))
#modisNDVImean$ID<-NULL
#colnames(modisNDVImean)[1]<-'ID'
df<-merge(modislai1,modisNDVI1,by='ID',all.x=F)
summary(lmlai<-lm(NDVI~log(Lai),df))
summary(nlmlai<-nls(NDVI~a*(1-exp(-b*(Lai-c))),start=list(a=1,b=0.6,c=0.1),df))
plot(lmlai)

#######################################detemining drought accumulation period#####################################################################
###################################and getting CGW and CCI from equation 3 in paper##########################
years2<-seq(1,215, by=1)# months used in the study
resmultiregtrend<-c()
for (y in Dbestdata2$ID){
  tryCatch({
  x<-get(y)$Trend
  #laitrans<-function(x){0.02-((log((-x/0.85)+1))/0.78)}# add commented lines to compute for LAI insteas of NDVI
  #x<-laitrans(x)
  xcoord<-Dbestdata2@coords[Dbestdata2$ID==y]
  green<-lm(x~index(x))
  ccftable<-c('cor','lag')
  
  #run loop to determine best SPEI scale creating a table with correlation number and spei scale evaluated
    for (l in c(1:48)){
    s<-get(paste('spei',l,sep=''))[,y]
    s<-lm(coredata(s)~index(s))$residuals
    x2<-lm(x~index(x))$residuals
    ccftable<-rbind(ccftable,Find_Max_CCF(ts(lowess(coredata(s),f=0.1)$y),ts(x2)))
    }#end of loop
  
  ccftable<-data.frame(ccftable[2:49,],row.names =c(1:48) ,stringsAsFactors = F)
  colnames(ccftable)<-c('cor','lag')
  ccftable$lag<-as.numeric(ccftable$lag)
  ccfresult<-ccftable[ccftable$cor==max(as.numeric(ccftable$cor)),]
  ccflag<-ccfresult[,2]
  ccfcor<-ccfresult[,1]
  ccfspei<-rownames(ccfresult)
  s1<-get(paste('spei',ccfspei,sep=''))[,y] #best scale for pixel
  s1<-lowess(coredata(s1),f=0.1)$y
  s2<-lm(s1~years2)# SPEI model with time #fitting equation 1 in paper SPEIx(t)=alphaC.t+betaC+epsilonC(t)
  betaC<-coef(summary(s2))[1,1]#intercept
  alphaC<-coef(summary(s2))[2,1]#trend coefficient in spei
  speifit<-summary(s2)$adj.r.squared
  speiaic<-AIC(s2)
  epsilonC<-s2$residuals#detrended spei 
  #s<-lag(s,k=ccflag,na.pad = T)
  ndvilm0<-lm(x~years2)#ndvi linear model with time #fitting equation 2 in paper NDVI(t)=alphaG.t+betaG+epsilonG(t)
  ndvifit<-summary(ndvilm0)$adj.r.squared
  ndviaic<-AIC(ndvilm0)
  ndvimeantrue<-mean(x)
  alphaG<-coef(summary(ndvilm0))[2,1]#observed greening coefficient
  betaG<-coef(summary(ndvilm0))[1,1]#observed greening intercept
  ndvilm1<-lm(x~years2+s1) # NDVI linear model with spei and time #fitting equation 3 in paper NDVI(t)=alphaP.t+betaP+epsilonP(t)+gammaG.SPEIx(t)
  bothfit<-summary(ndvilm1)$adj.r.squared
  bothaic<-AIC(ndvilm1)
  alphaP<-coef(summary(ndvilm1))[2,1]#potential greening if not affected by drought (multiplies t)
  gammaG<-coef(summary(ndvilm1))[3,1]# CCI short-term covariance between ndvi and spei (multplies spei)
  betaP<-coef(summary(ndvilm1))[1,1]#intercept
  epsilonP<-ndvilm1$residuals
  #alphaG=alphaP+gammaG*alphaC  #actual growth
  teffect<-alphaG-alphaP#long term responsiveness
  pteffect<-(alphaG-alphaP)/sqrt((coef(summary(ndvilm1))[2,2])^2+(coef(summary(ndvilm0))[2,2])^2)#accept null Ho if -1.96<z<1.96
  trendspei12<-lm(spei12[,y]~years2)$coefficients[2]
  trendspei1<-lm(lowess(coredata(spei1[,y]),f=0.1)$y~years2)$coefficients[2]
  time2<-(214/(1+betaG/alphaG))*100#in percent term
  estimgrowth2<-(214/(1+betaP/alphaP))*100#in percent term
  teffect2<-time2-estimgrowth2#long term responsiveness in %
  resmultiregtrend<-rbind(resmultiregtrend,data.frame(ID=y,cor=as.numeric(ccfcor),R2=summary(ndvilm1)$r.squared, pval=lmp(ndvilm1), 
                                             time=alphaG,time1=time2,time_p=coef(summary(ndvilm0))[2,4],estimgrowth=alphaP,estimgrowth2=estimgrowth2,
                                             cci=teffect,cci_p=pteffect,cci2=teffect2,spei_scale=ccfspei,cgw=gammaP,cgw_p=coef(summary(ndvilm1))[3,4],
                                             trendspei=alphaC,trendspei_p=coef(summary(s2))[2,4],trendspei_x12=trendspei12,trendspei_x1=trendspei1,spei_r2=speifit,spei_aic=speiaic,ndvi_r2=ndvifit,ndvi_aic=ndviaic,ndvimeantrue=ndvimeantrue,time_spei_r2=bothfit,time_spei_aic=bothaic,long=xcoord[1],lat=xcoord[2]))
},error=function(e){cat("ERROR:",conditionMessage(e),"/n")})}
#### see Suplmentary Information Table S1 for reference to models and methods for information about equations
#cor : correlation between greenness and SPEI
#R2: r2 of the linear regression between greenness and SPEI
#pval: p value of the previous linear regression
#time: coefficient of change in greenness as result of time in model M4 (also coef from M2)
#time2: change in greeness as a result of time in model M4 expressed in % from starting value
#time_p: p value of coefficient time in M4
#estimgrowth: potential change in greenness if there was not drought effect estimated from model M4
#estimgrowth2: potential change in greenness expressed in % from starting greeenness 
#cci: climate change influence on greenness as determined from model M4 or time-estimegrowth
#cci_p: pvalue of climate change influence on greenness 
#cci2: climate change influence expressed in % change in greenness
#spei_scale: SPEI scale determined which confer highest correlation between greeness and SPEI
#cgw: covariance betweeen greeness and spei as determined from cofficient of SPEI in model M4
#cgw_p: pvalue of coefficient SPEI in M4
#trendspei: linear trend of SPEI as determined from regression with time (model M1)
#trendspei_p: pvalue of M1
#trendspei_x12: trend of spei evaluted at scale 12
#trendspei_x1: trend of spei evaluated at scale 1
#ndvi_r2: adjusted r2 of model M2
#ndvi_aic: AIC of model M2
#spei_r2: adjusted r2 of model M1
#spei_aic: AIC of model M1
#time_spei_r2: adjusted r2 of model M4
#time_spei_aic: AIC of model M4
#ndvimeantrue: average ndvi
#long: longitude of pixel sampled
#lat: latitude of pixel sampled
coordinates(resmultiregtrend)<-c('long','lat')#transform into sptial points
proj4string(resmultiregtrend)<-Dbestdata2@proj4string
######################################################
ndvimetrics<-resmultiregtrend
laimetrics<-resmultiregtrend#after run with LAI
##################################koppen classes
rm(resmultiregtrend,Find_Max_CCF,lmp)
save.image('./Khoury&Coomes/processed/step3_timeseriesRegressionData.RData')
