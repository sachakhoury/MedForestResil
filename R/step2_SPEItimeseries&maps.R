##########################################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step2: Extraction of SPEI time-series from SPEI-scale dataset in .nc format.
#and mapping SPEI change over spain and NDVI time-series
##########################################
library(zoo)
library(ncdf4)
library(raster)
library(proj4)
library(SPEI)
library(Hmisc)
library(gdata)
############################################exploring drought index first
#rm(list = ls())#empty environment
################function to smooth timeseries of raster values
impute.lowess <- function(x, y = NULL, f =0.1, iter = 3, delta = 0.01 * diff(range(x))) {
  p<-lowess(x,f=0.1)
  y<-p$y
  return(y)
}
#################################
#creating raster layer of spei change 2000-2017 and getting relevant dates
n1 <-'./Khoury&Coomes_GCB/data/Spanish SPEI_1.1km_48scales/spei1.nc'
d1 <- nc_open(n1)
rbrick1 <- brick(n1,varname="value",level=1)
datenames<-rbrick1@data@names
datenames<-datenames[1877:2736]# dates 2000 to 2018
datenames3<-substr(datenames, 1, 8)
datenames3<-factor(datenames3)
rlayer<-subset(rbrick1,datenames)#subset dataset to get relevant dates
rlayer3<-stackApply(rlayer, indices=as.numeric(datenames3), fun=mean , na.rm=TRUE)#make sure values are numeric
##smooth SPEI values to determine linear trend
e22 <- rlayer3#smoothed spei maps
e22[] <- NA
for( rl in 1:nrow(rlayer3) ) { 
  v <- getValues(rlayer3, rl, 1)
  e22[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=impute.lowess)))
}
##############functions to apply linear model across time (216 months int total)
#and get coefficients a and b (y=ax+b with x being time 1:18 and y being SPEI)
fun2 <- function(x) {                      
   m <- NA 
   try( m <- lm(x~ time)$coefficients[2] ,silent=T) #get the time coefficient
   m } y
fun3 <- function(x) {
  m <- NA
  try( m <-  summary(lm(lowess(coredata(x),f=0.1)$y ~ time))$coefficients[2,4] ,silent=T) #get the intercept
  m }
##############################
e22x <- calc(e22, fun2)#get lm coefficient
e33x <- calc(e22, fun3)#get lm coefficient
############################## compute change over 18 years (216 month) and project
e22xx<-e22x*216#get SPEI change over 18 years
proj4string(e22xx)<-CRS('+proj=utm +zone=30 +ellps=intl +towgs84=-131,-100.3,-163.4,-1.244,-0.02,-1.144,9.39 +units=m +no_defs')
e22xx<-projectRaster(e22xx,crs=wgs)
e22xx<-setMinMax(e22xx)
##get minimum and maximum change over 18 years
max(values(e22xx),na.rm=T)
min(values(e22xx),na.rm=T)

##get average and standard deviation of spei over 18 years to plot average spei
spei1ts<-cellStats(rlayer3, stat='mean', na.rm=TRUE)
spei1sd<-cellStats(rlayer3, stat='sd', na.rm=TRUE)
zoospei1spainsd<-zoo(spei1sd,seq(as.Date('2000/02/01'),as.Date('2017/12/01'),by='month'))
zoospei1spain<-zoo(spei1ts,seq(as.Date('2000/02/01'),as.Date('2017/12/01'),by='month'))
zoospei1spainsd<-lowess(coredata(zoospei1spainsd),f=0.1)$y
zoospei1spain<-lowess(coredata(zoospei1spain),f=0.1)$y
summary(lm(zoospei1spain~years2))

###############################get average water balance (wb) for spain 
#to be compared with SPEI map
Temp1<-raster('./Khoury&Coomes_GCB/data/Tair/Tair0.tif')#spain average temperature over 18y period
plot(Temp1)
Temp1<-raster::mask(Temp1,spain)#clip to spain
latt <- init(Temp1, 'y')#get latitudes
temps<-preStack(pattern = "*.tif$",path='./Khoury&Coomes_GCB/data/Tair')#get monthly values in deg C
temps<-stack(temps)
Precip<-preStack(pattern = "*.tif$",path='./Khoury&Coomes_GCB/data/Precip')#get monthly precipitation values in kg/m^2/s
Precip<-stack(Precip)
#
tempdates=seq(as.Date("2013-01-01"), as.Date("2013-12-31"), by="month")
temps<- setZ(temps,tempdates)# add dates to teprature data
names(temps) <- as.yearmon(getZ(temps))# turn dates into month year values
##################apply thornthwaite equation from SPEI package to compute potential evapotranspitation (pet)
#function to apply thornthwaite ET function to determin pet for raster cells
th <- function(Tave, lat) {
  as.vector(SPEI::thornthwaite(as.vector(Tave), lat,na.rm = T))
} 

pet <- brick(temps, values=FALSE) #stack pet rasters
ncell<-c(1:69345)[is.na(values(pet$layer.1))==T]#remove empty cells
tcell<-c(1:69345)[is.na(values(Temp1))==F]#remove empty cells
fcell<-intersect(ncell,tcell)#get cells in common
#get pet per cell 
for (i in fcell) {
  pet[i] <- th(temps[i], latt[i])
}
####################compute water balance for each cell
Wbmap<-Precip-stack(pet)
######
precipmapsum<-stackApply(stack(Precip), indices=rep(1,12), fun=sum , na.rm=F)
wbmapmean<-stackApply(stack(Wbmap), indices=rep(1,12), fun=mean , na.rm=TRUE)
wbmapsum<-stackApply(stack(Wbmap), indices=rep(1,12), fun=sum , na.rm=F)
wbmapmean2<- focal(wbmapmean,w=matrix(1/9, nc=3, nr=3), fun=fillNoData,NAonly=T,pad=T,na.rm=T)
#
#73%signf signif
#22% no signif change
#5%signif increase
#############################save spei and wb maps for plotting later and keep spei dates
save.image('./Khoury&Commes/processed/step2_SPEI&waterbalance_Mapsdata.RData')
keep(datenames,sure=T)
#############################################################################################

load("./Khoury&Commes/processed/step1_NDVIdata.RData")
########################################################################################################
##extraction to SPEIs and all scale for time-series with NDVI and LAI
##creating dataframes for each SPEI scale (1:48) with pixel time-series 
for (i in c(1:48)){
  n1 <- paste('./Khoury&Commes/data/Spanish SPEI_1.1km_48scales',i,'.nc',sep='')# spei dataset location
  d1 <- nc_open(n1)
  rbrick1 <- brick(n1,varname="value",level=1)
  proj4string(rbrick1)<-CRS('+proj=utm +zone=30 +ellps=intl +towgs84=-131,-100.3,-163.4,-1.244,-0.02,-1.144,9.39 +units=m +no_defs')
  rlayer<-subset(rbrick1,datenames)
  spei <- extract(rlayer, Dbestdata2, method='bilinear')
  spei <- data.frame(spei)
  spei<-t(spei)
  spei<-data.frame(spei)
  ind<-rownames(spei)
  ind<-as.Date(ind,format='X%Y.%m.%d' )
  zoospei<-zoo(spei,ind)                                                                                      
  colnames(zoospei)<-Dbestdata2$ID
  zoospei<-na.approx(zoospei, rule=2)
  zoospei<-aggregate(zoospei,as.yearmon(index(zoospei)),mean)
  speinum<-paste('spei',i,sep='')#lag
  speinum<-assign(speinum,zoospei)
}
#################################################################################################
#save.image('./Khoury&Commes/processed/step2_SPEIdata.RData')
