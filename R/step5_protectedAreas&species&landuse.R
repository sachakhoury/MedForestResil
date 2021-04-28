##################################################################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Step5: Distinguishing protected areas from non-protected areas. 
#Extracting landcover types from Corine dataset 1990, 2000 and checking for land cover change. Looking at species groups.
#################################################################
library(geosphere)
library(maptools)
library(stringr)
library(raster)
library(rgdal)
library(spatialEco)
##################Loading shapefiles and rasters#####################################
#National borders
load("./Khoury&Coomes/data/Europe.rda")
spain=Europe[Europe$name=="Spain", ]
spain=as(spain,'SpatialPolygons')
spain<-spTransform(spain,wgs)
#red natura borders
rn2000<-readOGR('./protected areas/rn2000/RedNatura.shp')
rn2000<-spTransform(rn2000,ndvimetrics@proj4string)
enp<-readOGR('./Khoury&Coomes/data/protected areas/enp/enp.shp')
enp<-spTransform(enp,ndvimetrics@proj4string)
############################################get those protecte areas##############
pixels<-ndvimetrics[,1]
pixels$rn<-extract(rn2000,pixels)
pixels2<-extract(enp,pixels,fun=min)
pixels2<-pixels2[,c(1,7,9)]
pixels2$YEAR<-as.numeric(pixels2$YEAR)
pixels2<-aggregate(pixels2, by=list(pixels2$point.ID),FUN=min)
pixels$enp<-pixels2$ENP_FGRAL
pixels$rn2000<-pixels$rn$REDNATURA
pixels$enp<-as.numeric(as.character(pixels$enp))
pixels[is.na(pixels$enp),3]<-0
pixels$rn2000<-as.numeric(as.character(pixels$rn2000))
pixels[is.na(pixels$rn2000),4]<-0
pixels$pa<-0
pixels[pixels$rn2000==1,5]<-1
pixels[pixels$rn2000==2,5]<-2
pixels[pixels$rn2000==3,5]<-2#intersection
pixels$pa<-as.factor(pixels$pa)
summary(pixels$pa)#
ndvimetrics<-merge(ndvimetrics, pixels@data[,c(1,5)], by='ID',all.x=T,all.y=T)
ndvimetrics<-merge(ndvimetrics, pixels@data[,c(1,5)], by='ID',all.x=T,all.y=T)
###
keep(ndvimetrics,ndvimetrics,wgs,spain,sure=T)


##get the species again for more precision and sort in descending order from most present in spain to least as a factor.
domtree<-ndvimetrics[,1]
  
ndvimetrics$domspp<-as.character(extract(treesp, domtree, method='simple',factors=T))
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='20']<-'10'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='17']<-'10'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='8']<-'10'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='14']<-'11'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='1']<- 'Abies spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='10']<-'Other Broadleaves'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='5']<- 'Castanea spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='6']<- 'Eucalyptus spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='11']<-'Other Conifers'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='7']<- 'Fagus spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='8']<- 'Fraxinus domspp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='14']<-'Picea spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='15']<-'P. pinaster'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='12']<-'Other Pinus spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='13']<-'Other Quercus spp.'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='19']<-'Q. robur & Q. petraea'
ndvimetrics@data$domspp[ndvimetrics@data$domspp=='16']<-'P. sylvestris'

ndvimetrics@data$domspp<-as.factor(ndvimetrics@data$domspp,levels=c())  

summary(ndvimetrics$domspp)
###########################
laimetrics<-merge(laimetrics,ndvimetrics[,c('ID','domspp')],by='ID')


speciesnames<-c("Other Conifers","Other Pinus spp.","P. pinaster","Other Broadleaves","Eucalyptus spp.",
                "Castanea spp.","Other Quercus spp.","Abies spp.","P. sylvestris","Fagus spp.","Q. robur & Q. petraea")


ndvimetrics<-ndvimetrics[ndvimetrics@data$domspp!='Abies spp.']
########################## Get previous land cover and its change ####################################
ndvi2013<-raster('./Khoury&Coomes/data/ndvi_lai_waterbalance_maps/ndvi2013.tif')
corine2000<-raster('./Khoury&Coomes/data/clc2000/reproj2000.tif')
corine1990<-raster('./Khoury&Coomes/data/clc1990/reproj1990.tif')
resampcorine2000<-raster::resample(x=corine2000,y=ndvi2013,method='ngb')
resampcorine1990<-raster::resample(x=corine1990,y=ndvi2013,method='ngb')
#plot(reprojcorine2000)
######extract land cover
corine2000pts<- raster::extract(resampcorine2000, as(ndvimetrics,'SpatialPoints'), method='simple',fun=Mode,na.rm=T)
corine1990pts<-raster:: extract(resampcorine1990, as(ndvimetrics,'SpatialPoints'), method='simple',fun=Mode,na.rm=T)
######legend
legendclc<-read.delim('./Khoury&Coomes/data/clc1990/Legend/CLC1990_CLC1990_V2018_20_QGIS.txt', header = F, sep = ",", dec = ".")
colnames(legendclc)<-c('class#','r','g','b','alpha','classname')
######replace numbers by legend
corine2000pts2<-data.frame(corine2000pts,corine2000pts)
for (i in corine2000pts2$corine2000pts){corine2000pts2[corine2000pts2$corine2000pts==i,2]<-as.character(legendclc[legendclc$`class#`==i,'classname'])}
corine1990pts2<-data.frame(corine1990pts,corine1990pts)
for (i in corine1990pts2$corine1990pts){corine1990pts2[corine1990pts2$corine1990pts==i,2]<-as.character(legendclc[legendclc$`class#`==i,'classname'])}
#####turn land cover to factors
corine2000pts2$corine2000pts.1<-as.factor(corine2000pts2$corine2000pts.1)
corine1990pts2$corine1990pts.1<-as.factor(corine1990pts2$corine1990pts.1)
#######get summary
summary(corine1990pts2$corine1990pts.1)
summary(corine2000pts2$corine2000pts.1)
#####
corine1990pts2$ID<-ndvimetrics$ID
colnames(corine1990pts2)[2]<-'lc'
corine2000pts2$ID<-ndvimetrics$ID
colnames(corine2000pts2)[2]<-'lc'
###################################################
ndvimetrics<-merge(ndvimetrics,corine1990pts2,by='ID')
ndvimetrics$fil1<-F#loop to indicate whether pixels (that are not other broadleaves) are in agricutlural and other non-vegetation land cover in 1990
for (i in levels(ndvimetrics$domspp)){
  ndvimetrics@data[ndvimetrics$domspp==i,79]<-ndvimetrics@data[ndvimetrics$domspp==i,78]%in%levels(ndvimetrics$lc)[c(1,2,5,7,9,11,14,15,16,17,18,20,23,24,25)]}
summary(ndvimetrics$fil1)
ndvimetrics<-merge(ndvimetrics,corine2000pts2,by='ID')
ndvimetrics$fil2<-F #loop to indicate whether pixels (that are not other broadleaves) are in agricutlural and other non-vegetation land cover in 2000

for (i in levels(ndvimetrics$domspp)){
  ndvimetrics@data[ndvimetrics$domspp==i,82]<-ndvimetrics@data[ndvimetrics$domspp==i,81]%in%levels(ndvimetrics$lc.y)[c(1,2,5,7,8,11,12,16,17,18,19,20,22,25,26,27)]}
summary(ndvimetrics$fil2)
##filtering land covers
ndvimetrics$pafil<-ndvimetrics$pa==levels(ndvimetrics$pa)[1]# is pixel in non-protected area
ndvimetrics$sppfil<-ndvimetrics$domspp!=levels(ndvimetrics$domspp)[4]# is pixel other broadleaves category
ndvimetrics$deletecorine<-ndvimetrics$pafil+ndvimetrics$fil2+ndvimetrics$sppfil# column to filer non-broadleaved non-protected pixels that are agricuture in 2000
ndvimetrics<-ndvimetrics[ndvimetrics$deletecorine!=3,]#remove it
ndvimetrics$lcc<-(ndvimetrics$fil1+ndvimetrics$sppfil)==2#column to check if non-other broadleaves is agriculture in 1990
ndvimetrics[,c('fil1','sppfil','pafill','deletecorine','corine2000pts','corine1990pts')]<-NULL
summary(ndvimetrics$fil2)
summary(as.factor(ndvimetrics$lcc))
#[c(1,2,3,5,6,7,8,9,10)]#c(1,5,9,11,16,17,18,19,20,26)#[c(1,2,5,7,9,10,11,16,17,18,19,20,24,25,26)]
laimetrics<-merge(laimetrics,ndvimetrics[,c('ID','lc','lcc')],by='ID')
#########################
keep(ndvimetrics,laimetrics,spain,wgs,sure=T)
#####
save.image('./Khoury&Coomes/processed/step5_NDVIandfactorData.RData')
###########