##################################################################
#Khoury,S. and Coomes, D.A. (2020) Resilience of Spanish forests 
#to recent droughts and climate change. Global change biology. 
#Resilience metrics analysis alongside environmental data 
#using spatial auto-correlation regression models.
#################################################################
library(spdep)
library(spatialreg)
library(smatr)
library(doParallel)
library(foreach)
library(tidyverse)
#useful function 
'%!in%' <- function(x,y)!('%in%'(x,y))
###############################
#load("./Khoury&Coomes/processed/step5_NDVIandfactorData.RData")#ndvimetrics,laimetrics,spain
load("./Khoury&Coomes/processed/step3.5_regressionDataWithWB.RData")
#add important variables
ndvimetrics$time18<-ndvimetrics$time*215#change in greenness over the 18 years (ndvi unit/18years)
ndvimetrics<-merge(ndvimetrics,ndvimetrics_wb[,c('trendwb','WB')],by='ID')# get water balance change from waterbalance analysis
ndvimetrics$trendspei18<-ndvimetrics$trendspei1*215#change in spei over the 18 years (ndvi unit/18years)
ndvimetrics$trendwb18<-ndvimetrics$trendwb*215#change in greenness over the 18 years (mm/18years)

laimetrics$time18<-laimetrics$time*215#change in LAI over the 18 years (total m2/m2/18years)
laimetrics$trendspei18<-laimetrics$trendspei1*215#change in spei over the 18 years (spei unit/18years)

###############################set contrast used for analysing forest species and protected area effects
attr(ndvimetrics$domspp3,'contrasts')<-'contr.sum'(10)#species being compared to the average of all species
attr(ndvimetrics$pa,'contrasts')<-'contr.treatment'(3)#protected areas being compared to non-protected areas
attr(laimetrics$domspp3,'contrasts')<-'contr.sum'(10)#same as with NDVI
attr(laimetrics$pa,'contrasts')<-'contr.treatment'(3)#same as with NDVI

###################list of models being evaluated (listed in supporting information)
###list of models to run
lmlist<-c(trendspeitot~WB,
          trendspei18~trendwb18,
          ndvimean ~ WB,
          cci~ scale(WB)+scale(elevation),
          time18~scale(WB)+scale(elevation)+ scale(ndvimean)+domsp+ pa#+lcl,
          ndvimean~scale(WB)+scale(elevation)+ domspp+ pa#+lcl,
          time18~scale(WB)+scale(elevation)+scale(ndvimean)+ domspp3+ pa,
          ndvimean~scale(WB)+scale(elevation)+ domspp3+ pa
          cci~ scale(WB)+scale(elevation)+ domspp+ pa,#+lcl,
          cci~ scale(WB)+scale(elevation)+ domspp3+ pa,#+lcl,
          cgw~ scale(ndvimean)+scale(elevation), 
          cgw~ scale(ndvimean)+scale(elevation)+domspp+ pa,#+lcl, 
          cgw~ scale(ndvimean)+scale(elevation)+domspp3+ pa,#+lcl, 
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean),
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean)+ domspp+ pa,#+lcl,
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean)+ domspp3+ pa,#+lcl,
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean),
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean)+ domspp+ pa,#+lcl,
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean)+ domspp3+ pa,#+lcl,
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean),
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean) + domspp+ pa,#+lcl,#+ domspp+ pa#15
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(ndvimean) + domspp3+ pa#+lcl
          )
lailist<-c(laimean ~ WB,#35
          cci~ scale(WB)+scale(laimean),
          cci~ scale(WB)+scale(laimean)+ domspp+ pa,#+lcl,
          cci~ scale(WB)+scale(laimean)+ domspp2+ pa,#+lcl,
          cgw~ scale(laimean)+scale(elevation), 
          cgw~ scale(laimean)+scale(elevation)+domspp+ pa,#+lcl, 
          cgw~ scale(laimean)+scale(elevation)+domspp2+ pa,#+lcl, 
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean),
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp+ pa,#+lcl,
          log(senstv)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp2+ pa,#+lcl,
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean),
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp+ pa,#+lcl,
          log(recov)~ scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp2+ pa,#+lcl,
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean),
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp+ pa,#+lcl,#+ domspp+ pa#15
          log(resil)~scale(WB)+scale(minspei)+scale(maxspei)+scale(elevation)+scale(laimean)+ domspp2+ pa,#+lcl
          time~pa,
          decl~pa
          )

#############function that runs auto-correlation regression on models using different neighborhood distances keeping the best AIC model which minimise spatial auto-correlation ####################
sarlm_testvario<-function(data,form=form, ncores = 2){
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  dis<-c(15,20,25,30,35)#neighborhood distance in km with 5 km increment being considered here after checking residual auto-correlation with pgrimess Moran's I
  myres<-foreach(
    i = dis,
    # i = 1:2,
    .packages = c("spdep","spatialreg", "tidyverse"),
    .inorder = TRUE,
    .errorhandling = "pass"
  ) %dopar% { 
    neighbours.dist<-dnearneigh(as(data,'SpatialPoints'),d1=0,d2=i)
    neighbours.dist.w<-nb2listw(neighbours.dist,style="W",zero.policy = TRUE)
    testmodel<- errorsarlm(formula=as.character(form),data = data,tol.solve=1*10^(-34),listw = neighbours.dist.w,zero.policy = TRUE, quiet = FALSE,method = "eigen")
   # return(AIC(testmodel))
    l<-list(summary(testmodel,Nagelkerke=T)$N,testmodel,AIC(testmodel))
    return(l)
    }
  stopCluster(cl)
  minaic<-min(myres[[1]][[3]],myres[[2]][[3]],myres[[3]][[3]],myres[[4]][[3]],myres[[5]][[3]]) 
  bol<-c(myres[[1]][[3]],myres[[2]][[3]],myres[[3]][[3]],myres[[4]][[3]],myres[[5]][[3]])==minaic
  myres<-myres[bol]
  return(myres)
}
#testrun<-sarlm_testvario(data=testdf,form=lmlist[2])
#print("modelling complete")
#qqnorm(y=scale(testrun[[1]][[2]]$residuals))# check normality of residuals
#qqline(y=scale(testrun[[1]][[2]]$residuals)) 
######### run function on models
for (i in c(1:length(lmlist))){
  char<-str_split_fixed(as.character(lmlist[i]),'~',2)[1]
  name<-paste('mods',char,i,sep='')
  modelres<-sarlm_testvario(data=ndvimetrics,form=lmlist[i])
  name<-assign(name,modelres)
}
###################################################################################################################
save.image('./Khoury&Coomes/processed/step6_spatialAutocorrelationModelsResults.RData')
