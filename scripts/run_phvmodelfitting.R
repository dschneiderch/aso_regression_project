library(raster)
library(rgdal)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(broom)
library(purrr)
library(modelr)
library(mgcv)
library(glmnet)
library(glmnetUtils)
source('scripts/functions_modelfitting.R')

geoarea2='tuo'
ires='500m'

asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)
aso_stack=stack(asofn)
asoswedates=sapply(strsplit(names(aso_stack),'[.\\_]'),'[',3)

## phv variables
elev_rast=raster('data/elevation/tuo_dem500m.tif')
zness_rast=raster('data/elevation/tuo_zness500m.tif')
northness_rast=raster('data/elevation/tuo_northness500m.tif')
eastness_rast=raster('data/elevation/tuo_eastness500m.tif')

tpi_rast=raster('data/elevation/tuo_tpi500m.tif')
vrm_rast=raster('data/elevation/tuo_vrm1cell_dem500m.tif')

northing_rast=raster('data/elevation/tuo_northing500m.tif')
easting_rast=raster('data/elevation/tuo_easting500m.tif')

vegheight_rast=raster('data/veg/tuo_vegheight_500m.tif')

stdslope_rast=raster('data/elevation/tuo_stdslope1cell_dem500m.tif')
## ---

phv_stack=stack(elev_rast,zness_rast,stdslope_rast,vrm_rast,northness_rast,eastness_rast,northing_rast,easting_rast,vegheight_rast)
phvdata=get_phvdata(geoarea2,ires,phv_stack)



phvdata <- phvdata %>% filter(complete.cases(.))

foldid=sample(1:10,size=nrow(phvdata),replace=T)
y=data.matrix(phvdata[['swe']])
x=data.matrix(phvdata[,c(7:15)])
cv1=cv.glmnet(x,y,type.measure='mse',foldid=foldid,alpha=1,parallel = TRUE)
cv.5=cv.glmnet(x,y,type.measure='mse',foldid=foldid,alpha=.5,parallel = TRUE)
cv0=cv.glmnet(x,y,type.measure='mse',foldid=foldid,alpha=0,parallel=TRUE)
par(mfrow=c(2,2))
plot(cv1);plot(cv.5);plot(cv0)
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))

phvdata$cv0pred=predict(cv0,newx=x)
phvdata$cv1pred=predict(cv1,newx=x)

pctmae=function(yobs,yhat){
  mean(abs(yhat-yobs),na.rm=T)/mean(yobs,na.rm=T)
}

lambda_ind=which(cv1$lambda==cv1$lambda.1se)
cv1$cvm[lambda_ind]
lambda_ind=which(cv0$lambda==cv0$lambda.1se)
cv0$cvm[lambda_ind]
lambda_ind=which(cv.5$lambda==cv.5$lambda.1se)
cv.5$cvm[lambda_ind]

cv1l=cv1$lambda.1se
cv0l=cv0$lambda.1se
cv.5l=cv.5$lambda.1se

gnet1=glmnet(x,y,alpha=1)
gnet0=glmnet(x,y,alpha=0)
gnet.5=glmnet(x,y,alpha=0.5)
plot(gnet1,label=T)
coef(gnet1,s=cv1l)

plot(gnet0,label=T)
coef(gnet0,s=cv0l)

plot(gnet.5,label=T)
coef(gnet.5,s=cv.5l)


graphics.off()
plot(cv1)
coef(cv1)
coef(cv0)
phvdata %>%
  summarise(
    pctmae0=pctmae(swe,cv0pred),
    pctmae1=pctmae(swe,cv1pred)
  )


gnetfitting=function(dF){
  cvfit=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight,data=dF,type.measure='mse',alpha=0,parallel = TRUE)
  return(cvfit)
}
gnetcoefs=function(mdl){
  coef(mdl)
}

raster::beginCluster()
allmdls <-
  phvdata %>%
  filter(complete.cases(.)) %>%
  group_by(dte) %>%
  do({
    crossv_mc(.,2) %>%
      mutate(
        glmnet_mdl=map(train, gnetfitting),
        glmnet_myr2=map2_dbl(glmnet_mdl,test,myr2),
        glmnet_pctmae=map2_dbl(glmnet_mdl,test,mypctmae),
        glmnet_coefs=map(glmnet_mdl,gnetcoefs)
      )
  })
raster::endCluster()

allmdls %>%
  gather(metricvar,metricval,glmnet_myr2,glmnet_pctmae) %>%
ggplot(.)+
  geom_boxplot(aes(x=dte,y=metricval))+
  facet_wrap(~metricvar,scales='free')+
  theme(axis.text.x=element_text(angle=45,hjust=1))
