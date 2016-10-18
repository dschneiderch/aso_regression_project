library(tidyverse)
library(cowplot)
library(viridis)
library(broom)
library(modelr)
library(spdep)

source('scripts/functions_modelfitting.R')

geoarea2='tuo'
ires='500m'
fscasource='aso'#needs to be aso if ires < 500m
pathin='output/splitsample-modeling/20iterations/'
pathout='output/splitsample-modeling/20iterations/'
dir.create(pathout,recursive=TRUE)

# functions
extractAug <- function(allaug,mdlpreds,mdltype){
  coln <- paste(mdlpreds,'aug',mdltype,sep='_')
  allaug %>%
    group_by(mdldte) %>%
    unnest_(coln) %>%
    # gather(var,val,swe:resid) %>%
    mutate(mdlpreds,mdltype)
}

r2 <- function(yobs,yhat){
  cor(yobs,yhat,use='complete.obs')^2
}
pctmae <- function(yobs,yhat){
  mean(abs(yhat-yobs),na.rm=T)/mean(yobs,na.rm=T)*100
}
pctbias <- function(yobs,yhat){
  mean(yhat-yobs,na.rm=T)/mean(yobs,na.rm=T)*100
}
calc_error <- function(dF){
  dF %>%
    summarise(
      r2=r2(swe,swehat),
      pctmae=pctmae(swe,swehat),
      pctbias=pctbias(swe,swehat)
    )
}

# ## ----- Check spatial autocorrelation of regional performance
# ## are the counts spatially correlated?
#
# allphvmdls <- readRDS(paste0(pathout,'phvmdls_augment_',ires,'.rds'))
# phvresid=allphvmdls %>%
#   rename(mdldte=dte) %>%
#   extractAug(.,'phv','glmmdl') %>%
#   mutate(resid=swehat-swe,
#          x=as.numeric(x),
#          y=as.numeric(y)) %>%
#   dplyr::select(mdldte,.id,x,y,resid) %>%
#   ungroup
#
# dFsp=phvresid %>% filter(.id=='01',mdldte=='20130403')
# est_moranI <- function(dFsp){
#     print(paste(dFsp$dte,dFsp$.id))
#     coordinates(dFsp)=~x+y
#     proj4string(dFsp)='+proj=utm +zone=11 +datum=WGS84'
#     kd=dnearneigh(dFsp,d1=0,d2=50000)
#     dist=nbdists(kd,dFsp)
#     p=1#fassnacht found 0.5 to be optimal until mid-march, then up to 1.3 at end-april then back down to .5.
#     idw=lapply(dist,function(x) 1/(x^p))
#     kd.idw=nb2listw(kd,glist=idw,style='B')
#     mt=moran.test(dFsp$resid,listw=kd.idw)
#     tidy(mt)
#   }
#
# tmp=phvresid %>% group_by(.id,mdldte) %>% nest %>% mutate(mt=map_df(data,est_moranI))
#
# phvresid %>% split(.$dte,.$id)


# phvaso models error statistics
## calculate rmse and pct rmse to determine best aso date for each date simulated ----
fns=dir(pathin,pattern='phvasomdls_aug',full.names=TRUE)
phvasomdls=bind_rows(lapply(fns,readRDS))

stats <-
  phvasomdls %>%
  # filter(dte==asodte) %>%
  group_by(dte,asodte,yr) %>%
  unnest(phvaso_aug_glmmdl) %>%
  summarise(r2=r2(swe,swehat),
            rmse=rootmse(swe,swehat),
            pctrmse=rmse/mean(swe,na.rm=T)*100)


readr::write_tsv(stats,path=paste0(pathout,'errors_allasodates_phvaso_',ires,'.txt'))

stats %>%
  mutate(asoyr=substr(asodte,1,4)) %>%
  filter(yr!=asoyr) %>%
  group_by(dte) %>%
  summarise(bestasodte=asodte[which.min(rmse)],
            bestrmse=min(rmse,na.rm=T),
            bestpctrmse=pctrmse[which.min(rmse)]) %>%
  readr::write_tsv(.,path=paste0(pathout,'bestasodates_phvaso_',ires,'.txt'))

## calculate r2 and percent mae ----
# fns <- dir(pathin,glob2rx(paste0('^phvasomdls_aug*',ires,'*.rds')),full.names=T)
phvaso <-
  # bind_rows(
  #   lapply(fns,readRDS)
  # ) %>%
  phvasomdls %>%
  rename(mdldte=dte) %>%
  extractAug(.,'phvaso','glmmdl') %>%
  group_by(mdldte,asodte,mdlpreds) %>%
  calc_error()
rm(phvasomdls)


# phvasofsca models error statistics
## calculate rmse and pct rmse to determine best aso date for each date simulated ----
fns=dir(pathin,pattern='phvasofscamdls_aug',full.names=TRUE)
phvasofscamdls=bind_rows(lapply(fns,readRDS))

stats <-
  phvasofscamdls %>%
  # filter(dte!=asodte) %>%
  group_by(dte,asodte,yr) %>%
  unnest(phvasofsca_aug_glmmdl) %>%
  summarise(rmse=rootmse(swe,swehat),
            pctrmse=rmse/mean(swe,na.rm=T)*100)

readr::write_tsv(stats,path=paste0(pathout,'errors_allasodates_phvasofsca_',ires,'.txt'))

stats %>%
  mutate(asoyr=substr(asodte,1,4)) %>%
  filter(yr!=asoyr) %>%
  group_by(dte) %>%
  summarise(bestasodte=asodte[which.min(rmse)],
            bestrmse=min(rmse,na.rm=T),
            bestpctrmse=pctrmse[which.min(rmse)]) %>%
  readr::write_tsv(.,path=paste0(pathout,'bestasodates_phvasofsca_',ires,'.txt'))

## calculate r2 and percent mae ----
# fns  <- dir(pathin,glob2rx(paste0('^phvasofscamdls_aug*',ires,'*.rds')),full.names=T)
phvasofsca <-
  # bind_rows(
  #   lapply(fns,readRDS)
  # ) %>%
  phvasofscamdls %>%
  rename(mdldte=dte) %>%
  extractAug(.,'phvasofsca','glmmdl')%>%
  group_by(mdldte,asodte,mdlpreds) %>%
  calc_error()
rm(phvasofscamdls)


## calculate r2 and percent mae for phv(fsca)----
# saveRDS(alldata,paste0(pathout,'alldata_',ires,'.rds'))
allphvmdls <- readRDS(paste0(pathout,'phvmdls_augment_',ires,'.rds'))

phv <-
  # readRDS(paste0(pathin,'phvmdls_augment_',ires,'.rds')) %>%
  allphvmdls %>%
  rename(mdldte=dte) %>%
  extractAug(.,'phv','glmmdl') %>%
  group_by(mdldte,mdlpreds) %>%
  calc_error()

phvfsca <-
  # readRDS(paste0(pathin,'phvmdls_augment_',ires,'.rds')) %>%
  allphvmdls %>%
  rename(mdldte=dte) %>%
  extractAug(.,'phvfsca','glmmdl') %>%
  group_by(mdldte,mdlpreds) %>%
  calc_error()

# combine r2 and pct mae
predictions=bind_rows(phv,phvfsca,phvaso,phvasofsca)
saveRDS(predictions,paste0(pathout,'stats_combined_',ires,'.rds'))
