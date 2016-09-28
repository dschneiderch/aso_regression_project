library(raster)
library(rgdal)
library(tidyverse)
library(cowplot)
library(viridis)
library(broom)
library(modelr)
library(mgcv)
library(glmnetUtils)
library(rpart)

source('scripts/functions_modelfitting.R')

geoarea2='tuo'
ires='500m'
fscasource='aso'#needs to be aso if ires < 500m
pathout='output/splitsample-modeling/'
dir.create(pathout,recursive=TRUE)

## aso swe ----
asoswe=get_asoswe(geoarea2,ires)
asoswedates=unique(asoswe$dte)
# testdates=unique(asoswe$dte)[c(1,2,8,9)]
# asoswe %>% filter(dte %in% testdates) %>% group_by(dte) %>% summarise(min(swe),max(swe))
# # asoswe1=filter(asoswe,dte %in% testdates)
# # unique(asoswe1$dte)
# # range(asoswe1$swe)
# asoswe <- asoswe %>% filter(dte %in% testdates)
# unique(asoswe$dte)
# range(asoswe$swe)

## phv variables ----
elev_rast=raster(paste0('data/elevation/',geoarea2,'_dem',ires,'.tif'))
zness_rast=raster(paste0('data/elevation/',geoarea2,'_zness',ires,'.tif'))
northness_rast=raster(paste0('data/elevation/',geoarea2,'_northness',ires,'.tif'))
eastness_rast=raster(paste0('data/elevation/',geoarea2,'_eastness',ires,'.tif'))

tpi_rast=raster(paste0('data/elevation/',geoarea2,'_tpi',ires,'.tif'))
# tri_rast=raster(paste0('data/elevation/',basin,'_tri',ires,'.tif'))
vrm_rast=raster(paste0('data/elevation/',geoarea2,'_vrm1cell_dem',ires,'.tif'))

northing_rast=raster(paste0('data/elevation/',geoarea2,'_northing',ires,'.tif'))
easting_rast=raster(paste0('data/elevation/',geoarea2,'_easting',ires,'.tif'))

vegheight_rast=raster(paste0('data/veg/',geoarea2,'_vegheight',ires,'.tif'))

stdslope_rast=raster(paste0('data/elevation/',geoarea2,'_stdslope1cell_dem',ires,'.tif'))
# ---

phv_stack=stack(elev_rast,zness_rast,stdslope_rast,vrm_rast,northness_rast,eastness_rast,northing_rast,easting_rast,vegheight_rast)
phvdata=get_phvdata(geoarea2,ires,phv_stack,asoswe,asoswedates,fscasource) %>%
  filter(complete.cases(.)) %>%
  mutate(yr=substr(dte,1,4))

# change names for iterating over patterns
asoswe <- asoswe %>%
  rename(asoswe=swe,
         asodte=dte)

sens=0
# for(sens in seq_len(7)){

if(fscasource=='modscag'){
  imageyrs=unique(phvdata$yr)
  avgdate_rast=stack(paste0('data/snowoffdate/tuo_snowoff_X',imageyrs,'.tif'))
  adrast=as.data.frame(avgdate_rast,xy=T) %>%
    tbl_df %>%
    gather(yr,snowoff,-x,-y) %>%
    # separate(yr,into=c('geoarea2','var','yr')) %>%
    mutate(yr=as.character(readr::parse_number(yr)))
  #
  alldata=inner_join(phvdata,adrast,by=c('x','y','yr'))
} else {
  alldata=phvdata
}

# Calculate the number of cores
num_cores <- parallel::detectCores() - 1
#if(num_cores>22) num_cores=22
# Initiate cluster
#if not fork then need to export libraries and variables to nodes. unix only
cl <- parallel::makeCluster(num_cores,type='FORK')

allphvmdls <-
  alldata %>%
  group_by(dte) %>%
  nest() %>%
  mutate(straps=map(data,crossv_mc,60)) %>%
  unnest(straps) %>%
  mutate(
    phv_obj_glmmdl=map(train,gnet_phv,cl=cl),
    phv_aug_glmmdl=map2(phv_obj_glmmdl,test,augmentEnet),
    phv_coef_glmmdl=map(phv_obj_glmmdl,coef),
    #
    phvfsca_obj_glmmdl=map(train,gnet_phvfsca,cl),
    phvfsca_aug_glmmdl=map2(phvfsca_obj_glmmdl,test,augmentEnet),
    phvfsca_coef_glmmdl=map(phvfsca_obj_glmmdl,coef)
  )

print('phv models finished.')

saveRDS(alldata,paste0(pathout,'alldata_',ires,'.rds'))
saveRDS(allphvmdls %>% select(-contains('coef'),-contains('obj')), paste0(pathout,'phvmdls_augment_',ires,'.rds'))
saveRDS(allphvmdls %>% select(-contains('aug'),-contains('obj')), paste0(pathout,'phvmdls_coef_',ires,'.rds'))

#asoswe=asoswe %>% filter(asodte %in% unique(asodte)[5:6])
# allphvmdls <- allphvmdls %>% filter(.id=='01'|.id=='02')
## need to join aso swe data for asoswe models. doing phv separately before so it isn't fitting models on sets of the same data

allmdldata <- allphvmdls %>%
  select(dte,.id,train,test)  # %>% filter(dte==unique(dte)[1:2])
# rm(allphvmdls)

allasomdls <-
  allmdldata %>%
  group_by(dte,.id) %>%
  mutate(#need to join asodata separately because it takes too muchmemory to join with all train elements first
    phvaso_obj_glmmdl=map(train,gnet_phvaso,join_asodata,asoswe,cl),
    phvasofsca_obj_glmmdl=map(train,gnet_phvasofsca,join_asodata,asoswe,cl)
  ) %>%
  unnest(phvaso_obj_glmmdl,phvasofsca_obj_glmmdl,.drop=FALSE)

if(!identical(allasomdls[[5]],allasomdls[[7]])) stop()

print('phvaso models finished.')

allasomdls <-
  allasomdls[!duplicated(as.list(allasomdls))] %>%
  mutate(
    phvaso_aug_glmmdl=pmap(list(phvaso_obj_glmmdl,test,asodte),augmentEnet,asoswe),
    phvaso_coef_glmmdl=map(phvaso_obj_glmmdl,coef),
    phvasofsca_aug_glmmdl=pmap(list(phvasofsca_obj_glmmdl,test,asodte),augmentEnet,asoswe),
    phvasofsca_coef_glmmdl=map(phvasofsca_obj_glmmdl,coef)
  )

parallel::stopCluster(cl); cl=NULL

allasomdls <-
  allasomdls %>%
  mutate(yr=substr(dte,1,4))

for(iyr in unique(allasomdls$yr)){
  saveRDS(allasomdls %>% select(dte,yr,.id,asodte,phvaso_aug_glmmdl) %>% filter_(~yr==iyr), paste0(pathout,'phvasomdls_augment_',iyr,'_',ires,'.rds'))
  saveRDS(allasomdls %>% select(dte,yr,.id,asodte,phvaso_coef_glmmdl) %>% filter_(~yr==iyr), paste0(pathout,'phvasomdls_coef_',iyr,'_',ires,'.rds'))
  #
  saveRDS(allasomdls %>% select(dte,yr,.id,asodte,phvasofsca_aug_glmmdl) %>% filter_(~yr==iyr), paste0(pathout,'phvasofscamdls_augment_',iyr,'_',ires,'.rds'))
  saveRDS(allasomdls %>% select(dte,yr,.id,asodte,phvasofsca_coef_glmmdl) %>% filter_(~yr==iyr), paste0(pathout,'phvasofscamdls_coef_',iyr,'_',ires,'.rds'))
}

allasomdls <- allasomdls %>% dplyr::select(dte,yr,.id,asodte,contains('aug'))

# phvaso models error statistics
## calculate rmse and pct rmse to determine best aso date for each date simulated ----
if(exists('allasomdls')){#'if/else' mostly here in case r crashes after saving
  phvasomdls=allasomdls %>% dplyr::select(-phvasofsca_aug_glmmdl)
} else {
  fns=dir('output/splitsample-modeling',pattern='phvasomdls_aug',full.names=TRUE)
  phvasomdls=bind_rows(lapply(fns,readRDS))
}

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
if(exists('allasomdls')){
  phvasofscamdls=allasomdls %>% dplyr::select(-phvaso_aug_glmmdl)
} else {
  fns=dir('output/splitsample-modeling',pattern='phvasofscamdls_aug',full.names=TRUE)
  phvasofscamdls=bind_rows(lapply(fns,readRDS))
}

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


## calculate r2 and percent mae ----
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
saveRDS(predictions,paste0('output/splitsample-modeling/stats_combined_splitsample',ires,'.rds'))

# }    # <------ uncomment to run sensitivity
