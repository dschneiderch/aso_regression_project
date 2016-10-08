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
pathout='output/asooptimize-modeling/'
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
cl <- NULL
cl <- parallel::makeCluster(num_cores,type='FORK')

nsamples=c(50,100,500,1000)
# nsamp=nsamples[1]
# plusallphvmdls=data_frame()
# for(nsamp in nsamples){
# optdata <-
#   alldata %>%
#   # filter(dte==unique(dte)[1:3]) %>%
#   mutate(dte2=dte) %>%
#   group_by(dte2) %>%
#   nest() %>%
#   mutate(straps=map(data,crossv_mc,20,test=.2)) %>%
#   unnest(straps) %>%
#   mutate(trainopt=map2(train,.id,function(x,.id) {
#     xdf=as.data.frame(x) %>% mutate(.id)
#     # print(fr*nr)
#     list(replicate(20,sample_n(xdf,nsamp),simplify=FALSE))
#       } %>% flatten
#     ))
#
#   allphvmdls <-
#     map_df(optdata$trainopt,unlist,rec=F) %>% dplyr::select(dte,.id,everything()) %>%
#     left_join(optdata %>% select(dte2,.id),by=c('dte'='dte2','.id')) %>%
#     group_by(dte,.id) %>%
#     nest(.key='train') %>%
#     left_join(optdata %>% select(dte2,.id,test),by=c('dte'='dte2','.id')) %>%
#   mutate(
#     phv_obj_glmmdl=map(train,gnet_phv,cl=cl),
#     phv_aug_glmmdl=map2(phv_obj_glmmdl,test,augmentEnet),
#     phv_coef_glmmdl=map(phv_obj_glmmdl,coef),
#     #
#     phvfsca_obj_glmmdl=map(train,gnet_phvfsca,cl),
#     phvfsca_aug_glmmdl=map2(phvfsca_obj_glmmdl,test,augmentEnet),
#     phvfsca_coef_glmmdl=map(phvfsca_obj_glmmdl,coef),
#     #
#     nsamp
#   )
#
#   plusallphvmdls <- bind_rows(plusallphvmdls,allphvmdls)
#   print(nsamp)
# }
# rm(allphvmdls)
# print('phv models finished.')
#
# saveRDS(alldata,paste0(pathout,'alldata_',ires,'.rds'))
# saveRDS(plusallphvmdls %>% select(-contains('coef'),-contains('obj')), paste0(pathout,'phvmdls_augment_',ires,'.rds'))
# saveRDS(plusallphvmdls %>% select(-contains('aug'),-contains('obj')), paste0(pathout,'phvmdls_coef_',ires,'.rds'))

alldata <- readRDS(paste0(pathout,'alldata_',ires,'.rds'))
plusallphvmdls <- readRDS(paste0(pathout,'phvmdls_augment_',ires,'.rds'))

#asoswe=asoswe %>% filter(asodte %in% unique(asodte)[5:6])
# allphvmdls <- allphvmdls %>% filter(.id=='01'|.id=='02')
## need to join aso swe data for asoswe models. doing phv separately before so it isn't fitting models on sets of the same data

allmdldata <- plusallphvmdls %>%
  select(dte,.id,nsamp,train,test) %>%
  mutate(yr=substr(dte,1,4))# %>% filter(dte==unique(dte)[1:2])
rm(plusallphvmdls)

iyr=2013
for(iyr in unique(allmdldata$yr)){
if(iyr==2013) next
print(iyr)
  allasomdls <-
    allmdldata %>%
    filter_(~yr==iyr) %>%
    group_by(dte,yr,.id,nsamp) %>%
    mutate(#need to join asodata separately because it takes too muchmemory to join with all train elements first
      phvaso_obj_glmmdl=map(train,gnet_phvaso,join_asodata,asoswe,cl),
      phvasofsca_obj_glmmdl=map(train,gnet_phvasofsca,join_asodata,asoswe,cl)
    ) %>%
    unnest(phvaso_obj_glmmdl,phvasofsca_obj_glmmdl,.drop=FALSE)

  # if(!identical(allasomdls[[5]],allasomdls[[7]])) stop()

  print('phvaso models finished.')
print(pryr::mem_used())
  allasomdls <-
    allasomdls[!duplicated(as.list(allasomdls))] %>%
    mutate(
      phvaso_aug_glmmdl=pmap(list(phvaso_obj_glmmdl,test,asodte),augmentEnet,asoswe),
      phvaso_coef_glmmdl=map(phvaso_obj_glmmdl,coef),
      phvasofsca_aug_glmmdl=pmap(list(phvasofsca_obj_glmmdl,test,asodte),augmentEnet,asoswe),
      phvasofsca_coef_glmmdl=map(phvasofsca_obj_glmmdl,coef)
    )
print(pryr::mem_used())
  parallel::stopCluster(cl); cl=NULL

  saveRDS(allasomdls %>% select(dte,yr,.id,nsamp,asodte,phvaso_aug_glmmdl), paste0(pathout,'phvasomdls_augment_',iyr,'_',ires,'.rds'))
  saveRDS(allasomdls %>% select(dte,yr,.id,nsamp,asodte,phvaso_coef_glmmdl), paste0(pathout,'phvasomdls_coef_',iyr,'_',ires,'.rds'))
  #
  saveRDS(allasomdls %>% select(dte,yr,.id,nsamp,asodte,phvasofsca_aug_glmmdl), paste0(pathout,'phvasofscamdls_augment_',iyr,'_',ires,'.rds'))
  saveRDS(allasomdls %>% select(dte,yr,.id,nsamp,asodte,phvasofsca_coef_glmmdl), paste0(pathout,'phvasofscamdls_coef_',iyr,'_',ires,'.rds'))

}
