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
ires='100m'
fscasource='aso'#needs to be aso if ires < 500m
pathout='output/xdate-modeling'
dir.create(pathout,recursive=TRUE)


## aso swe ----
asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)
aso_stack=stack(asofn)
asoswedates=sapply(strsplit(names(aso_stack),'[.\\_]'),'[',3)

asoswe=as.data.frame(aso_stack,xy=T) %>%
  tbl_df %>%
  gather(dte,swe,-x,-y) %>%
  separate(dte,into=c('basin','res','dte')) %>%
  filter(!is.na(swe))

# which asodate should be used for each simdate?  Calculate from splitsample
# bestasodte=read_csv('output/splitsample-modeling/bestasodate.tsv')

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
phvdata=get_phvdata(geoarea2,ires,phv_stack,asoswe,fscasource) %>%
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

  alldata=inner_join(phvdata,adrast,by=c('x','y','yr'))
} else {
  alldata=phvdata
}


cl <- parallel::makePSOCKcluster(4)
allphvmdls <-
  alldata %>%
  rename(mdldte=dte) %>%
  group_by(mdldte) %>%
  nest() %>%
  mutate(
    phv_obj_glmmdl=map(data,gnet_phv,cl),
    phv_aug_glmmdl=map(phv_obj_glmmdl,augmentEnet,alldata),
    phv_coef_glmmdl=map(phv_obj_glmmdl,coef),
    #
    phvfsca_obj_glmmdl=map(data,gnet_phvfsca,cl),
    phvfsca_aug_glmmdl=map(phvfsca_obj_glmmdl,augmentEnet,alldata),
    phvfsca_coef_glmmdl=map(phvfsca_obj_glmmdl,coef)
  ) %>%
  dplyr::select(-contains('obj'))

print('phv models finished.')

saveRDS(alldata,paste0('output/alldata_',ires,'.rds'))
saveRDS(allphvmdls %>% dplyr::select(contains('aug'),-data), paste0(pathout,'phvmdls_augment_',ires,'.rds'))
saveRDS(allphvmdls %>% dplyr::select(contains('coef'),-data), paste0(pathout,'phvmdls_coef_',ires,'.rds'))

# asoswe=asoswe %>% filter(asodte %in% unique(asodte)[5:9])
# allphvmdls <- allphvmdls %>% filter(.id=='01'|.id=='02')

## need to join aso swe data for asoswe models. doing phv separately before so it isn't fitting models on sets of the same data
allmdldata <- allphvmdls %>% select(mdldte,data)
rm(allphvmdls)

allasomdls <-
  allmdldata %>%
  group_by(mdldte) %>%
  mutate(
    phvaso_obj_glmmdl=map(data,gnet_phvaso,join_asodata,asoswe,cl,pickASOdte),
    phvasofsca_obj_glmmdl=map(data,gnet_phvasofsca,join_asodata,asoswe,cl,pickASOdte)
  ) %>%
  unnest(phvaso_obj_glmmdl,phvasofsca_obj_glmmdl,.drop=FALSE)

if(!identical(allasomdls[[5]],allasomdls[[7]])) stop()

print('phvaso models finished.')

allasomdls <-
  allasomdls[!duplicated(as.list(allasomdls))] %>%
  mutate(
    phvaso_aug_glmmdl=map(phvaso_obj_glmmdl,augmentEnet,alldata),
    phvaso_coef_glmmdl=map(phvaso_obj_glmmdl,coef),
    phvasofsca_aug_glmmdl=map(phvasofsca_obj_glmmdl,augmentEnet,alldata),
    phvasofsca_coef_glmmdl=map(phvasofsca_obj_glmmdl,coef)
  )

parallel::stopCluster(cl)

saveRDS(allasomdls %>% select(-contains('coef'),-contains('obj')), paste0(pathout,'asomdls_augment_',ires,'.rds'))
saveRDS(allasomdls %>% select(-contains('aug'),-contains('obj')), paste0(pathout,'asomdls_coef_',ires,'.rds'))
