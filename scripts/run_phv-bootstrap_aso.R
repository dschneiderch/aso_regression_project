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

asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)
aso_stack=stack(asofn)
asoswedates=sapply(strsplit(names(aso_stack),'[.\\_]'),'[',3)

## phv variables ---
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
phvdata=get_phvdata(geoarea2,ires,phv_stack) %>%
  filter(complete.cases(.)) %>%
  mutate(yr=substr(dte,1,4))

sens=0
# for(sens in seq_len(7)){

imageyrs=unique(phvdata$yr)
avgdate_rast=stack(paste0('data/snowoffdate/tuo_snowoff_X',imageyrs,'.tif'))
adrast=as.data.frame(avgdate_rast,xy=T) %>%
  tbl_df %>%
  gather(yr,snowoff,-x,-y) %>%
  # separate(yr,into=c('geoarea2','var','yr')) %>%
  mutate(yr=as.character(readr::parse_number(yr)))

alldata=inner_join(phvdata,adrast,by=c('x','y','yr'))

asoswe=as.data.frame(aso_stack,xy=T) %>%
  tbl_df %>%
  gather(asodte,asoswe,-x,-y) %>%
  separate(asodte,into=c('basin','res','asodte')) %>%
  filter(!is.na(asoswe))

# cl=parallel::makePSOCKcluster(4)
# test1=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff+fsca,data=alldata[1:500,],type.measure='mse',alpha=0.3)
# str(test1)
# test1$alpha=.3
# predict(test1,head(alldata))
# predict

# rpartfitting_phv(alldata)
# rpartfitting_phv=function(dF){
#   fit=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca,data=alldata,method='anova',cp=0.03,xval=10)
# }


cl <- parallel::makePSOCKcluster(8)
allphvmdls <-
  alldata %>%
  group_by(dte) %>%
  nest() %>%
  mutate(straps=map(data,crossv_mc,2)) %>%
  unnest(straps) %>%
  mutate(
    phv_obj_glmmdl=map(train,gnet_phv,cl),
    phv_r2_glmmdl=map2_dbl(phv_obj_glmmdl,test,myr2),
    phv_pctmae_glmmdl=map2_dbl(phv_obj_glmmdl,test,mypctmae),
    #
    phvfsca_obj_glmmdl=map(train,gnet_phvfsca,cl),
    phvfsca_r2_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,myr2),
    phvfsca_pctmae_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,mypctmae)
  )

# asoswe=asoswe %>% filter(asodte %in% unique(asodte)[5:9])
# allphvmdls <- allphvmdls %>% filter(.id=='01'|.id=='02')
## need to join aso swe data for asoswe models. doing phv separately before so it isn't fitting models on sets of the same data
resamp_data <-
  allphvmdls %>%
  group_by(dte,.id) %>%
  dplyr::select(dte,.id,train,test) %>%
  mutate(
    train=map(train,join_asodata,asoswe),
    test=map(test,join_asodata,asoswe)
  )


allasomdls <-
  resamp_data %>%
  group_by(dte,.id) %>%
  mutate(
    phvaso_obj_glmmdl=map(train,gnet_phvaso,cl),
    phvasofsca_obj_glmmdl=map(train,gnet_phvasofsca,cl)
  ) %>%
  unnest(phvaso_obj_glmmdl,phvasofsca_obj_glmmdl,.drop=FALSE)

if(!identical(allasomdls[[5]],allasomdls[[7]])) stop()

allasomdls <-
  allasomdls[!duplicated(as.list(allasomdls))] %>%
  mutate(
    phvaso_r2_glmmdl=map2_dbl(phvaso_obj_glmmdl,test,myr2),
    phvaso_pctmae_glmmdl=map2_dbl(phvaso_obj_glmmdl,test,mypctmae),
    phvasofsca_r2_glmmdl=map2_dbl(phvasofsca_obj_glmmdl,test,myr2),
    phvasofsca_pctmae_glmmdl=map2_dbl(phvasofsca_obj_glmmdl,test,mypctmae)
  )

parallel::stopCluster(cl)

saveRDS(allphvmdls,paste0('output/allphvmdls_phv-bootstrap_aso_',ires,'.rds'))
saveRDS(allasomdls,paste0('output/allasomdls_phv-bootstrap_aso_',ires,'.rds'))



# }    # <------ uncomment to run sensitivity

