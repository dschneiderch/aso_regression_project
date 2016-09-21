library(raster)
library(rgdal)
library(tidyverse)
library(cowplot)
library(viridis)
library(broom)
library(purrr)
library(modelr)
library(mgcv)
library(glmnetUtils)


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
phvdata=get_phvdata(geoarea2,ires,phv_stack) %>% filter(complete.cases(.))

sens=0
# for(sens in seq_len(7)){
# avgdate_rast=raster(paste0('data/snowoffdate/sensitivity/',sens,'/tuo_avgdate.tif'))
avgdate_rast=raster(paste0('data/snowoffdate/tuo_avgdate.tif'))
spdata=get_spdata(geoarea2,ires,avgdate_rast) %>% filter(complete.cases(.))

alldata <- inner_join(spdata,phvdata)

# cl=parallel::makePSOCKcluster(4)
# test1=cvAlpha.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+avgdate+fsca,data=alldata,type.measure='mse',outerParallel = cl)
# predict

gnetfitting=function(dF){
  cvfit=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight,data=dF,type.measure='mse',alpha=1,parallel = TRUE)
  return(cvfit)
}
gnetfitting_fsca=function(dF){
  cvfit=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca,data=dF,type.measure='mse',alpha=1,parallel = TRUE)
  return(cvfit)
}
gnetfitting_sp=function(dF){
  cvfit=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+avgdate,data=dF,type.measure='mse',alpha=1,parallel = TRUE)
  return(cvfit)
}
gnetfitting_spfsca=function(dF){
  cvfit=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+avgdate+fsca,data=dF,type.measure='mse',alpha=1,parallel = TRUE)
  return(cvfit)
}
#



cl=parallel::makePSOCKcluster(8)
rm(swe,avgdate)
allmdls <-
  alldata %>%
  group_by(dte) %>%
  do({
    datsplit=crossv_mc(.,10) %>%
      mutate(
        phv_obj_glmmdl=map(train, gnetfitting),
        phv_r2_glmmdl=map2_dbl(phv_obj_glmmdl,test,myr2),
        phv_pctmae_glmmdl=map2_dbl(phv_obj_glmmdl,test,mypctmae),
        #
        phvfsca_obj_glmmdl=map(train,gnetfitting_fsca),
        phvfsca_r2_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,myr2),
        phvfsca_pctmae_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,mypctmae),
        #
        phvaso_obj_glmmdl=map(train,gnetfitting_aso),
        phvaso_r2_glmmdl=map2_dbl(phvaso_obj_glmmdl,test,myr2),
        phvaso_pctmae_glmmdl=map2_dpl(phvaso_obj_glmmdl,test,mypctmae)
      )
  })
parallel::stopCluster(cl)

# gam models aren't working with resample obj so extract train/test data and fit separately
resamp_data <-
  allmdls %>%
  group_by(dte,.id) %>%
  do({
    data_frame(train=list(as.data.frame(.$train)),
               test=list(as.data.frame(.$test)))
  })

gammdls <-
  resamp_data %>%
  mutate(
    phv_obj_gammdl=map(train,~gam(swe~s(dem)+s(easting)+s(northing)+s(eastness)+s(northness)+s(zness)+s(vrm1cell)+s(stdslope1cell)+s(vegheight),select=TRUE,data=.)),
    phv_r2_gammdl=map2_dbl(phv_obj_gammdl,test,myr2),
    phv_pctmae_gammdl=map2_dbl(phv_obj_gammdl,test,mypctmae),
    #
    phvfsca_obj_gammdl=map(train,~gam(swe~s(dem)+s(easting)+s(northing)+s(eastness)+s(northness)+s(zness)+s(vrm1cell)+s(stdslope1cell)+s(vegheight)+s(fsca),select=TRUE,data=.)),
    phvfsca_r2_gammdl=map2_dbl(phvfsca_obj_gammdl,test,myr2),
    phvfsca_pctmae_gammdl=map2_dbl(phvfsca_obj_gammdl,test,mypctmae),
    #
    sp_obj_gammdl=map(train,~gam(swe~s(avgdate),data=.)),
    sp_r2_gammdl=map2_dbl(sp_obj_gammdl,test,myr2),
    sp_pctmae_gammdl=map2_dbl(sp_obj_gammdl,test,mypctmae),
    #
    spfsca_obj_gammdl=map(train, ~gam(swe~s(avgdate)+s(fsca),data=.)),
    spfsca_r2_gammdl=map2_dbl(spfsca_obj_gammdl,test,myr2),
    spfsca_pctmae_gammdl=map2_dbl(spfsca_obj_gammdl,test,mypctmae)
  )

moremdls=full_join(allmdls %>% ungroup %>% dplyr::select(-train,-test),
                   gammdls %>% ungroup )


# mdltypes=c("sp","spgam","spfsca","spfscagam","spxfsca","spxfscagam")
mdlstats <-
  moremdls %>%
  dplyr::select(-dplyr::contains('obj')) %>%
  gather(metricvar,metricval,-dte,-.id,-train,-test) %>%
  separate(metricvar,into=c("mdlpreds",'metricvar','mdltype')) %>%
  mutate(yr=substr(dte,1,4))

for(iyr in unique(mdlstats$yr)){
g1 <- mdlstats %>%
  filter_(~yr==iyr) %>%
  filter(metricvar=='r2') %>%
  # filter( !(metricvar=='r2' & metricval>1)) %>%
  # filter(metricvar!='r2') %>%
ggplot(data=.) +
  geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdltype,fill=mdlpreds))+
  scale_fill_brewer(palette='Set1',guide=F)+
  scale_colour_manual(values=c('black','grey80','grey30'),guide=F)+
# facet_grid(metricvar~.,scales='free')+
  labs(y=expression(r^2))+
  coord_cartesian(ylim=c(0,.75))+
  labs(x='')+
  facet_wrap(~yr,scales='free_x')+
  theme(plot.margin=unit(c(.05,.05,0,.05),"npc"),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank())#element_text(angle=45,hjust=1))


g2 <- mdlstats %>%
  filter_(~yr==iyr) %>%
  filter(metricvar=='pctmae') %>%
  # filter( !(metricvar=='r2' & metricval>1)) %>%
  # filter(metricvar!='r2') %>%
  ggplot(data=.) +
  geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdltype,fill=mdlpreds))+
  scale_fill_brewer(palette='Set1')+
  scale_colour_manual(values=c('black','grey80','grey30'))+
  labs(y='Mean Absolute Error [%]')+
  facet_wrap(~yr,scales='free_x')+
  # facet_grid(metricvar~.,scales='free')+
  coord_cartesian(ylim=c(20,100))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.margin=unit(c(0,.05,0.05,.05),"npc"),
        legend.direction='horizontal',
        legend.position='bottom')


pg=plot_grid(g1,g2,align='hv',nrow=2,rel_heights = c(.8,1))
pg
save_plot(plot=pg,filename=paste0('figs/mdlstats_sensitivity',sens,'-',iyr,'.pdf'),nrow=2,base_width=20,base_height=6)
}

mdlstats %>%
  # filter(metricvar!='r2') %>%
  group_by(metricvar,mdlpreds,mdltype) %>%
  summarise(avg=mean(metricval)) %>%
  spread(metricvar,avg) %>%
  readr::write_tsv(path=paste0('output/avgmdlstats_sensitivity',sens,'.txt'))

# }    # <------ uncomment to run sensitivity

