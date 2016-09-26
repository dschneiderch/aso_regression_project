library(raster)
library(rgdal)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(broom)
library(purrr)
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


# cl=parallel::makePSOCKcluster(4)
# test1=cv.glmnet(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff+fsca,data=alldata[1:500,],type.measure='mse',alpha=0.3)
# str(test1)
# test1$alpha=.3
# predict(test1,head(alldata))
# predict

gnet_phv=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight)
  bestalpha <- pickAlpha(dF,myformula,cl)
  cvfit=cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)
  cvfit$alpha=bestalpha
  return(cvfit)
}
gnet_phvfsca=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca)
  bestalpha <- pickAlpha(dF,myformula,cl)
  cvfit=cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)
  cvfit$alpha=bestalpha
  return(cvfit)
}
gnet_phvsp=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff)
  bestalpha <- pickAlpha(dF,myformula,cl)
  cvfit=cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)
  cvfit$alpha=bestalpha
  return(cvfit)
}
gnet_phvfscasp=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff+fsca)
  bestalpha <- pickAlpha(dF,myformula,cl)
  cvfit=cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)
  cvfit$alpha=bestalpha
  return(cvfit)
}
#
# rpartfitting_phv(alldata)
# rpartfitting_phv=function(dF){
#   fit=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca,data=alldata,method='anova',cp=0.03,xval=10)
# }

cl <- parallel::makePSOCKcluster(6)
print(cl)

allmdls <-
  alldata %>%
  group_by(dte) %>%
  do({
    datsplit=crossv_mc(.,10) %>%
      mutate(
        phv_obj_glmmdl=map(train, gnet_phv,cl),
        phv_r2_glmmdl=map2_dbl(phv_obj_glmmdl,test,myr2),
        phv_pctmae_glmmdl=map2_dbl(phv_obj_glmmdl,test,mypctmae),
        #
        phvfsca_obj_glmmdl=map(train,gnet_phvfsca,cl),
        phvfsca_r2_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,myr2),
        phvfsca_pctmae_glmmdl=map2_dbl(phvfsca_obj_glmmdl,test,mypctmae),
        #
        phvsp_obj_glmmdl=map(train,gnet_phvsp,cl),
        phvsp_r2_glmmdl=map2_dbl(phvsp_obj_glmmdl,test,myr2),
        phvsp_pctmae_glmmdl=map2_dbl(phvsp_obj_glmmdl,test,mypctmae),
        #
        phvspfsca_obj_glmmdl=map(train,gnet_phvfscasp,cl),
        phvspfsca_r2_glmmdl=map2_dbl(phvspfsca_obj_glmmdl,test,myr2),
        phvspfsca_pctmae_glmmdl=map2_dbl(phvspfsca_obj_glmmdl,test,mypctmae),
        #
        sp_obj_glmmdl=map(train, ~glm(swe~snowoff,data=.)),
        sp_r2_glmmdl=map2_dbl(sp_obj_glmmdl,test,myr2),
        sp_pctmae_glmmdl=map2_dbl(sp_obj_glmmdl,test,mypctmae),
        #
        spfsca_obj_glmmdl=map(train, ~glm(swe~snowoff+fsca,data=.)),
        spfsca_r2_glmmdl=map2_dbl(spfsca_obj_glmmdl,test,myr2),
        spfsca_pctmae_glmmdl=map2_dbl(spfsca_obj_glmmdl,test,mypctmae),
        #
        phv_obj_treemdl=map(train, ~rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight, data=.,method='anova',cp=0.03,xval=10)),
        phv_r2_treemdl=map2_dbl(phv_obj_treemdl,test,myr2),
        phv_pctmae_treemdl=map2_dbl(phv_obj_treemdl,test,mypctmae),
        #
        phvfsca_obj_treemdl=map(train, ~rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca,data=.,method='anova',cp=0.03,xval=10)),
        phvfsca_r2_treemdl=map2_dbl(phvfsca_obj_treemdl,test,myr2),
        phvfsca_pctmae_treemdl=map2_dbl(phvfsca_obj_treemdl,test,mypctmae),
        #
        phvspfsca_obj_treemdl=map(train, ~rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff+fsca,data=.,method='anova',cp=0.03,xval=10)),
        phvspfsca_r2_treemdl=map2_dbl(phvspfsca_obj_treemdl,test,myr2),
        phvspfsca_pctmae_treemdl=map2_dbl(phvspfsca_obj_treemdl,test,mypctmae)
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
    sp_obj_gammdl=map(train,~gam(swe~s(snowoff),data=.)),
    sp_r2_gammdl=map2_dbl(sp_obj_gammdl,test,myr2),
    sp_pctmae_gammdl=map2_dbl(sp_obj_gammdl,test,mypctmae),
    #
    spfsca_obj_gammdl=map(train, ~gam(swe~s(snowoff)+s(fsca),data=.)),
    spfsca_r2_gammdl=map2_dbl(spfsca_obj_gammdl,test,myr2),
    spfsca_pctmae_gammdl=map2_dbl(spfsca_obj_gammdl,test,mypctmae)
  )

moremdls=full_join(allmdls %>% ungroup %>% dplyr::select(-train,-test),
                   gammdls %>% ungroup )

saveRDS(moremdls,'output/modelfitting_moremdls.rds')

# mdltypes=c("sp","spgam","spfsca","spfscagam","spxfsca","spxfscagam")
mdlstats <-
  moremdls %>%
  dplyr::select(-contains('obj')) %>%
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

