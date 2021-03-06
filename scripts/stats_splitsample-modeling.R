library(raster)
library(rgdal)
library(tidyverse)
library(cowplot)
library(viridis)
library(broom)


```{r params}
ires='500m'
basin='tuo'
pathin='output/splitsample-modeling/'
```

```{r load_data}
alldata=readRDS(paste0(pathin,'alldata_',ires,'.rds'))
phv=readRDS(paste0(pathin,'phvmdls_augment_',ires,'.rds'))
phvfsca=readRDS(paste0(pathin,'phvfscamdls_aug_',ires,'.rds'))
fns=dir(pathin,glob2rx(paste0('^phvasomdls_aug*',ires,'*.rds')),full.names=T)
phvaso =
  bind_rows(
    do.call(readRDS,as.list(fns))
    )
fns=dir(pathin,glob2rx(paste0('^phvasofscamdls_aug*',ires,'*.rds')),full.names=T)
phvasofsca =
  bind_rows(
    lapply(fns,readRDS)#do.call(readRDS,as.list(fns))
    )
allaug=phvmdls

extractAug <- function(allaug,mdlpreds,mdldtype){
  coln <- paste(mdlpreds,'aug',mdltype,sep='_')
  allaug %>%
    group_by(dte) %>%
    unnest_(coln) %>%
    # gather(var,val,swe:resid) %>%
    mutate(mdlpreds,mdltype)

}

mdlpreds='phvasofsca'
mdltype='glmmdl'
phv=extractAug(allaug,'phv','glmmdl')
phvasofsca=extractAug(allaug,'phvasofsca','glmmdl')

r2=function(dF){
  cor(dF$swe,dF$swehat)^2
}
pctmae=function(dF){
  mean(abs(dF$swehat-dF$swe),na.rm=T)/mean(dF$swe,na.rm=T)*100
}

predictions=bind_rows(phv,phvfsca)
mdlstats <-
  predictions %>%
  group_by(dte,mdlpreds,mdltype) %>%
  nest() %>%
  mutate(r2=map_dbl(data, r2),
        pctmae=map_dbl(data, pctmae),
        yr=substr(dte,1,4)
)

alldata
mdlstats

mdlstats2 <- mdlstats %>%
  gather(metricvar,metricval,r2,pctmae)
  # mdlstats %>%
  # filter(dte!=asodte)

mdlstats2 %>%
  group_by(mdlpreds,metricvar) %>%
  summarise(
    avg=median(metricval),
    N=n()
  )




mdlstats2 %>%
  # filter(metricvar!='r2') %>%
  group_by(metricvar,mdlpreds,mdltype) %>%
  summarise(avg=mean(metricval)) %>%
  spread(metricvar,avg) %>%
  readr::write_tsv(path=paste0('output/avgmdlstats_phv-bootstrap_aso_',ires,'.txt'))

g <-
  ggplot(mdlstats2)+
  geom_boxplot(aes(x=mdlpreds,y=metricval))+
  facet_wrap(~metricvar,scales='free')+
  theme(
    axis.text.x=element_text(angle=45,hjust=1)
  )
g
ggsave(plot=g,filename=paste0('figs/boxplot_phv-bootstrap_aso_',ires,'.pdf'),width=12,height=6)


wtest <-
  mdlstats %>%
  filter(dte!=asodte) %>%
  spread(mdlpreds,metricval)

wtest %>%
  filter(metricvar=='r2') %>%
  group_by(metricvar) %>%
  do(
    bind_rows(
      tidy(with(.,t.test(phvasofsca,phvaso,alternative='g',conf.int=T))) %>% mutate(test='t'),
      tidy(with(.,wilcox.test(phvasofsca,phvaso,alternative='g',conf.int=T))) %>% mutate(test='wilcox')
    ))
wtest %>%
  filter(metricvar=='r2') %>%
  group_by(metricvar) %>%
  do(
    bind_rows(
      tidy(with(.,t.test(phvfsca,phv,alternative='g',conf.int=T))) %>% mutate(test='t'),
      tidy(with(.,wilcox.test(phvfsca,phv,alternative='g',conf.int=T))) %>% mutate(test='wilcox')
    ))
wtest %>%
  filter(metricvar=='pctmae') %>%
  group_by(metricvar) %>%
  do(
    bind_rows(
      tidy(with(.,t.test(phvasofsca,phvaso,alternative='l',conf.int=T))) %>% mutate(test='t'),
      tidy(with(.,wilcox.test(phvasofsca,phvaso,alternative='l',conf.int=T))) %>% mutate(test='wilcox')
    ))
wtest %>%
  filter(metricvar=='pctmae') %>%
  group_by(metricvar) %>%
  do(
    bind_rows(
      tidy(with(.,t.test(phvfsca,phv,alternative='l',conf.int=T))) %>% mutate(test='t'),
      tidy(with(.,wilcox.test(phvfsca,phv,alternative='l',conf.int=T))) %>% mutate(test='wilcox')
    ))


iyr=2014
for(iyr in unique(mdlstats2$yr)){
  g1 <- mdlstats2 %>%
    filter_(~yr==iyr) %>%
    filter(metricvar=='r2') %>%
    # filter( !(metricvar=='r2' & metricval>1)) %>%
    # filter(metricvar!='r2') %>%
    ggplot(data=.) +
    geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdlpreds,fill=mdlpreds))+
    scale_fill_brewer(palette='Set1',guide=F)+
    scale_colour_manual(values=c('black','grey80','grey30'),guide=F)+
    # facet_grid(metricvar~.,scales='free')+
    labs(y=expression(r^2))+
    coord_cartesian(ylim=c(0,1))+
    labs(x='')+
    facet_wrap(~yr,scales='free_x')+
    theme(plot.margin=unit(c(.05,.05,0,.05),"npc"),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())#element_text(angle=45,hjust=1))


  g2 <- mdlstats2 %>%
    filter_(~yr==iyr) %>%
    filter(metricvar=='pctmae') %>%
    # filter( !(metricvar=='r2' & metricval>1)) %>%
    # filter(metricvar!='r2') %>%
    ggplot(data=.) +
    geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdlpreds,fill=mdlpreds))+
    scale_fill_brewer(palette='Set1')+
    scale_colour_manual(values=c('black','grey80','grey30'))+
    labs(y='Mean Absolute Error [%]')+
    facet_wrap(~yr,scales='free_x')+
    # facet_grid(metricvar~.,scales='free')+
    coord_cartesian(ylim=c(0,100))+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          plot.margin=unit(c(0,.05,0.05,.05),"npc"),
          legend.direction='horizontal',
          legend.position='bottom')


  pg=plot_grid(g1,g2,align='hv',nrow=2,rel_heights = c(.8,1))
  pg
  save_plot(plot=pg,filename=paste0('figs/mdlstats_phv-bootstrap_aso-',iyr,'_',ires,'.pdf'),nrow=2,base_width=20,base_height=6)
}
