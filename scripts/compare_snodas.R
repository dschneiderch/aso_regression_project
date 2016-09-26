library(raster)
library(rgdal)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(broom)

inpath='data/swe'
dte=20150430
geoarea2='tuo'
ires='500m'
asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)

aso_stack=stack(asofn)
avgdate_rast=raster('data/snowoffdate/tuo_avgdate.tif')

asoswe <- as.data.frame(aso_stack,xy=T) %>%
  tbl_df %>%
  bind_cols(as.data.frame(avgdate_rast)) %>%
  rename(avgdate=tuo_avgdate) %>%
  gather(dte,swe,-x,-y,-avgdate) %>%
  separate(dte,into=c('basin','res','dte')) %>%
  filter(!is.na(avgdate),swe>0)


# cvdte=as.data.frame(mask(cvdate,aso))
ggplot(asoswe,aes(x=avgdate,y=swe))+
  geom_point(aes(colour=swe))+
  # geom_vline(aes(xintercept=130))+
  # geom_hline(aes(yintercept=0.061))+
  facet_wrap(~dte,scales='free')

##
fscafn <- dir('data/fsca/modscag/utm',pattern=glob2rx(paste0('fsca_',ires,'_',geoarea2,'_*.tif$')),full.names=T)
fscadtes <- strftime(strptime(sapply(strsplit(basename(fscafn),'[.\\_]'),'[',4),'%Y%j'),'%Y%m%d')
fscaind <- which(fscadtes %in% unique(asoswe$dte))

fscadf=as.data.frame(stack(fscafn[fscaind]),xy=T) %>%
  gather(dte,fsca,-x,-y) %>%
  separate(dte,into=c('var','ires','basin','dte')) %>%
  select(-var) %>%
  mutate(dte=strftime(strptime(dte,'%Y%j'),'%Y%m%d')) %>%
  tbl_df

dat=asoswe %>%
  inner_join(fscadf) %>%
  filter(fsca>0)

dat %>%
  arrange(fsca) %>%
  ggplot(.,aes(x=avgdate,y=swe))+
  geom_point(aes(color=fsca))+
  facet_wrap(~dte,scales='free_y')+
  scale_color_viridis(guide='colorbar')



