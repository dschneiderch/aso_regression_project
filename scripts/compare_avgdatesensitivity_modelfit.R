## compare sensitivity avgdate regressions ----
library(tidyverse)
library(cowplot)

avgmdlstats_fn <- dir('output',pattern='avgmdlstats',full.names=T)
avgmdlstats <- tibble()
for(f in avgmdlstats_fn){
  ams <- readr::read_tsv(f) %>%
    mutate(id=sapply(strsplit(basename(f),'[_\\.]'),'[',2))
  avgmdlstats <- bind_rows(avgmdlstats,ams)
}


r2sens <- avgmdlstats %>%
  select(mdlpreds,mdltype,r2,id) %>%
  spread(id,r2) %>%
  gather(sensitivity,val,sensitivity1:sensitivity7) %>%
  rename(baseline=sensitivity0) %>%
  mutate(diff=baseline-val,
         metricvar='r2')

pctmaesens <- avgmdlstats %>%
  select(mdlpreds,mdltype,pctmae,id) %>%
  spread(id,pctmae) %>%
  gather(sensitivity,val,sensitivity1:sensitivity7) %>%
  rename(baseline=sensitivity0) %>%
  mutate(diff=baseline-val,
         metricvar='pctmae')

avgmdlstats_sens <- bind_rows(r2sens,pctmaesens)
avgmdlstats_sens %>%
  group_by(mdltype,mdlpreds,metricvar) %>%
  summarise(
    avgdiff=mean(diff)) %>%
  spread(metricvar,avgdiff)

avgmdlstats_sens %>%
  group_by(mdltype,mdlpreds,metricvar) %>%
  summarise(
    summdiff=min(diff)) %>%
  spread(metricvar,summdiff)

avgmdlstats_sens %>%
  group_by(mdltype,mdlpreds,metricvar) %>%
  summarise(
    summdiff=max(diff)) %>%
  spread(metricvar,summdiff)

plotsens <-
  avgmdlstats_sens %>%
  gather(modelspec,metricval,baseline,val,diff) %>%
  mutate(modelspec=factor(modelspec,levels=c('baseline','val','diff')))

for(im in unique(plotsens$mdltype)){
  gsens <-
    filter_(plotsens,~mdltype==im) %>%
    ggplot(.)+
    geom_bar(aes(x=mdlpreds,y=metricval,fill=modelspec),stat='identity',position='dodge')+
    scale_fill_manual(values=c('blue','forestgreen','black'))+
    facet_grid(metricvar~sensitivity,scales='free')+
    theme(axis.text.x=element_text(angle=45,hjust=1))

  ggsave(gsens,filename=paste0('figs/mdlstat_sensitivity_differences_',im,'.pdf'),width=12,height=4)
}
