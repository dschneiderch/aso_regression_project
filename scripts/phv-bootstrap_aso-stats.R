ires='500m'
basin='tuo'

allphvmdls=readRDS(paste0('output/allphvmdls_phv-bootstrap_aso_',ires,'.rds'))
allasomdls=readRDS(paste0('output/allasomdls_phv-bootstrap_aso_',ires,'.rds'))


totalmdls <- full_join(allphvmdls %>% dplyr::select(-train,-test),allasomdls %>% dplyr::select(-train,-test))
# saveRDS(totalmdls,'output/phv-bootstrap_aso.rds')

mdlstats <-
  totalmdls %>%
  dplyr::select(dte,asodte,.id,contains('r2'),contains('pctmae')) %>%
  gather(metricvar, metricval,-dte,-asodte,-.id) %>%
  separate(metricvar,into=c('mdlpreds','metricvar','mdltype')) %>%
  mutate(yr=substr(dte,1,4))

mdlstats2 <-
  mdlstats %>%
  filter(dte!=asodte)

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
    geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdltype,fill=mdlpreds))+
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
    geom_boxplot(aes(x=as.factor(dte),y=metricval,colour=mdltype,fill=mdlpreds))+
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
