# Data Collation ----
get_phvdata=function(geoarea2,ires,phv_stack,fscasource){
  if(fscasource=='modscag') {
    fscafn <- dir('data/fsca/modscag/utm',pattern=glob2rx(paste0('fsca_',ires,'_',geoarea2,'_*.tif$')),full.names=T)
    fscadtes <- strftime(strptime(sapply(strsplit(basename(fscafn),'[.\\_]'),'[',4),'%Y%j'),'%Y%m%d')
  }
  if(fscasource=='aso') {
    fscafn <- dir(paste0('data/fsca/aso/utm/',ires),pattern=glob2rx(paste0('fsca_',ires,'_',geoarea2,'_*.tif$')),full.names=T)
    fscadtes <- sapply(strsplit(basename(fscafn),'[.\\_]'),'[',4)
  }
  fscaind <- which(fscadtes %in% asoswedates)
  fsca_stack <- stack(fscafn[fscaind])
  #replace names in stack with %Y%m%d (in case they aren't like from modscag)
  names(fsca_stack)  <- sapply(seq_along(names(fsca_stack)),function(x) gsub('[0-9]+$',names(fsca_stack)[x],fscadtes[x]))

  a1=as.data.frame(aso_stack,xy=T) %>%
    tbl_df %>%
    gather(dte,swe,-x,-y) %>%
    separate(dte,into=c('basin','res','dte')) %>%
    filter(!is.na(swe))

  p1=as.data.frame(phv_stack,xy=T) %>%
    tbl_df %>%
    gather(phvvar,phvval,-x,-y) %>%
    mutate(phvvar=sapply(strsplit(phvvar,'_',fixed=T),'[',2),
           phvvar=gsub('[0-9]+m$','',phvvar)) %>%
    spread(phvvar,phvval)

  asoswe <-
    inner_join(a1,p1)

  fscadf=as.data.frame(fsca_stack,xy=T) %>%
    gather(dte,fsca,-x,-y) %>%
    separate(dte,into=c('var','ires','basin','dte')) %>%
    dplyr::select(-var,-ires) %>%
    tbl_df

  dat=asoswe %>%
    inner_join(fscadf) %>%
    filter(fsca>0)

  return(dat)
}

join_asodata=function(dF,asoswe2){
  dfout <-
    as.data.frame(dF) %>%
    inner_join(asoswe2,by=c('x','y','basin','res'))
  dfout
}

# error functions ----

myr2=function(model,data){
  yvar=all.vars(formula(model)[[2]])
  stats::cor(predict(model, data), as.data.frame(data)[[yvar]], use='complete.obs')^2
}


mypctmae=function(model,data) {
  yvar=all.vars(formula(model)[[2]])
  mae(model,data)/mean(as.data.frame(data)[[yvar]],na.rm=T)*100
}


# Elastic-Net Regularized Generalized Linear Models ----

pickAlpha=function(dF,myformula,cl=cl){
  cvalpha=cvAlpha.glmnet(myformula,data=dF,type.measure='mse',outerParallel=cl)

  maxlambda=100
  alpha=cvalpha$alpha
  ialpha=1
  cvdf=data_frame(alpha) %>%
    bind_cols(data.frame(matrix(NA,ncol=maxlambda,nrow=length(alpha)))) %>%
    bind_cols(as.data.frame(matrix(NA,ncol=maxlambda,nrow=length(alpha))))#
  for(ialpha in seq_along(alpha)){
    cvlambdas=cvalpha$modlist[[ialpha]]$lambda
    numlambda=length(cvlambdas)
    cvdf[ialpha,2:(numlambda+1)] <- cvlambdas
    cvdf[ialpha,maxlambda+1+(1:numlambda)] <- cvalpha$modlist[[ialpha]]$cvm
  }
  cvresults <-
    cvdf %>%
    mutate(alpha=alpha) %>%
    gather(lambdaid,lambdaval,num_range('X',1:100)) %>%
    gather(mseid,mseval,num_range('V',1:100))

  med_mse <-
    cvresults %>%
    group_by(alpha) %>%
    summarise(medmse=median(mseval,na.rm=T))

  med_se <- med_mse %>%
    summarise(
      se=sd(medmse,na.rm=T)/sqrt(n())
    ) %>% as.numeric

  alphaind=which.min(abs(med_mse$medmse-(min(med_mse$medmse,na.rm=T)+med_se)))[1]
  bestalpha=alpha[alphaind]

  return(bestalpha)
}

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

gnet_phvaso=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+asoswe)
  wrap_cvnet=function(dF,bestalpha,myformula) cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)

 dF %>%
    group_by(asodte) %>%
    nest %>%
    mutate(
      bestalpha=map_dbl(data,pickAlpha,myformula,cl),
      phvaso_obj_glmmdl=map2(data,bestalpha,wrap_cvnet,myformula)
    ) %>%
   dplyr::select(asodte,phvaso_obj_glmmdl)

}

gnet_phvasofsca=function(dF,cl=cl){
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+asoswe+fsca)
  wrap_cvnet=function(dF,bestalpha,myformula) cv.glmnet(myformula,data=dF,type.measure='mse',alpha=bestalpha,parallel = TRUE)
  dF %>%
    group_by(asodte) %>%
    nest %>%
    mutate(
      bestalpha=map_dbl(data,pickAlpha,myformula,cl),
      phvasofsca_obj_glmmdl=map2(data,bestalpha,wrap_cvnet,myformula)
    ) %>%
    dplyr::select(asodte,phvasofsca_obj_glmmdl)
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



# regression tree ----
pruneRpart=function(rp){
  cptable=rp$cptable
  minxerror.ind=which.min(cptable[,"xerror"])
  minerr=cptable[minxerror.ind,'xerror']
  minstd=cptable[minxerror.ind,'xstd']
  rpind=which(cptable[,'xerror'] < (minerr+minstd))[1]
  cp=cptable[rpind,'CP']
  prune(rp,cp=cp)
}
rpart_phvfsca=function(dF){
  rp=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca,data=dF,method='anova',cp=0.01,xval=10)
  pruneRpart(rp)
}
rpart_phv=function(dF){
  rp=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight,data=dF,method='anova',cp=0.01,xval=10)
  pruneRpart(rp)
}
rpart_phvsp=function(dF){
  rp=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff,data=dF,method='anova',cp=0.01,xval=10)
  pruneRpart(rp)
}
rpart_phvspfsca=function(dF){
  rp=rpart(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+snowoff+fsca,data=dF,method='anova',cp=0.01,xval=10)
  pruneRpart(rp)
}
