# Data Collation ----
get_asoswe=function(geoarea2,ires){
  asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)
  aso_stack=stack(asofn)
  # asoswedates=sapply(strsplit(names(aso_stack),'[.\\_]'),'[',3)

  # plot(aso_stack[[22:33]],colNA='black')
  # aso_stack[aso_stack<=0] <- NA
  # plot(aso_stack[[c(1,2,8,9)]],colNA='black')
  # cellStats(aso_stack[[c(1,2,8,9)]],'range')

  asoswe=as.data.frame(aso_stack,xy=T) %>%
    tbl_df %>%
    mutate(x=as.character(x),
           y=as.character(y)) %>%
    gather(dte,swe,-x,-y) %>%
    separate(dte,into=c('basin','res','dte')) %>%
    filter(!is.na(swe),swe>=0)#need to remvoe some random -9999

  return(asoswe)
}


get_phvdata=function(geoarea2,ires,phv_stack,a1,asoswedates1,fscasource){

  p1=as.data.frame(phv_stack,xy=T) %>%
    tbl_df %>%
    mutate(x=as.character(x),
           y=as.character(y)) %>%
    gather(phvvar,phvval,-x,-y) %>%
    mutate(phvvar=sapply(strsplit(phvvar,'_',fixed=T),'[',2),
           phvvar=gsub('[0-9]+m$','',phvvar)) %>%
    spread(phvvar,phvval)

  phvswe <-
    inner_join(a1,p1)

  if(fscasource=='modscag') {
    fscafn <- dir('data/fsca/modscag/utm',pattern=glob2rx(paste0('fsca_',ires,'_',geoarea2,'_*.tif$')),full.names=T)
    fscadtes <- strftime(strptime(sapply(strsplit(basename(fscafn),'[.\\_]'),'[',4),'%Y%j'),'%Y%m%d')
  }
  if(fscasource=='aso') {
    fscafn <- dir(paste0('data/fsca/aso/utm/',ires),pattern=glob2rx(paste0('fsca_',ires,'_',geoarea2,'_*.tif$')),full.names=T)
    fscadtes <- sapply(strsplit(basename(fscafn),'[.\\_]'),'[',4)
  }
  fscaind <- which(fscadtes %in% asoswedates1)
  fsca_stack <- stack(fscafn[fscaind])
  #replace names in stack with %Y%m%d (in case they aren't like from modscag)
  names(fsca_stack)  <- sapply(seq_along(names(fsca_stack)),function(x) gsub('[0-9]+$',names(fsca_stack)[x],fscadtes[x]))

  fscadf=as.data.frame(fsca_stack,xy=T) %>%
    mutate(x=as.character(x),
           y=as.character(y)) %>%
    gather(dte,fsca,-x,-y) %>%
    separate(dte,into=c('var','res','basin','dte')) %>%
    dplyr::select(-var) %>%
    tbl_df

  dat=phvswe %>%
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

add_asodata=function(data1,bestasodates,asoswe2){

  newdata=left_join(data1,bestasodates,by=c('dte'))
  newdata=newdata %>% #use left join so we don't lose data. need rows to stay the same as in alldata
    left_join(asoswe2,by=c('x','y','basin','res','bestasodte'='asodte')) %>%
    rename(asodte=bestasodte) %>%
    dplyr::select(-dte)

  return(newdata)
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

r2=function(yobs,yhat){
  cor(yobs,yhat,use='complete.obs')^2
}

rootmse = function(yobs,yhat){
  sqrt(mean((yobs-yhat)^2,na.rm=T))
}


# Elastic-Net Regularized Generalized Linear Models ----

augmentEnet=function(model,data,asodte=NULL,asoswe2=NULL){

  realdata=as.data.frame(data)
  # print(str(realdata))
  if(!is.null(asodte)) realdata=join_asodata(realdata,asoswe2 %>% filter_(~asodte==asodte))
  x=realdata[['x']]
  y=realdata[['y']]

  yvar=all.vars(formula(model)[[2]])
  yobs=realdata[[yvar]]

  yhat=as.numeric(predict(model,newdata=realdata))
  data_frame(x,y,yobs,yhat) %>% setNames(c('x','y',yvar,paste0(yvar,'hat')))
}


pickAlpha=function(dF,myformula,nfolds,cl=cl){
  cvalpha=cvAlpha.glmnet(myformula,data=dF,nfolds=nfolds,type.measure='mse',outerParallel=cl)

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
    summarise(
      medmse=median(mseval,na.rm=T)
    )

  med_se <- med_mse %>%
    summarise(
      se=sd(medmse,na.rm=T)/sqrt(n())
    ) %>% as.numeric

  alphaind=which.min(abs(med_mse$medmse-(min(med_mse$medmse,na.rm=T)+med_se)))[1]
  bestalpha=alpha[alphaind]
  # bestalpha=1
  return(bestalpha)
}

gnet_phv=function(dF,cl=NULL){
  if(is.null(cl)) yesParallel=FALSE else yesParallel=TRUE
  # print(unique(dF$dte2))
  # print(str(dF))
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight)
  nfolds=floor(nrow(dF)/30)
  # print(nfolds)
  if(nfolds<4) nfolds=4
  if(nfolds>10) nfolds=10
  # print(paste('update:',nfolds))
  bestalpha <- pickAlpha(dF,myformula,nfolds,cl)

  cvfit=cv.glmnet(myformula,data=dF,nfolds=nfolds,type.measure='mse',alpha=bestalpha,parallel = yesParallel)
  cvfit$alpha=bestalpha
  return(cvfit)
}

gnet_phvfsca=function(dF,cl=NULL){
  if(is.null(cl)) yesParallel=FALSE else yesParallel=TRUE
  nfolds=floor(nrow(dF)/30)
  if(nfolds<4) nfolds=4
  if(nfolds>10) nfolds=10
  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+fsca)
  bestalpha <- pickAlpha(dF,myformula,nfolds,cl)
  cvfit=cv.glmnet(myformula,data=dF,nfolds=nfolds,type.measure='mse',alpha=bestalpha,parallel = yesParallel)
  cvfit$alpha=bestalpha
  return(cvfit)
}


gnet_phvaso=function(dF,join_asodata=NULL,asoswe2=NULL,cl=NULL,bestasodte=NULL){
  if(!is.null(bestasodte)) {asoswe2=filter(asoswe2,asodte %in% bestasodte)}
  if(!is.null(asoswe2)) {dF=join_asodata(dF,asoswe2)}
  if(is.null(cl)) yesParallel=FALSE else yesParallel=TRUE
  dFsize <- dF %>% split(.$asodte) %>% map_dbl(nrow) %>% min
  nfolds=floor(dFsize/30)
  if(nfolds<4) nfolds=4
  if(nfolds>10) nfolds=10

  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+asoswe)
  wrap_cvnet=function(dF,bestalpha,myformula,nfolds) cv.glmnet(myformula,data=dF,nfolds=nfolds,type.measure='mse',alpha=bestalpha,parallel = yesParallel)

  dFout=dF %>%
    group_by(asodte) %>%
    nest %>%
    mutate(
      bestalpha=map_dbl(data,pickAlpha,myformula,nfolds,cl),
      phvaso_obj_glmmdl=map2(data,bestalpha,wrap_cvnet,myformula,nfolds)
    ) %>%
    dplyr::select(asodte,phvaso_obj_glmmdl)

  return(dFout)
}

gnet_phvasofsca=function(dF,join_asodata=NULL,asoswe2=NULL,cl=NULL,bestasodte=NULL){
  if(!is.null(bestasodte)) asoswe2=filter(asoswe2,asodte %in% bestasodte)
  if(!is.null(asoswe2)) dF=join_asodata(dF,asoswe2)
  if(is.null(cl)) yesParallel=FALSE else yesParallel=TRUE
  dFsize <- dF %>% split(.$asodte) %>% map_dbl(nrow) %>% min
  nfolds=floor(dFsize/30)
  if(nfolds<4) nfolds=4
  if(nfolds>10) nfolds=10

  myformula=as.formula(swe~dem+easting+northing+eastness+northness+zness+vrm1cell+stdslope1cell+vegheight+asoswe+fsca)
  wrap_cvnet=function(dF,bestalpha,myformula,nfolds) cv.glmnet(myformula,data=dF,nfolds=nfolds,type.measure='mse',alpha=bestalpha,parallel = TRUE)

  dFout=dF %>%
    group_by(asodte) %>%
    nest %>%
    mutate(
      bestalpha=map_dbl(data,pickAlpha,myformula,nfolds,cl),
      phvasofsca_obj_glmmdl=map2(data,bestalpha,wrap_cvnet,myformula,nfolds)
    ) %>%
    dplyr::select(asodte,phvasofsca_obj_glmmdl)

  return(dFout)
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
