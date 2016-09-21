library(raster)

outpath='data/snowoffdate/sensitivity/'
ires='500m'
geoarea2='tuo'
basinmask=raster(paste0('data/gis/',geoarea2,'_basinmask_',ires,'.tif'))

fn_fsca_all=dir('data/fsca/modscag/utm/',pattern=glob2rx(paste0('fsca*',geoarea2,'*.tif')),full.names=TRUE)
dte_str_all=sapply(strsplit(x=basename(fn_fsca_all),split='[_//.]'),'[',4)
dtes_all=as.Date(strptime(dte_str_all,'%Y%j'))
yrs_all=strftime(dtes,'%Y')
yrind=which(yrs_all<2013)
fn_fsca_all <- fn_fsca_all[yrind]
dte_str_all <- dte_str_all[yrind]
dtes_all <- dtes_all[yrind]
yrs_all <- yrs_all[yrind]


dayi=1
for(dayi in seq_len(7)){
  outpath2=file.path(outpath,dayi)
  if(!dir.exists(outpath2)) dir.create(outpath2,recursive=TRUE)

sensitivity_days=rep(seq(dayi,245,7),length(unique(yrs_all)))
sensfn_ind=which(substr(dte_str_all,5,7) %in% sensitivity_days)

dte_str=dte_str_all[sensfn_ind]
yrs=yrs_all[sensfn_ind]
fn_fsca=fn_fsca_all[sensfn_ind]
fsca_stack <- stack(fn_fsca)
fs_names=names(fsca_stack)

# snow off date each year ----
statfun=function(x,na.rm=TRUE){
  which(x==0)[1]
}
beginCluster()
date_ind <- clusterR(fsca_stack,stackApply,args=list(indices=as.integer(yrs),fun=statfun))#stackApply returns index layers in order of the first appearance in the index value.
endCluster()
names(date_ind) <- unique(yrs)

iyr=2000
date_stack <- date_ind
for(iyr in yrs){
  yrlayers <- dte_str[grep(iyr,yrs)]
  ilayer <- grep(iyr,names(date_ind))
  sdd <- substr(yrlayers[getValues(date_ind[[ilayer]])],5,7)
  date_stack[[ilayer]]=setValues(date_ind[[ilayer]],values=as.integer(sdd))
}
writeRaster(date_stack,file.path(outpath2,paste0(geoarea2,'_snowoff.tif')),NAflag=-99,bylayer=TRUE,suffix='names',overwrite=T)


cvdate <- calc(date_stack,fun = function(x) sd(x)/mean(x))*basinmask
writeRaster(cvdate,filename=paste0(outpath2,'/',geoarea2,'_cvdate.tif'),NAflag=-99,overwrite=T)
avgdate <- mean(date_stack)*basinmask
writeRaster(avgdate,filename=paste0(outpath2,'/',geoarea2,'_avgdate.tif'),NAflag=-99,overwrite=T)

}
