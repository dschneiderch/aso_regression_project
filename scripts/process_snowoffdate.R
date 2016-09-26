library(raster)

outpath='data/snowoffdate'
if(!dir.exists(outpath)) dir.create(outpath,recursive=TRUE)

ires='500m'
geoarea2='tuo'
basinmask=raster(paste0('data/gis/',geoarea2,'_basinmask_',ires,'.tif'))

fn_fsca=dir('data/fsca/modscag/utm/',pattern=glob2rx(paste0('fsca*',geoarea2,'*.tif')),full.names=TRUE)
dte_str=sapply(strsplit(x=basename(fn_fsca),split='[_//.]'),'[',4)
dtes=as.Date(strptime(dte_str,'%Y%j'))
yrs=strftime(dtes,'%Y')
yrind=which(yrs<2016 & yrs>2012)
fn_fsca <- fn_fsca[yrind]
dte_str <- dte_str[yrind]
dtes <- dtes[yrind]
yrs <- yrs[yrind]


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
writeRaster(date_stack,file.path(outpath,paste0(geoarea2,'_snowoff.tif')),NAflag=-99,bylayer=TRUE,suffix='names',overwrite=T)


# cvdate <- calc(date_stack,fun = function(x) sd(x)/mean(x))*basinmask
# writeRaster(cvdate,filename=paste0('data/snowoffdate/',geoarea2,'_cvdate.tif'),NAflag=-99,overwrite=T)
# avgdate <- mean(date_stack)*basinmask
# writeRaster(avgdate,filename=paste0('data/snowoffdate/',geoarea2,'_avgdate.tif'),NAflag=-99,overwrite=T)
