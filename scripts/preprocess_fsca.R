library(raster)
geoarea2='tuo'
ires='500m'

## - modscag-historic fsca - always 500m ----
fn_fsca=dir('data/fsca/modscag/raw/sinu',pattern=glob2rx(paste0('^',geoarea2,'*.tif$')),full.names=TRUE)
dte_str=sapply(strsplit(x=basename(fn_fsca),split='[_//.]'),'[',2)
datesOI <- which(substr(dte_str,5,7)<245)#245 is through end of august
dtes=as.Date(strptime(dte_str,'%Y%j'))[datesOI]
yrs=strftime(dtes,'%Y')[datesOI]

fsca_stack <- stack(fn_fsca[datesOI])
fs_names=names(fsca_stack)
# fsca_stack[fsca_stack>100]=NA
beginCluster()
fsca_stack <- clusterR(fsca_stack,fun=function(x) {x[x>100] <- NA; return(x)})
endCluster()

basinmask=raster(paste0('data/gis/tuo_basinmask_500m.tif'))

beginCluster()
fsca_stack <- projectRaster(fsca_stack,basinmask,method='bilinear')
endCluster()
names(fsca_stack) <- fs_names
writeRaster(fsca_stack,filename='data/fsca/modscag/utm/fsca_500m.tif',bylayer=TRUE,suffix='names',overwrite=T,NAflag=-99)


## - aso 51m, 99m and 501m fsca created in snowdepth-predict_project from 3m snowdepth ----
# resample to 500 m
oldres=99
newres=100
dir.create(paste0('data/fsca/aso/utm/',newres,'m'),recursive = TRUE)
basinmask=raster(paste0('data/gis/tuo_basinmask_',newres,'m.tif'))
asofsca_fn=dir(paste0('data/fsca/aso/utm/',oldres,'m'),pattern=glob2rx('fsca*.tif$'),full.names = T)
asofsca=stack(asofsca_fn)
asofsca_new <- resample(asofsca,basinmask,method='ngb')
names(asofsca_new) <- gsub('fsca','tuo',names(asofsca_new))
writeRaster(asofsca_new,filename=paste0('data/fsca/aso/utm/',newres,'m/fsca_',newres,'m.tif'),bylayer=T,suffix='names',NAflag=-99,overwrite=T)
