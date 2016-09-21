library(raster)

## - modscag-historic fsca
geoarea2='tuo'
ires='500m'

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


## symlinks of modscag images for aso flight dates
asofn=dir('data/swe',pattern=glob2rx(paste0(geoarea2,'_',ires,'_*.tif$')),full.names=T)
asoswedates=sapply(strsplit(names(aso_stack),'[.\\_]'),'[',3)
asd_yj=strftime(strptime(asoswedates,'%Y%m%d'),'%Y%j')
asodteind=which(dte_str %in% asd_yj)
for(f in fn_fsca[asodteind]){
  dte=strftime(strptime(unlist(strsplit(basename(f),'[_\\.]'))[4],'%Y%j'),'%Y%m%d')
  newf=paste0('fsca_500m_tuo_',dte,'.tif')
  file.symlink(f,paste0('data/fsca/modscag/',newf))
}


## - aso 501m fsca created in snowdepth-predict_project from 3m snowdepth
# resample to 500 m
dir.create('data/fsca/aso/500m',recursive = TRUE)
basinmask=raster('data/gis/tuo_basinmask_500m.tif')
asofsca_fn=dir('data/fsca/aso/501m',pattern=glob2rx('fsca*.tif$'),full.names = T)
asofsca=stack(asofsca_fn)
asofsca500 <- resample(asofsca,basinmask,method='ngb')
writeRaster(asofsca500,filename='data/fsca/aso/500m/tuo.tif',bylayer=T,suffix='names',NAflag=-99,overwrite=T)
