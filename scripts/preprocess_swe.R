library(raster)
pathin='data/swe/raw'

basinFN <- function(x){
  switch(x,
         'TB'='tuo',
         'UB'='uncom')
}
basin='TB'
geoarea2=basinFN(basin)
newres='100m'
aggfactor=2#from 50m

basinmask=raster(paste0('data/gis/',geoarea2,'_basinmask_50m.tif'))

swefn=dir(pathin,pattern=glob2rx(paste0(basin,'*.tif$')),full.names = T)
# asodtes=substr(sapply(strsplit(basename(swefn),'[.\\_]'),'[',1),3,12)
# write.csv(data.frame(dte=asodtes),file='data/tuo_aso_flightdates.txt',row.names=F,quote=F)

# f=swefn[20]
# dte=20150427
for(f in swefn){
  dte=substr(sapply(strsplit(basename(f),'[.\\_]'),'[',1),3,10)
  r=raster(f)
  res(r)=c(50,50)
  if(!compareRaster(r,basinmask,extent=T,stopiffalse = F))  r=raster::resample(r,basinmask,method='bilinear')
  swe=r*basinmask
  dte=substr(sapply(strsplit(basename(f),'[.\\_]'),'[',1),3,10)
  writeRaster(swe,filename=paste0('data/swe/',geoarea2,'_50m_',dte,'.tif'),NAflag=-99,overwrite=T)
  # if(geoarea2=='uncom') swe[swe==0]=NA
  aggregate(swe,fact=aggfactor,method='mean',filename=paste0('data/swe/',geoarea2,'_',newres,'_',dte,'.tif'),NAflag=-99,overwrite=T)
}

r=raster('data/swe/raw/TB20130429_SUPERswe_50p0m_agg')
res(r)=c(50,50)
if(!compareRaster(r,basinmask,extent=T,stopiffalse = F))  r=raster::resample(r,basinmask,method='bilinear')
swe=r*basinmask
dte='20130429'
writeRaster(swe,filename=paste0('data/swe/',geoarea2,'_50m_',dte,'.tif'),NAflag=-99,overwrite=T)
# if(geoarea2=='uncom') swe[swe==0]=NA
aggregate(swe,fact=aggfactor,method='mean',filename=paste0('data/swe/',geoarea2,'_',newres,'_',dte,'.tif'),NAflag=-99,overwrite=T)

