library(raster)
library(dplyr)

ires='100m'
basin='tuo'

basinmask=raster(paste0('data/gis/tuo_basinmask_',ires,'.tif'))

r=raster('data/elevation/tuo_dem99m.tif')
resample(r,basinmask,method='ngb',filename='data/elevation/tuo_dem100m.tif',NAflag=-99,overwrite=T)
elev_rast=raster(paste0('data/elevation/',basin,'_dem',ires,'.tif'))

slope_rast=raster(paste0('data/elevation/',basin,'_slope',ires,'.tif'))
zness_rast=sin(slope_rast*pi/180)
writeRaster(zness_rast,paste0('data/elevation/',basin,'_zness',ires,'.tif'),NAflag=-99,overwrite=T)

aspect_rast=raster(paste0('data/elevation/',basin,'_aspect',ires,'.tif'))
eastness_rast=mask(sin(aspect_rast*pi/180),elev_rast,filename=paste0('data/elevation/',basin,'_eastness',ires,'.tif'),NAflag=-99)
northness_rast=mask(cos(aspect_rast*pi/180),elev_rast,filename=paste0('data/elevation/',basin,'_northness',ires,'.tif'),NAflag=-99)

tpi_rast=raster(paste0('data/elevation/',basin,'_tpi',ires,'.tif'))
tri_rast=raster(paste0('data/elevation/',basin,'_tri',ires,'.tif'))
vrm_rast=raster(paste0('data/elevation/',basin,'_vrm1cell_dem',ires,'.tif'))

northing=as.data.frame(elev_rast,xy=T) %>% select(y) %>% tbl_df
northing_rast=mask(setValues(elev_rast,northing$y),elev_rast,filename=paste0('data/elevation/',basin,'_northing',ires,'.tif'),NAflag=-99)

easting=as.data.frame(elev_rast,xy=T) %>% select(x)
easting_rast=mask(setValues(elev_rast,easting$x),elev_rast,filename=paste0('data/elevation/',basin,'_easting',ires,'.tif'),NAflag=-99)

v=raster('data/veg/raw/tuo_vegheight_3m.tif')
va=aggregate(v,33,method='median',filename='data/veg/tuo_vegheight_99m.tif',NAflag=-99)
resample(va,elev_rast,method='ngb',filename='data/veg/tuo_vegheight_100m.tif',NAflag=-99,overwrite=T)
vegheight_rast=raster(paste0('data/veg/',basin,'_vegheight_',ires,'.tif'))

stdslope_rast=focal(slope_rast,w=matrix(1,3,3), fun=sd, na.rm=TRUE, filename=paste0('data/elevation/',basin,'_stdslope1cell_dem',ires,'.tif'))



