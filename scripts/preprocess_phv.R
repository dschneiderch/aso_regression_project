library(raster)
library(dplyr)

ires='500m'
basin='tuo'

elev_rast=raster('data/elevation/tuo_dem500m.tif')
slope_rast=raster('data/elevation/tuo_slope500m.tif')
zness_rast=sin(slope_rast*pi/180)
writeRaster(zness_rast,'data/elevation/tuo_zness500m.tif',NAflag=-99,overwrite=T)

aspect_rast=raster('data/elevation/tuo_aspect500m.tif')
eastness_rast=mask(sin(aspect_rast*pi/180),elev_rast,filename='data/elevation/tuo_eastness500m.tif',NAflag=-99)
northness_rast=mask(cos(aspect_rast*pi/180),elev_rast,filename='data/elevation/tuo_northness500m.tif',NAflag=-99)

tpi_rast=raster('data/elevation/tuo_tpi500m.tif')
tri_rast=raster('data/elevation/tuo_tri500m.tif')
vrm_rast=raster('data/elevation/tuo_vrm1cell_dem500m.tif')

northing=as.data.frame(elev_rast,xy=T) %>% select(y)
northing_rast=mask(setValues(elev_rast,northing$y),elev_rast,filename='data/elevation/tuo_northing500m.tif',NAflag=-99)

easting=as.data.frame(elev_rast,xy=T) %>% select(x)
easting_rast=mask(setValues(elev_rast,easting$x),elev_rast,filename='data/elevation/tuo_easting500m.tif',NAflag=-99)

vegheight_rast=raster('data/veg/tuo_vegheight_500m.tif')

stdslope_rast=focal(slope_rast,w=matrix(1,3,3), fun=sd, na.rm=TRUE, filename='data/elevation/tuo_stdslope1cell_dem500m.tif')



