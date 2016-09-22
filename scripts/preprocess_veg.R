library(raster)

elev_rast=raster('data/elevation/tuo_dem50m.tif')

vegheight_rast=raster('data/veg/raw/tuo_vegheight_3m.tif')
aggregate(vegheight_rast,fact=17,method='median',filename='data/veg/tuo_vegheight_51m.tif',NAflag=-99,overwrite=T)
resample(raster('data/veg/tuo_vegheight_51m.tif'), elev_rast,method='ngb',filename='data/veg/tuo_vegheight_50m.tif',NAflag=-99,overwrite=T)


