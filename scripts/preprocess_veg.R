library(raster)

vegheight_rast=raster('data/veg/raw/tuo_vegheight_3m.tif')
aggregate(vegheight_rast,fact=167,method='median',filename='data/veg/tuo_vegheight_501m.tif',NAflag=-99,overwrite=T)
resample(raster('data/veg/tuo_vegheight_501m.tif'), elev_rast,method='ngb',filename='data/veg/tuo_vegheight_500m.tif',NAflag=-99,overwrite=T)
