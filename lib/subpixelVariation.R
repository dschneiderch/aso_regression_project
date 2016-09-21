subpixelVariation <- function(x,y,fxn='stats::sd'){
# ----- Get Stats of subpixels inside larger pixel
# x is small scale raster
# y is large scale raster
# fxn is function to apply to aggregate the smaller pixels
# returns a raster with cell size of y

## Create a Raster with the same resolution as small scale raster, whose values are
## the cell numbers of larger scale raster
xyDEM <- xyFromCell(x, cell = seq_len(ncell(x)))
zones <- setValues(x, cellFromXY(y, xyDEM))
## Compute the metric for each zone. `zoneVariation` has the same
## resolution as y
zoneVariation <- setValues(y, raster::zonal(x, zones, fun = fxn, na.rm=T)[,2])
return(zoneVariation)
}

