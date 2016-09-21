library(raster)
library(rgdal)
# this is for uncom and any other basin within the recon UCO domain.
# Sierra basins are currently already cropped small so just copy them from snowserver

## Make sure you have hydroProjects and hydroData mounted using smb
water.mask=raster('/Volumes/hydroProjects/SWE/Rockies/data/masks/cs_NHD_MOD44_water_mask_geo.tif')
projection(water.mask)=CRS('+proj=longlat +datum=NAD83')
outfile='data/fsca/raw'
if(!dir.exists(outfile)) dir.create(outfile,recursive=TRUE)

proc=function(fn,bd,dname){
  r=raster(fn)
  projection(r)=CRS('+proj=longlat +datum=NAD83')
  r <- mask(r,water.mask,maskvalue=NA,updatevalue=253)
  # r[is.na(water.mask)]=253
  # r=crop(r,bd)
  f=unlist(strsplit(fn,'/'))
  f=f[length(f)]
  print(paste0(outfile,dname,'_',f))
  writeRaster(r,filename=paste0(outfile,dname,f),format='GTiff',overwrite=T)
}




##Fish
#domain=readOGR('/Volumes/shareProjects/WWA/data/gis_dhsvm/boundary','fsh_domain')
#domain_fsh=spTransform(domain,CRS('+proj=longlat +datum=NAD83'))
#
##Snake
#domain=readOGR('/Volumes/shareProjects/WWA/data/gis_dhsvm/boundary','snk_domain')
#domain_snk=spTransform(domain,CRS('+proj=longlat +datum=NAD83'))
#
##Boulder
#domain=readOGR('/Volumes/shareProjects/WWA/data/gis_dhsvm/boundary','bld_domain')
#domain_bld=spTransform(domain,CRS('+proj=longlat +datum=NAD83'))
#
##Uncom
#domain=readOGR('/Volumes/shareProjects/WWA/data/gis_dhsvm/boundary','unc_domain')
domain=readOGR(paste0('data/gis/',geoarea2,'_extent.gpkg'),geoarea2)
domain_unc=spTransform(domain,CRS('+proj=longlat +datum=NAD83'))
yr=2012
for (yr in seq(2000,2013)){
  #yr=2010
  print(sprintf('getting modscag: %d',yr))
  ## pn.r='/Volumes/hydroData/WestUS_Data/MODSCAG/UpperColoradoRiver/Geographic/'
  ## ffns=list.files(path=pn.r,pattern=glob2rx(paste0(yr,'*.tif')),full.names=T)
  pn.r='/Volumes/hydroData/WestUS_Data/UCO_FSCA/'
  fdr=list.dirs(pn.r,full.names=T,recursive=F)
  dirs=grep(yr,fdr)
  if(dirs != ''){
    ffns=list.files(path=fdr[dirs],pattern=glob2rx('*.tif'),full.names=T)
    #
    #lapply(ffns,proc,domain_fsh,'fish')#
    #lapply(ffns,proc,domain_snk,'snake')
    #lapply(ffns,proc,domain_bld,'boulder')
    lapply(ffns,proc,domain_unc,'uncom')
  }
}

## f2s=list.files('~/Documents/wwa',pattern=glob2rx('fish2001*.tif'),full.names=T)
## fshstck=stack(f2s)

## findZero=function(lyr){
##     which(getValues(lyr)==0)
## }

## for(lyr in seq(1,nlayers(fshstck))){
## ans=which(getValues(fshstck[[lyr]])==0)
## print(ans)
## }

## findZero(fshstck[[10]])

## laply(fshstck,findZero)


## require(raster)
## r1 <- raster(nrows=10, ncols=10)
## r2=r3=r4=r1
## r1[]= runif(ncell(r1),75,100)
## r2[]= runif(ncell(r1),0,100)
## r3[]= runif(ncell(r1),0,100)
## r4[]= runif(ncell(r1),0,100)
## rs=stack(r1,r2,r3,r4)
## names(rs) <- c('2001060','2001073','2001084','2001100')


## test=calc(fshstck,fun=function(x,na.rm) x[order(x)])
## plot(test)
