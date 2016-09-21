# This script exclusively works with data downloaded from JPL Snow.
# It will mosaic and subset modis tiles from modscag-historic oon snow server
# Run this from your computer in scripts folder. Currently set up for Mac with hydroData mounted as a server.

pnmain=/Volumes/hydroData/WestUS_Data/MODSCAG
domain=Tuolumne
pin=modscag-historic
pnout=data/fsca/raw/sinu
mkdir -p ../$pnout

for yr in {2000..2015}; do
pn=$pnmain/$pin/$yr


for doy in {001..366}; do
doy=`printf "%003d" ${doy}`

if [ -f $pnout/${yr}${doy}.tif ]; then
rm $pnout/${yr}${doy}.tif
fi

fn=`ls $pnmain/$pin/$yr/$doy/*snow_fraction.tif`

case "$domain" in
  Tuolumne)
    #echo $domain
    gdalwarp -te -10547146 4189710 -10405042 4252215 -dstnodata -99 -ot Int16 ${fn} ../$pnout/tuo_${yr}${doy}.tif
    ;;
  Uncompahgre)
    echo $domain
    ;;

esac
#rm $doy/*.1 $doy/*.2

done # doy loop

done # year loop
