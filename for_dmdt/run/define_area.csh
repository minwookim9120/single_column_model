#!/bin/csh -f

set slon = 129  
set elon = 129.5
set slat = 37  
set elat = 37.5

set obj_data = "../input/e5.mnth.mean.an.pl.128_129_z.regn320sc.2017010100_2017120100.nc"


cat > find_index.ncl << EOF
fh     = addfile("${obj_data}", "r")

lat_1  = fh->latitude
lon_1  = fh->longitude

rm = new( (/1,2/), "float")

lat = (/ ${slat}, ${elat} /)
lon = (/ ${slon}, ${elon} /)

i_rlat = ind_nearest_coord(lat, lat_1, 0)
i_rlon = ind_nearest_coord(lon, lon_1, 0)

   print("slat= "+i_rlat(0) )
   print("elat= "+i_rlat(1) )
   print("slon= "+i_rlon(0) )
   print("elon= "+i_rlon(1) )
   print("i="+i_rlat+"  want_lat="+lat+"  lat(i)="+lat_1(i_rlat))
   print("i="+i_rlon+"  want_lon="+lon+"  lon(i)="+lon_1(i_rlon)) 

EOF

ncl find_index.ncl > index.txt
ncl find_index.ncl 

set xstart_index = `cat index.txt | awk '/slat/{print $3}'`
set xend_index = `cat index.txt | awk '/elat/{print $3}'`
set ystart_index = `cat index.txt | awk '/slon/{print $3}'`
set yend_index = `cat index.txt | awk '/elon/{print $3}'`

sed "s/@slat@/${xstart_index}/g" namelist.info_in > namelist.1
sed "s/@elat@/${xend_index}/g"   namelist.1 > namelist.2
sed "s/@slon@/${ystart_index}/g" namelist.2 > namelist.3 
sed "s/@elon@/${yend_index}/g"   namelist.3 > namelist.info 


rm -f namelist.1 namelist.2 namelist.3
rm -f find_index.ncl index.txt 
