'reinit'
'set display white'
'clear'
'set grads off'

taxis=off
var = q 

ncfile1='output_fv_geo_real_w.nc'
ncfile2='output_ppm_geo_real_w.nc'
* ncfile='output_fv_geo_ideal_w.nc'
* ncfile='output_ppm_geo_ideal_w.nc'

if ( taxis = on )
'sdfopen 'ncfile1'' 
'sdfopen 'ncfile2'' 
'set z 1 100'
'set t 1 last'
'd 'var'.1'
'd 'var'.2'
'd 'var'.3'
endif

if ( taxis = off )
'sdfopen 'ncfile1'' 
'sdfopen 'ncfile2'' 

'set t last'
'q dims'
st=sublin(result,5);last=subwrd(st,9)
t=1
while(t<=last)
'c'
'set grads off'
'set t 't
say t
'set z 1 30'
'set vrange 0 0.007'
'set xlint 0.001'
'set digsiz 0.05'
'set grid off' 
'd 'var'.1'
'set digsiz 0.0'
'd 'var'.2'
'set xlpos 0 t'
'set grid off' 

'set vrange -0.00001 0.00001'
'set xlint 0.000005'
'd 'var'.1-'var'.2'

'set string 1 l 1 0';'set strsiz 0.13 0.13';'d sum('var'.1,z=1,z=30)';a=sublin(result,2);b=subwrd(a,4)
                    'draw string 6.1 7.5 `5TOTAL-FEFV: `5'substr(b,1,5)'' ;
'set string 1 l 1 0';'set strsiz 0.13 0.13';'d sum('var'.2,z=1,z=30)';a=sublin(result,2);b=subwrd(a,4)
                    'draw string 6.1 7.2 `5TOTAL-PPM: `5'substr(b,1,5)'' ;

'gxprint x.png'
'!convert x.png -trim -quality 100 real_'var'_'t'.png'
'!rm -f x.png'

t=t+1
endwhile


endif
