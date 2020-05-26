'reinit'
'set display white'
'clear'
'set grads off'


ncfile='output.nc'
*ncfile='output_fv_geo_real_w.nc'
* ncfile='output_ppm_geo_real_w.nc'
* ncfile='output_fv_geo_ideal_w.nc'
* ncfile='output_ppm_geo_ideal_w.nc'

taxis=off
var =  t

if ( taxis = on )
'sdfopen 'ncfile 

'set z 1 30'
'set t 1 last'
'set gxout grfill'
'set cmin 0'
'd 'var
endif

if ( taxis = off )
'sdfopen 'ncfile'' 

'set t last'
'q dims'
st=sublin(result,5);last=subwrd(st,9)
t=1
while(t<=last)
'c'
'set t 't
say t
'set z 1 50'
* 'set vrange 0.0 0.01'
* 'set vrange -3 3'   
tt=t+1 
* 'd 'var'(t='tt')-'var'(t='t')'
'd 'var
* 'gxprint x.png'
* '!convert x.png -trim -quality 100 'ncfile'_'t'.png'
* '!rm -f x.png'


t=t+1
endwhile


endif
