'reinit'
'set display white'
'clear'
'set grads off'

'open q_check_interp.ctl'
'sdfopen output.nc'

'set gxout scatter' 

'set vrange2 0 6000'
'set vrange 0 0.008'

'set ccolor 6'
'set cmark 3'
'set digsiz 0.2'
'd t;lev'
'set dfile 2'
'set ccolor 1'
'set cthick 10'
'set digsiz 0.3'
'set cmark 1'
'd q;lev'

'gxprint x.png'
'!convert x.png -trim -quality 100 q.png'
'!rm -f x.png'
