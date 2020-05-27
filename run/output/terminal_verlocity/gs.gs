'reinit'
'set display white'
'clear'
'set grads off'

'open 1976.ctl'
'open 1977.ctl'
'open r.ctl'

* 'set vrange 0 0.4'
'set gxout scatter'
* 'set vrange 0 0.001'
* 'd r.3;((vt.2/vt.1)-1)'
'd vt'


'gxprint x.png'
'!convert x.png -trim -quality 100 vt.png'
'!rm -f x.png'
