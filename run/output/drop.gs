'reinit'
'set display white'
'clear'
'set grads off'

'sdfopen drop.nc'

'set t last'
'q dims'
st=sublin(result,5);last=subwrd(st,9)
t=1
while(t<=last)
'set x 1 50'
'set z 1'
'set t 't

'd num'

t=t+100
endwhile
