'reinit'
*'set display white'
'clear'
'set grads off'

'sdfopen drop.nc'

'set t last'
'q dims'
st=sublin(result,5);last=subwrd(st,9)
t=1
while(t<=last)
say t

'set z 1 100'
'set t 't
'set vrange 0 2000'
'set cmark 0'
'd sum(num,x=1,x=100)'

t=t+1
endwhile
