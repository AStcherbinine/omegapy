function testin,x0,y0,x1,y1
  nb=(size(x1))(1)
  x2=[x1,x1(0)]
  y2=[y1,y1(0)]
  dx=x2-x0
  dy=y2-y0
  atot=0.
  for n=0,nb-1 do begin
    ps=dx(n)*dx(n+1)+dy(n)*dy(n+1)
    pv=dx(n)*dy(n+1)-dx(n+1)*dy(n)
    atot=atot+atan(pv,ps)
  endfor
  if(abs(atot) gt 3) then goto, yes
  return, 0
yes:
  return, 1
end
    
