openr,10,'data/newmag.dat'
nx = 0l
ny = 0l
nz = 0l
readu,10,nx,ny,nz
bgrid = dblarr(nx,ny,nz,3)
x=dblarr(nx)
y=dblarr(ny)
z=dblarr(nz)
readu,10,bgrid,x,y,z
close,10

nulls = getnullsgordon(/me)

nx = 10
ny = 10

startpts = []

for i = 0, nx-1 do begin
  for j = 0, ny-1 do begin
    xp = 1.01
    yp = !dpi*(i+0.5)/nx
    zp = 2*!dpi*(j+0.5)/ny
    
    startpts = [[startpts],[xp,yp,zp]]
    print, xp,yp,zp
  endfor
endfor
help, startpts

for j = 0, n_elements(startpts)/3-1 do begin

  h = 1d-3
  line = fieldline3d(startpts[*,j], bgrid, x,y,z, h, 0.1d*h, 10d*h, 0.01d*h)

  line = transpose(line)

  if j eq 0 then overplot = 0 else overplot = 1
  plt = plot3d(line[*,0],line[*,1],line[*,2],overplot=overplot)
  
endfor

help, nulls

pos = transpose(nulls.pos)
plt = plot3d(pos[*,0],pos[*,1],pos[*,2],'r +',/overplot)
plt.xr = [1,1.01]
;plt


end
