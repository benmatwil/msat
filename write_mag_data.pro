PRO write_mag_data,a,b,c,aa,bb,cc,dd,tt,box=box,ngrid=ngrid

;write_mag_data,1.,1.,-2.,1.,1.,-2.,5.,0.16
;write_mag_data,2.,-1.,3.,1.,-1.,-1.,0.,0.2

If keyword_set(ngrid) then n = ngrid else n = 200l
If keyword_set(box) then bbox=box else bbox=[[-5,5],[-5,5],[-5,5]]

nx = n
ny = n
nz = n

x = dblarr(nx)
y = dblarr(ny)
z = dblarr(nz)

x = bbox[0,0] + (bbox[1,0]-bbox[0,0])*dindgen(nx)/double(nx-1) 
y = bbox[0,1] + (bbox[1,1]-bbox[0,1])*dindgen(ny)/double(ny-1) 
z = bbox[0,2] + (bbox[1,2]-bbox[0,2])*dindgen(nz)/double(nz-1) 

bx = dblarr(nx,ny,nz)
by = dblarr(nx,ny,nz)
bz = dblarr(nx,ny,nz)

;nt = 9
;tim = dblarr(nt)

;for t=0,nt - 1 do begin
;  tim = 0.02*(t+1)
for k=0,n-1 do begin &$
  for j=0,n-1 do begin &$
    bx[*,j,k] = (a-z[k])*x + b*y[j] + aa*tt &$
    by[*,j,k] = c*x - (a+z[k])*y[j] + bb*tt &$
    bz[*,j,k] = dd*tt*y[j] + z[k]*z[k] + cc*tt &$
 endfor &$
endfor
;fnum = strmid(strcompress(string(tt*1000),/remove_all),1,3)
get_lun,lun
openw,lun,"data/bmag000.dat"
writeu,lun, nx,ny,nz
writeu,lun, bx,by,bz
writeu,lun, x,y,z
close,lun
free_lun,lun
;endfor

end
