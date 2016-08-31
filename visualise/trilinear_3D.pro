FUNCTION trilinear_3D,xp,yp,zp,bgrid,xx,yy,zz,periodic=periodic,jacobi=jacobi

;Finds magnetic field components a given point (xp,yp,zp)

;INPUTS
;(xp,yp,zp) point at which magfield required in x,y,z units
;bgrid[nx,ny,nz,3] - grid of magnetic field
;xx[nx],yy[ny],zz[nz] - x,y,z values on which magnetic field grid defined

if keyword_set(periodic) then begin
;SORT THIS OUT!!!!!!! 
endif

nx = n_elements(xx)
ny = n_elements(yy)
nz = n_elements(zz)

ii = where(xx lt xp,nii)
if (ii[0] eq -1) then ix = 0 else ix = ii(nii-1)
if (ix eq nx-1) then ix = nx-2 
jj = where(yy lt yp,njj)
if (jj[0] eq -1) then iy = 0 else iy = jj(njj-1)
if (iy eq ny-1) then iy = ny-2 
kk = where(zz lt zp,nkk)
if (kk[0] eq -1) then iz = 0 else iz = kk(nkk-1)
if (iz eq nz-1) then iz = nz-2 

dx = (xp-xx(ix))/(xx(ix+1)-xx(ix))
dy = (yp-yy(iy))/(yy(iy+1)-yy(iy))
dz = (zp-zz(iz))/(zz(iz+1)-zz(iz))

bp = dblarr(3)
jacobi = dblarr(3,3)

for k=0,2 do begin
  b000 = bgrid(ix,iy,iz,k)
  b100 = bgrid(ix+1,iy,iz,k)
  b010 = bgrid(ix,iy+1,iz,k)
  b110 = bgrid(ix+1,iy+1,iz,k)
  b001 = bgrid(ix,iy,iz+1,k)
  b101 = bgrid(ix+1,iy,iz+1,k)
  b011 = bgrid(ix,iy+1,iz+1,k)
  b111 = bgrid(ix+1,iy+1,iz+1,k)
  
  a = b000
  b = b100-b000
  c = b010-b000
  d = b110-b100-b010+b000
  e = b001-b000
  f = b101-b100-b001+b000
  g = b011-b010-b001+b000
  h = b111-b110-b101-b011+b100+b010+b001-b000

  bp(k) = a+b*dx+c*dy+d*dx*dy+e*dz+f*dx*dz+g*dy*dz+h*dx*dy*dz
  jacobi(0,k) = b+d*dy+f*dz+h*dy*dz
  jacobi(1,k) = c+d*dx+g*dz+h*dx*dz
  jacobi(2,k) = e+f*dx+g*dy+h*dx*dy
endfor

RETURN,bp

END
