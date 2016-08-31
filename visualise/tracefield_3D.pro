FUNCTION tracefield_3D,startpt,bgrid,x,y,z,h,hmin,hmax,epsilon,mxline

;startpt[3,nl] - start points for field line
;bgrid[nx,ny,nz,3] - magnetic field 
;x[nx],y[ny],z[nz] - grid of points on which mag field given 

;h - initial step length
;hmin - minimum step length
;hmax - maximum step length
;epsilon - tolerance to which we require point on field line known

;define edges of box
xmin = min(x)
ymin = min(y)
zmin = min(z)
xmax = max(x)
ymax = max(y)
zmax = max(z)

a2 = 0.25
b2 = 0.25
a3 = 3./8. 
b3 = 3./32.
c3 = 9./32.
a4 = 12./13. 
b4 = 1932./2197.
c4 = -7200./2197.
d4 = 7296./2197.
a5 = 1.
b5 = 439./216. 
c5 = -8.
d5 = 3680./513.
e5 = -845./4104.
a6 = 0.5 
b6 = -8./27
c6 = 2. 
d6 = -3544./2565.
e6 = 1859./4104.
f6 = -11./40.

;used to determine y_i+1 from y_i if using rkf45 (4th order)
n1 = 25./216.
n3 = 1408./2565.
n4 = 2197./4104.
n5 = -1./5.
  
;used to determine y_i+1 from y_i if using rkf54 (5th order)
nn1 = 16./135.
nn3 = 6656./12825.
nn4 = 28561./56430.
nn5 = -9./50
nn6 = 2./55.

cont = 1   ;flag that if 1 means determine next point on field line
           ;          if 0 means that field line has reached some user
           ;          defined stopping criteria,
           ;          e.g. field line has left the domain 

If (h lt 0) then begin
  minh = -hmax
  maxh = -hmin
  hmin = minh
  hmax = maxh
endif
;mxline = 5000
line = dblarr(3,mxline)-9999.
s = dblarr(mxline)
s[0] = 0.
line[*,0] = startpt
count = 0
while(cont eq 1) do begin
t = s[count]
xx = line[0,count]
yy = line[1,count]
zz = line[2,count]

;print,count,xx,yy

xx1 = xx
yy1 = yy
zz1 = zz
bbf = trilinear_3D(xx1,yy1,zz1,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k1 = h*bbf[0]/bbf2
l1 = h*bbf[1]/bbf2
m1 = h*bbf[2]/bbf2
xx2 = xx+b2*k1
yy2 = yy+b2*l1
zz2 = zz+b2*m1
if (xx2 gt xmin and xx2 lt xmax and yy2 gt ymin and yy2 lt ymax and zz2 gt zmin and zz2 lt zmax) then begin &$
  cont = 1 &$
   endif else begin &$
  cont = 0 &$
  goto,checkcont &$
endelse
bbf = trilinear_3D(xx2,yy2,zz2,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k2 = h*bbf[0]/bbf2
l2 = h*bbf[1]/bbf2
m2 = h*bbf[2]/bbf2
xx3 = xx+b3*k1+c3*k2
yy3 = yy+b3*l1+c3*l2
zz3 = zz+b3*m1+c3*m2
if (xx3 gt xmin and xx3 lt xmax and yy3 gt ymin and yy3 lt ymax and zz3 gt zmin and zz3 lt zmax) then begin &$
  cont = 1 &$
   endif else begin &$
  cont = 0 &$
  goto,checkcont &$
endelse
bbf = trilinear_3D(xx3,yy3,zz3,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k3 = h*bbf[0]/bbf2
l3 = h*bbf[1]/bbf2
m3 = h*bbf[2]/bbf2
xx4 = xx+b4*k1+c4*k2+d4*k3
yy4 = yy+b4*l1+c4*l2+d4*l3
zz4 = zz+b4*m1+c4*m2+d4*m3
if (xx4 gt xmin and xx4 lt xmax and yy4 gt ymin and yy4 lt ymax and zz4 gt zmin and zz4 lt zmax) then begin &$
  cont = 1 &$
   endif else begin &$
  cont = 0 &$
  goto,checkcont &$
endelse
bbf = trilinear_3D(xx4,yy4,zz4,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k4 = h*bbf[0]/bbf2
l4 = h*bbf[1]/bbf2
m4 = h*bbf[2]/bbf2
xx5 = xx+b5*k1+c5*k2+d5*k3+e5*k4
yy5 = yy+b5*l1+c5*l2+d5*l3+e5*l4
zz5 = zz+b5*m1+c5*m2+d5*m3+e5*m4
if (xx5 gt xmin and xx5 lt xmax and yy5 gt ymin and yy5 lt ymax and zz5 gt zmin and zz5 lt zmax) then begin &$
  cont = 1 &$
   endif else begin &$
  cont = 0 &$
  goto,checkcont &$
endelse
bbf = trilinear_3D(xx5,yy5,zz5,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k5 = h*bbf[0]/bbf2
l5 = h*bbf[1]/bbf2
m5 = h*bbf[2]/bbf2
xx6 = xx+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5
yy6 = yy+b6*l1+c6*l2+d6*l3+e6*l4+f6*l5
zz6 = zz+b6*m1+c6*m2+d6*m3+e6*m4+f6*m5
if (xx6 gt xmin and xx6 lt xmax and yy6 gt ymin and yy6 lt ymax and zz6 gt zmin and zz6 lt zmax) then begin &$
  cont = 1 &$
   endif else begin &$
  cont = 0 &$
  goto,checkcont &$
endelse
bbf = trilinear_3D(xx6,yy6,zz6,bgrid,x,y,z)
bbf2 = sqrt(bbf[0]*bbf[0]+bbf[1]*bbf[1]+bbf[2]*bbf[2])
k6 = h*bbf[0]/bbf2
l6 = h*bbf[1]/bbf2
m6 = h*bbf[2]/bbf2

xxnew4 = xx + n1*k1+n3*k3+n4*k4+n5*k5
yynew4 = yy + n1*l1+n3*l3+n4*l4+n5*l5
zznew4 = zz + n1*m1+n3*m3+n4*m4+n5*m5
xxnew5 = xx + nn1*k1+nn3*k3+nn4*k4+nn5*k5+nn6*k6
yynew5 = yy + nn1*l1+nn3*l3+nn4*l4+nn5*l5+nn6*l6
zznew5 = zz + nn1*m1+nn3*m3+nn4*m4+nn5*m5+nn6*m6

;optimum stepsize
err= sqrt((xxnew5-xxnew4)^2.+(yynew5-yynew4)^2.+(zznew5-zznew4)^2.)
t = (epsilon*abs(h)/(2.*err))^0.25

;print,'err,  t',err,t

if (t lt 0.1) then t = 0.1
if (t gt 4) then t = 4 

h = t*h
if (h lt hmin) then h = hmin
if (h gt hmax) then h = hmax
count = count+1
if (count ge mxline) then begin
   cont = 0
   goto,checkcont
endif
s[count] = t
line[*,count] = [xxnew4,yynew4,zznew4]

;check that line is still in domain
if (xxnew4 gt xmin and xxnew4 lt xmax and yynew4 gt ymin and yynew4 lt ymax and zznew4 gt zmin and zznew4 lt zmax) then cont=1 else cont=0
;checkline is still moving
;dlx1 = line[0,count-1]-line[0,count-2]
;dly1 = line[1,count-1]-line[1,count-2]
;dlz1 = line[2,count-1]-line[2,count-2]
;line_spc1 = sqrt(dlx1*dlx1+dly1*dly1+dlz1*dlz1)
;if (line_spc1 gt hmin*0.5) then cont=1 else cont=0
if (count ge 2) then begin
  dlx2 = line[0,count-1]-line[0,count-2]
  dly2 = line[1,count-1]-line[1,count-2]
  dlz2 = line[2,count-1]-line[2,count-2]
  line_spc2 = sqrt(dlx2*dlx2+dly2*dly2+dlz2*dlz2)
  if (line_spc2 gt hmin*0.5) then cont=1 else cont=0
;  if (line_spc1 gt hmin*0.5 and line_spc2 gt hmin*0.5) then cont=1 else cont=0
endif
checkcont: 

;print,'x range',xmin,xmax
;print,xx1,xx2,xx3,xx4,xx5,xx6
;print,'y range',ymin,ymax
;print,yy1,yy2,yy3,yy4,yy5,yy6

;print,'check',count,h,t,xxnew4,yynew4,cont

end

;print,count,xxnew4,yynew4

If (h lt 0) then begin
  minh = -hmax
  maxh = -hmin
  hmin = minh
  hmax = maxh
endif

return,line

END
