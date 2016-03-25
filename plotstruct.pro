openr,10,'fandata.dat' & readu,10,n & r = dblarr(3,n) & readu,10,r & close,10 & help,r
p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),lines=' ',sym='x', xr=[-1,1], yr=[-1,1],zr=[-1,1], aspect_ratio=1)

openr,10,'spinedata.dat' & readu,10,n & r = dblarr(3,n) & readu,10,r & close,10 & help,r
p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),'r',/overplot)
