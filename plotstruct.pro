n=0l
spine = dblarr(3)
maxvec = spine
minvec = spine

foreach name, ['fandata.dat', 'spinedata.dat'] do begin
  openr,10,name
  readu,10,n
  r = dblarr(3,n)
  readu,10,r
  help,r
  if name eq 'fandata.dat' then begin
    p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),lines=' ',sym='x', xr=[-1,1], yr=[-1,1],zr=[-1,1], aspect_ratio=1)
    readu,10,maxvec,minvec
    p = plot3d([maxvec[0]],[maxvec[1]],[maxvec[2]],'cyan',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
    p = plot3d([minvec[0]],[minvec[1]],[minvec[2]],'g',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
    ;readu, 10, n & crossfan = dblarr(3,n) & readu, 10, crossfan & help, crossfan
    ;p = plot3d(reform(crossfan[0,*]),reform(crossfan[1,*]),reform(crossfan[2,*]),'orange',lines=' ',sym='x',/overplot, sym_thick=3, sym_size=2)
    ;cross = dblarr(3,n_elements(r)-1)
    ;for i = 1, n_elements(r)/3-1 do begin
    ;  cross[*,i-1] = crossp(r[*,0],r[*,i])
    ;  cross[*,i-1] = cross[*,i-1]/sqrt(cross[0,i-1]^2+cross[1,i-1]^2+cross[2,i-1]^2)
    ;endfor
    ;print, cross
    ;p = plot3d(reform(cross[0,*]),reform(cross[1,*]),reform(cross[2,*]),lines=' ',sym='x',/overplot, sym_thick=4, sym_size=3)
  endif
  if name eq 'spinedata.dat' then begin
    readu, 10, spine
    p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),'r',lines=' ',sym='x',/overplot)
    p = plot3d([spine[0]],[spine[1]],[spine[2]],'blue',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
  endif
  close,10
endforeach

end
