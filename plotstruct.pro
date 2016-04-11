n=0l
nfan = 0l

foreach name, ['fandata.dat', 'spinedata.dat'] do begin
  openr,10,name
  readu,10,n
  r = dblarr(3,n)
  readu,10,r
  help,r
  if name eq 'fandata.dat' then begin
    p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),lines=' ',sym='x', xr=[-1,1], yr=[-1,1],zr=[-1,1], aspect_ratio=1)
    readu,10,nfan
    nfan--
    p = plot3d(reform(r[0,0]),reform(r[1,0]),reform(r[2,0]),'cyan',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
    p = plot3d(reform(r[0,nfan]),reform(r[1,nfan]),reform(r[2,nfan]),'g',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
    nfan=-1
    p = plot3d(reform(r[0,nfan]),reform(r[1,nfan]),reform(r[2,nfan]),'g',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
  endif
  if name eq 'spinedata.dat' then begin
    p = plot3d(reform(r[0,*]),reform(r[1,*]),reform(r[2,*]),'r',/overplot)
  endif
  close,10
endforeach

end
