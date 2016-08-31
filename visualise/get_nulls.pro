pro get_nulls,nnulls,signs,r,spine,nulls

  get_lun, null
  openr, null, 'output/nulls.dat', /f77_unformatted
  nnulls = 0l
  readu, null, nnulls
  nulls = 0
  if (nnulls gt 0) then begin
    signs = lonarr(nnulls)
    r = dblarr(3,nnulls)
    spine = dblarr(3,nnulls)
    readu, null, signs
    readu, null, r, spine
    
    signs = signs*(-1)
    
    print, 'NOTE: HARDWIRED TO USE GRIDCELL COORDINATES'
    nulls = 1
  endif
  close, null
  free_lun, null
  
end
