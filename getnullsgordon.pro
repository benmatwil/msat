function getnullsgordon, gordon=gordon, me=me
  nullfile = '../../../crgordon/analysis/SPHR3/output/'
  if keyword_set(gordon) then nullfile = nullfile + 'nulls.dat'
  if keyword_set(me) then nullfile = nullfile + 'nullsben.dat'
  nullfile = "output/"
  openr,null,nullfile+"null.dat",/f77_unformatted, /get_lun;'../../../crgordon/analysis/SPHR3/output/null.dat'
  nnulls=0l
  readu,null,nnulls
  ;print,nnulls

  realpos = dblarr(3,nnulls)
  pos = realpos
  ;readu,null,realpos
  ;stop
  rp = dblarr(3)
  for i = 1, nnulls do begin
    readu,null,rp
    pos[*,i-1] = rp
  endfor
  for i = 1, nnulls do begin
    readu,null,rp
    realpos[*,i-1] = rp
  endfor

  close,null
  free_lun,null
  
  openr,null,nullfile+"nullsben.dat",/f77_unformatted,/get_lun;nullfile
  nnulls=0l
  readu,null,nnulls

  signs=lonarr(nnulls)
  r=dblarr(3,nnulls)
  spine=dblarr(3,nnulls)
  readu,null,signs
  readu,null,r,spine
  signs = -1*signs

  close,null
  free_lun,null

  nulls={nulldata,pos:dblarr(3),spine:dblarr(3),type:long(0)}
  nulls=replicate({nulldata},nnulls)

  for i = 0, nnulls-1 do begin
    ;nulls[i].pos = [realpos[0,i],realpos[2,i],!dpi/2-realpos[1,i]]
    nulls[i].pos = [pos[0,i],pos[1,i],pos[2,i]]
    nulls[i].spine = spine[*,i]
    nulls[i].type = signs[i]
  endfor
  
  return, nulls
end
