function read_nulls

  openr, null, 'output/null.dat', /f77_unformatted, /get_lun
  nnulls = 0L
  readu, null, nnulls
  if nnulls gt 0 then begin
    realpos = dblarr(3,nnulls)
    pos = realpos
    rp = dblarr(3)
    for i = 1, nnulls do begin
      readu, null, rp
      pos[*,i-1] = rp
    endfor
    for i = 1, nnulls do begin
      readu, null, rp
      realpos[*,i-1] = rp
    endfor
    close,null
    free_lun, null

    openr, null, 'output/nulls.dat', /f77_unformatted, /get_lun
    nnulls = 0L
    readu, null, nnulls

    signs = lonarr(nnulls)
    r = dblarr(3,nnulls)
    spine = dblarr(3,nnulls)
    fan = spine
    readu,null,signs
    readu,null,r,spine,fan
    signs = -signs

    close, null
    free_lun, null

    nulls = {nulldata,pos:dblarr(3),gridpos:dblarr(3),spine:dblarr(3),fan:dblarr(3),sign:long(0)}
    nulls = replicate({nulldata},nnulls)

    for i = 0, nnulls-1 do begin
      nulls[i].pos = realpos[*,i]
      nulls[i].gridpos = pos[*,i]
      nulls[i].spine = spine[*,i]
      nulls[i].fan = fan[*,i]
      nulls[i].sign = signs[i]
    endfor
  endif

  return, nulls
  
end
