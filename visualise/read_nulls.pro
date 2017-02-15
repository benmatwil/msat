function read_nulls, file, simple=simple

  file1 = strmid(file, strlen('data/'), strlen(file))
  file1 = strmid(file1, 0, strlen(file1)-strlen('.dat'))
  openr, null, 'output/'+file1+'-nullpos.dat', /get_lun
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

    if not keyword_set(simple) then begin
      openr, null, 'output/'+file1+'-nulldata.dat', /get_lun
      nnulls = 0L
      readu, null, nnulls

      signs = lonarr(nnulls)
      spine = dblarr(3,nnulls)
      fan = spine
      readu,null,signs,spine,fan
      signs = -signs

      close, null
      free_lun, null
    endif

    nulls = {nulldata,pos:dblarr(3),gridpos:dblarr(3),spine:dblarr(3),fan:dblarr(3),sign:long(0)}
    nulls = replicate({nulldata},nnulls)

    for i = 0, nnulls-1 do begin
      nulls[i].pos = realpos[*,i]
      nulls[i].gridpos = pos[*,i]
      if not keyword_set(simple) then begin
        nulls[i].spine = spine[*,i]
        nulls[i].fan = fan[*,i]
        nulls[i].sign = signs[i]
      endif
    endfor
  endif

  return, nulls
  
end
