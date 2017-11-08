function prefix, filename

    file1 = strmid(filename, strlen('data/'), strlen(filename))
    file1 = strmid(file1, 0, strlen(file1)-strlen('.dat'))

    return, file1

end

function read_nulls, filename, simple=simple

  openr, null, 'output/'+prefix(filename)+'-nullpos.dat', /get_lun
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
      openr, null, 'output/'+prefix(filename)+'-nulldata.dat', /get_lun
      nnulls = 0L
      readu, null, nnulls

      signs = lonarr(nnulls)
      spine = dblarr(3,nnulls)
      fan = spine
      readu,null,signs,spine,fan

      close, null
      free_lun, null
    endif

    nulls = {nulldata, pos:dblarr(3), gridpos:dblarr(3), spine:dblarr(3), fan:dblarr(3), sign:0, number:0}
    nulls = replicate({nulldata},nnulls)

    for i = 0, nnulls-1 do begin
      nulls[i].pos = realpos[*,i]
      nulls[i].gridpos = pos[*,i]
      nulls[i].number = i+1
      if not keyword_set(simple) then begin
        nulls[i].spine = spine[*,i]
        nulls[i].fan = fan[*,i]
        nulls[i].sign = signs[i]
      endif
    endfor
  endif

  return, nulls
  
end

function read_separators, filename, conlist=conlist, null_list=null_list, hcs=hcs

  nulldata = read_nulls(filename, /simple)

  if not keyword_set(hcs) then begin
    openr, con, 'output/'+prefix(filename)+'-connectivity.dat', /get_lun
    openr, sep, 'output/'+prefix(filename)+'-separators.dat', /get_lun
    allnulls_list = nulldata.number
  endif else begin
    openr, con, 'output/'+prefix(filename)+'-hcs-connectivity.dat', /get_lun
    openr, sep, 'output/'+prefix(filename)+'-hcs-separators.dat', /get_lun
    allnulls_list = [1, 2]
  endelse

  if not keyword_set(null_list) then null_list = allnulls_list

  conlist = list()
  seplist = list()

  start = 0L
  endnum = 0L
  length = 0L

  readu, con, start
  foreach inull, allnulls_list do begin
    coni = list()
    sepi = list()
    while start eq inull do begin
      readu, sep, length
      if where(inull eq null_list, /null) ne !null then begin
        readu, con, endnum
        coni.add, endnum
        skip_lun, con, 8
        separator = dblarr(3, length)
        readu, sep, separator
        sepi.add, separator
      endif else begin
        skip_lun, con, 12
        skip_lun, sep, 3LL*8LL*length
      endelse
      readu, con, start
    endwhile
    conlist.add, coni
    seplist.add, sepi
  endforeach

  close, sep, con
  free_lun, sep, con

  return, seplist

end

function read_spines, filename, null_list=null_list
  
  nulldata = read_nulls(filename, /simple)

  if not keyword_set(null_list) then null_list = nulldata.number

  spinelist = list()
  length = 0L

  openr, spi, 'output/'+prefix(filename)+'-spines.dat', /get_lun

  foreach inull, nulldata.number do begin
    spinelisti = list()
    for ispine = 1, 2 do begin
      readu, spi, length
      if where(inull eq null_list, /null) ne !null then begin
        spine = dblarr(3, length)
        readu, spi, spine
        spinelisti.add, spine
      endif else begin
        skip_lun, spi, 3LL*8LL*length
      endelse
    endfor
    spinelist.add, spinelisti
  endforeach

  close, spi
  free_lun, spi
  
  return, spinelist

end

function read_rings, filename, nskip=nskip, breaks=breaklist, null_list=null_list

  nulldata = read_nulls(filename, /simple)

  if not keyword_set(nskip) then nskip = 1

  if not keyword_set(null_list) then null_list = nulldata.number

  ringlist = list()
  breaklist = list()
  assoclist = list()

  openr, rinfo, 'output/'+prefix(filename)+'-ringinfo.dat', /get_lun
  openr, ring, 'output/'+prefix(filename)+'-rings.dat', /get_lun

  ringsmax = 0L

  readu, rinfo, ringsmax
  foreach inull, nulldata.number do begin
    breaklisti = list()
    ringlisti = list()
    assoclisti = list()
    lengths = lonarr(ringsmax)
    readu, rinfo, lengths
    if where(inull eq null_list, /null) ne !null then begin
      foreach length, lengths[where(lengths ne 0)], iring do begin
        if iring mod nskip eq 0 then begin
          assocs = lonarr(length)
          breaks = lonarr(length)
          rings = dblarr(3, length)
          readu, ring, assocs, breaks, rings
          breaklisti.add, breaks
          ringlisti.add, rings
          assoclisti.add, assocs
        endif else begin
          skip_lun, ring, 32LL*length
        endelse
      endforeach
    endif else begin
      skip_lun, ring, 32LL*total(lengths, /integer)
    endelse
    ringlist.add, ringlisti
    breaklist.add, breaklisti
    assoclist.add, assoclisti
  endforeach

  close, rinfo, ring
  free_lun, rinfo, ring
  
  return, ringlist

end
