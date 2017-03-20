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

function read_separators, filename, conlist=conlist

  nulldata = read_nulls(filename, /simple)

  openr, con, 'output/'+prefix(filename)+'-connectivity.dat', /get_lun
  openr, sep, 'output/'+prefix(filename)+'-separators.dat', /get_lun

  conlist = list()
  seplist = list()
  inull = 1
  coni = list()
  sepi = list()

  start = 0L
  end1 = 0L
  length = 0L

  readu, con, start
  while start ge 0 or inull le max(nulldata.number) do begin
    if inull eq start then begin
      readu, con, end1
      skip_lun, con, 8
      readu, sep, length
      separator = dblarr(3, length)
      readu, sep, separator
      coni = coni + list(end1)
      sepi = sepi + list(separator)
      readu, con, start
    endif
    if inull ne start then begin
      inull++
      conlist = conlist + list(coni)
      seplist = seplist + list(sepi)
      sepi = list()
      coni = list()
    endif
  endwhile

  close, sep, con
  free_lun, sep, con

  return, seplist

end

function read_spines, filename
  
  nulldata = read_nulls(filename, /simple)

  spinelist = list()
  length = 0L

  openr, spi, 'output/'+prefix(filename)+'-spines.dat', /get_lun

  foreach inull, nulldata.number do begin
    spinelisti = list()
    for ispine = 1, 2 do begin
      readu, spi, length
      spine = dblarr(3, length)
      readu, spi, spine
      spinelisti = spinelisti + list(spine)
    endfor
    spinelist = spinelist + list(spinelisti)
  endforeach

  close, spi
  free_lun, spi
  
  return, spinelist

end

function read_rings, filename, nskip=nskip, breaks=breaklist

  nulldata = read_nulls(filename, /simple)

  if not keyword_set(nskip) then nskip = 1

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
    iring = 0
    while iring lt ringsmax do begin
      if lengths[iring] eq 0 then break
      assocs = lonarr(lengths[iring])
      breaks = lonarr(lengths[iring])
      rings = dblarr(3,lengths[iring])
      readu, ring, assocs, breaks, rings
      iskip = 1
      while iskip + iring lt ringsmax and iskip lt nskip do begin
        skip_lun, ring, lengths[iring+iskip]*32L
        iskip++
      endwhile
      iring = iring + nskip
      assoclisti = assoclisti + list(assocs)
      breaklisti = breaklisti + list(breaks)
      ringlisti = ringlisti + list(rings)
    endwhile
    assoclist = assoclist + list(assoclisti)
    ringlist = ringlist + list(ringlisti)
    breaklist = breaklist + list(breaklisti)
  endforeach

  close, rinfo, ring
  free_lun, rinfo, ring
  
  return, ringlist

end
