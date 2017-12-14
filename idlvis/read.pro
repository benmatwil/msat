function prefix, filename

    file1 = strmid(filename, strlen('data/'), strlen(filename))
    file1 = strmid(file1, 0, strlen(file1)-strlen('.dat'))

    return, file1

end

function read_field, filename

  openr, field, filename, /get_lun
  nx = 0L
  ny = 0L
  nz = 0L
  readu, field, nx, ny, nz
  bx = dblarr(nx, ny, nz)
  by = bx
  bz = bx
  x = dblarr(nx)
  y = dblarr(ny)
  z = dblarr(nz)
  readu, field, bx, by, bz, x, y, z
  close, field
  free_lun, field

  return, {bx:bx, by:by, bz:bz, x:x, y:y, z:z}

end

function read_nulls, filename, simple=simple

  openr, null, 'output/'+prefix(filename)+'-nullpos.dat', /get_lun
  nnulls = 0L
  readu, null, nnulls
  if nnulls gt 0 then begin
    realpos = dblarr(3, nnulls)
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
    close, null
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
    allnulls_list = [1]
  endelse

  if not keyword_set(null_list) then null_list = allnulls_list

  conlist = list()
  seplist = list()

  endnum = 0L
  length = 0L
  nseps = 0L

  foreach inull, allnulls_list do begin
    coni = list()
    sepi = list()
    readu, con, nseps
    for i = 1, nseps do begin
      skip_lun, con, 4
      readu, con, endnum
      coni.add, endnum
      skip_lun, con, 8
      readu, sep, length
      separator = dblarr(3, length)
      readu, sep, separator
      sepi.add, separator
    endfor
    if where(inull eq null_list, /null) ne !null then begin
      conlist.add, coni
      seplist.add, sepi
    endif else begin
      conlist.add, list()
      seplist.add, list()
    endelse
    
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

function read_rings, filename, nskip=nskip, breaks=breaklist, assocs=assoclist, null_list=null_list, hcs=hcs

  nulldata = read_nulls(filename, /simple)

  if not keyword_set(nskip) then nskip = 1

  if keyword_set(hcs) then begin
    assocs_filename = 'output/'+prefix(filename)+'-hcs-assocs.dat'
    info_filename = 'output/'+prefix(filename)+'-hcs-ringinfo.dat'
    ring_filename = 'output/'+prefix(filename)+'-hcs-rings.dat'
    break_filename = 'output/'+prefix(filename)+'-hcs-breaks.dat'
  endif else begin
    assocs_filename = 'output/'+prefix(filename)+'-assocs.dat'
    info_filename = 'output/'+prefix(filename)+'-ringinfo.dat'
    ring_filename = 'output/'+prefix(filename)+'-rings.dat'
    break_filename = 'output/'+prefix(filename)+'-breaks.dat'
  endelse

  if keyword_set(assocs) and not file_test(assocs_filename) then print, 'Not reading in associations: file does not exist'
  if keyword_set(assocs) and file_test(assocs_filename) then do_assocs = 1 else do_assocs = 0

  openr, rinfo, info_filename, /get_lun
  openr, ring, ring_filename, /get_lun
  openr, brks, break_filename, /get_lun
  if do_assocs then openr, ass, 'output/'+prefix(filename)+'-assocs.dat', /get_lun

  ringsmax = 0L
  writeskip = 0L
  bytesize = 0L
  stepsize = 0D
  n_hcs = 0L

  readu, rinfo, ringsmax, writeskip, bytesize, stepsize
  if keyword_set(hcs) then readu, rinfo, n_hcs
  if writeskip gt 1 then ringsmax = ciel(double(ringsmax)/writeskip)

  if keyword_set(hcs) then begin
    to_do = indgen(n_hcs/2) + 1
    ndir = 2
    name = 'hcs'
  endif else begin
    to_do = nulldata.number
    ndir = 1
    name = 'null'
  endelse

  if not keyword_set(null_list) and not keyword_set(hcs) then begin
    null_list = nulldata.number
  endif else if keyword_set(hcs) then begin
    null_list = to_do
  endif

  ringlist = list()
  breaklist = list()
  if do_assocs then assoclist = list()

  foreach inull, to_do do begin
    for dir = 1, ndir do begin
      breaklisti = list()
      ringlisti = list()
      if do_assocs then assoclisti = list()
      lengths = lonarr(ringsmax)
      readu, rinfo, lengths
      if where(inull eq null_list, /null) ne !null then begin
        foreach length, lengths[where(lengths ne 0)], iring do begin
          if iring mod nskip eq 0 then begin
            if do_assocs then assocs = lonarr(length)
            breaks = lonarr(length)
            if bytesize eq 4 then begin
              rings = fltarr(3, length)
            endif else if bytesize eq 8 then begin
              rings = dblarr(3, length)
            endif
            readu, ring, rings
            if do_assocs then readu, ass, assocs
            readu, brks, breaks
            breaklisti.add, breaks
            ringlisti.add, rings
            if do_assocs then assoclisti.add, assocs
          endif else begin
            skip_lun, ring, 3LL*bytesize*length
            skip_lun, brks, 4LL*length
            if do_assocs then skip_lun, ass, 4LL*length
          endelse
        endforeach
      endif else begin
        total_len = total(lengths, /integer)
        skip_lun, ring, 3LL*bytesize*total_len
        skip_lun, brks, 4LL*total_len
        if do_assocs then skip_lun, ass, 4LL*total_len
      endelse
      ringlist.add, ringlisti
      breaklist.add, breaklisti
      if do_assocs then assoclist.add, assoclisti
    endfor
  endforeach

  close, rinfo, ring, brks
  free_lun, rinfo, ring, brks
  if do_assocs then begin
    close, ass
    free_lun, ass
  endif
  
  return, ringlist

end
