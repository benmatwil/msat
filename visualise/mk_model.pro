common shared_var, bgrid, xx, yy, zz, oModel, null, inull, ds, filename, file1

pro model_add_sepsurf, nskip
  common shared_var

  print, 'Adding Separatrix surface rings'
  
  rings = read_rings(filename, breaks=breaks, nskip=nskip)

  foreach i, inull do begin

    if null[i].sign gt 0 then colour = [255,128,128] else begin
      if null[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
    endelse

    for iring = 0, n_elements(rings[i])-1 do begin
      brks = [-1, where(breaks[i, iring] eq 1, /null), n_elements(breaks[i, iring])-1]
      for ib = 0, n_elements(brks)-2 do begin
        if brks[ib] ne brks[ib+1] then oModel -> add, obj_new("IDLgrPolyline", rings[i, iring, 0, brks[ib]+1:brks[ib+1]], rings[i, iring , 1, brks[ib]+1:brks[ib+1]], rings[i, iring, 2, brks[ib]+1:brks[ib+1]], color=colour)
      endfor
    endfor
    
  endforeach

end

pro model_add_spines, nskip
  common shared_var

  print, 'Adding Spines'

  spines = read_spines(filename)

  foreach i, inull do begin
    if null[i].sign gt 0 then col = [250,0,0] else col = [0,0,250]
    for dir = 0, 1 do begin
      
      oModel -> add, obj_new("IDLgrPolyline", spines[i, dir, 0, 0:-1:nskip], spines[i, dir, 1, 0:-1:nskip], spines[i, dir, 2, 0:-1:nskip], color=col, thick=4)

    endfor
  endforeach
end

pro model_add_fanlines, nskip
  common shared_var

  print, 'Adding Separatrix surface field lines'
  mxline = 2000

  rings = read_rings(filename, nskip=nskip)

  foreach i, inull do begin
    if null[i].sign gt 0 then colour = [255,128,128] else begin
      if null[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
    endelse

    ring = rings[i, 10, *, *]
    n = n_elements(ring)/3

    if n le 150 then skip = 2 else skip = fix((n-2)/150.)

    xf = ring[0,*]
    yf = ring[1,*]
    zf = ring[2,*]

    for k = 0, n-1, n/20 do begin
      print, k
      startpt = [xf[k],yf[k],zf[k]]

      h = 0.01*ds
      hmin = 0.001*ds
      hmax = 0.5*ds
      epsilon = 1.0d-5

      line = fieldline3d(startpt,bgrid,xx,yy,zz,h,hmin,hmax,epsilon)

      dist = (line[0,*]-null[i].pos[0])^2+(line[1,*]-null[i].pos[1])^2+(line[2,*]-null[i].pos[2])^2
      imin = where(dist eq min(dist))
      print, imin, n_elements(line)
      if null[i].sign lt 0 then begin
        minloc = min(imin)
        line = line[*, 0:minloc]
      endif else begin
        minloc = max(imin)
        line = line[*, minloc:-1]
      endelse
      
      oModel -> add, obj_new("IDLgrPolyline",line[0,*],line[1,*],line[2,*],color=colour,thick=2)
    endfor
  endforeach

end

pro model_add_separators, nskip
  common shared_var

  print, 'Adding separators'

  seps = read_separators(filename)

  foreach i, inull do begin
    
    for isep = 0, n_elements(seps[i])-1 do begin
      oModel -> add,obj_new("IDLgrPolyline", seps[i, isep, 0, 0:-1:nskip], seps[i, isep, 1, 0:-1:nskip], seps[i, isep, 2, 0:-1:nskip], color=[0,160,60], thick=4)
    endfor

  endforeach
end

pro model_add_nulls
  common shared_var

  radius = 10*ds
  
  foreach i, inull do begin
    mesh_obj, 4, vert, poly, replicate(radius,21,21)
    for j = 0, 2 do vert[j,*] = vert[j,*] + null[i].pos[j]
    if null[i].sign gt 0 then colour = [255,0,0] else begin
      if null[i].sign lt 0 then colour = [0,0,255] else colour = [0,255,0]
    endelse
    oModel -> add, obj_new("IDLgrPolygon",data=vert,polygons=poly,color=colour,/shading)
  endforeach
  
end

pro model_add_box
  common shared_var
  
  print, "Adding box"
  box = dblarr(3,2)
  box[*,0] = [min(xx),min(yy),min(zz)] 
  box[*,1] = [max(xx),max(yy),max(zz)]
  line = double(transpose([[0,0,0],[1,0,0],[1,0,1],[1,0,0],[1,1,0],[1,1,1],[1,1,0], $
    [0,1,0],[0,1,1],[0,1,0],[0,0,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[0,0,1]]))
  
  line[*,0] = line[*,0]*(box[0,1] - box[0,0]) + box[0,0]
  line[*,1] = line[*,1]*(box[1,1] - box[1,0]) + box[1,0]
  line[*,2] = line[*,2]*(box[2,1] - box[2,0]) + box[2,0]

  dist = min([n_elements(xx), n_elements(yy), n_elements(zz)])*ds/5
  oModel -> add, obj_new('idlgrtext', 'x', locations=[box[0,0]+dist,box[1,0],box[2,0]], /onglass)
  oModel -> add, obj_new('idlgrtext', 'y', locations=[box[0,0],box[1,0]+dist,box[2,0]], /onglass)
  oModel -> add, obj_new('idlgrtext', 'z', locations=[box[0,0],box[1,0],box[2,0]+dist], /onglass)
  
  oModel -> add, obj_new('IDLgrPolyline',line[*,0],line[*,1],line[*,2],color=[0,0,0])
end

function mk_model_mag_sep,fname,nulls=nulls,separators=separators,$
sepsurf=sepsurf,spines=spines,box=box,fanlines=fanlines,nskip=nskip
  common shared_var

  filename = fname

  get_lun, lun
  
  nx = 0L
  ny = 0L
  nz = 0L
  openr, lun, fname
  readu, lun, nx, ny, nz
  
  xx = dblarr(nx,ny,nz)
  yy = dblarr(nx,ny,nz)
  zz = dblarr(nx,ny,nz)
  readu, lun, xx, yy, zz
  bgrid = dblarr(nx,ny,nz,3)
  bgrid[*,*,*,0] = xx
  bgrid[*,*,*,1] = yy
  bgrid[*,*,*,2] = zz

  xx = dblarr(nx)
  yy = dblarr(ny)
  zz = dblarr(nz)
  readu, lun, xx, yy, zz
  close, lun
  free_lun, lun

  ds = min([min(xx[1:-1]-xx[0:-2]), min(yy[1:-1]-yy[0:-2]), min(zz[1:-1]-zz[0:-2])])

  file1 = strmid(fname, strlen('data/'), strlen(fname))
  file1 = strmid(file1, 0, strlen(file1)-strlen('.dat'))
  
  null = read_nulls(filename)
  inull = indgen(n_elements(null))

  oModel = obj_new("IDLgrModel")
  if keyword_set(box)        then model_add_box
  if keyword_set(nulls)      then model_add_nulls
  if keyword_set(fanlines)   then model_add_fanlines, nskip
  if keyword_set(spines)     then model_add_spines, nskip
  if keyword_set(sepsurf)    then model_add_sepsurf, nskip
  if keyword_set(separators) then model_add_separators, nskip
  
  return, oModel
end
