common shared_var_model3d, bgrid, xx, yy, zz, oModel, nulldata, ds, filename, nulllist, csystem_model3d

function sphr2cart, pts
  
  carts = pts
  carts[0, *] = pts[0, *]*sin(pts[1, *])*cos(pts[2, *])
  carts[1, *] = pts[0, *]*sin(pts[1, *])*sin(pts[2, *])
  carts[2, *] = pts[0, *]*cos(pts[1, *])

  return, carts

end

pro model_add_sepsurf, nskip
  common shared_var_model3d

  print, 'Adding Separatrix surface rings'
  
  rings = read_rings(filename, breaks=breaks, nskip=nskip, null_list=nulllist)

  foreach inull, nulllist-1 do begin

    if nulldata[inull].sign gt 0 then colour = [255, 128, 128] else begin
      if nulldata[inull].sign lt 0 then colour = [128, 128, 255] else colour = [128, 255, 128]
    endelse

    for iring = 0, n_elements(rings[inull])-1 do begin
      brks = [-1, where(breaks[inull, iring] eq 1, /null), n_elements(breaks[inull, iring])-1]
      ring = rings[inull, iring, *, *]
      if csystem_model3d eq 'spherical' then ring = sphr2cart(ring)
      for ib = 0, n_elements(brks)-2 do begin
        if brks[ib] ne brks[ib+1] then oModel.add, obj_new("IDLgrPolyline", ring[0, brks[ib]+1:brks[ib+1]], ring[1, brks[ib]+1:brks[ib+1]], ring[2, brks[ib]+1:brks[ib+1]], color=colour)
      endfor
    endfor
    
  endforeach

end

pro model_add_hcs, nskip
  common shared_var_model3d

  rings = read_rings(filename, breaks=breaks, nskip=nskip, /hcs)

  for inull = 0, n_elements(rings)-1 do begin

    colour = [0, 255, 0]

    for iring = 0, n_elements(rings[inull])-1 do begin
      brks = [-1, where(breaks[inull, iring] eq 1, /null), n_elements(breaks[inull, iring])-1]
      ring = rings[inull, iring, *, *]
      if csystem_model3d eq 'spherical' then ring = sphr2cart(ring)
      for ib = 0, n_elements(brks)-2 do begin
        if brks[ib] ne brks[ib+1] then oModel.add, obj_new("IDLgrPolyline", ring[0, brks[ib]+1:brks[ib+1]], ring[1, brks[ib]+1:brks[ib+1]], ring[2, brks[ib]+1:brks[ib+1]], color=colour)
      endfor
    endfor
    
  endfor

end

pro model_add_spines
  common shared_var_model3d

  print, 'Adding Spines'

  spines = read_spines(filename, null_list=nulllist)

  foreach inull, nulllist-1 do begin
    if nulldata[inull].sign gt 0 then col = [250, 0, 0] else col = [0, 0, 250]
    for dir = 0, 1 do begin
      
      spine = spines[inull, dir, *, *]
      if csystem_model3d eq 'spherical' then spine = sphr2cart(spine)
      oModel.add, obj_new("IDLgrPolyline", spine[0, *], spine[1, *], spine[2, *], color=col, thick=4)

    endfor
  endforeach
end

pro model_add_fanlines, nlines, nring=nring
  common shared_var_model3d

  print, 'Adding Separatrix surface field lines'

  rings = read_rings(filename, nskip=nskip, null_list=nulllist)

  foreach inull, nulllist-1 do begin
    if nulldata[inull].sign gt 0 then colour = [255, 128, 128] else begin
      if nulldata[inull].sign lt 0 then colour = [128, 128, 255] else colour = [128, 255, 128]
    endelse

    if not keyword_set(nring) then begin
      ring = rings[inull, n_elements(rings[inull])/5]
    endif else begin
      ring = rings[inull, nring]
    endelse

    n = n_elements(ring)/3/nlines

    xf = ring[0, *]
    yf = ring[1, *]
    zf = ring[2, *]

    for k = 0, n_elements(ring)/3 - 1, n do begin
      startpt = [xf[k], yf[k], zf[k]]

      h = 0.01
      hmin = 0.001
      hmax = 0.5
      epsilon = 1.0d-5

      line = fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem_model3d)

      dist = (line[0, *] - nulldata[inull].pos[0])^2 + $
        (line[1, *] - nulldata[inull].pos[1])^2 + $
        (line[2, *] - nulldata[inull].pos[2])^2
      imin = where(dist eq min(dist))
      if nulldata[inull].sign lt 0 then begin
        line = line[*, 0:min(imin)]
      endif else begin
        line = line[*, max(imin):-1]
      endelse

      if csystem_model3d eq 'spherical' then line = sphr2cart(line)
      
      oModel.add, obj_new("IDLgrPolyline", line[0, *], line[1, *], line[2, *], color=colour, thick=2)
    endfor
  endforeach

end

pro model_add_separators, hcs=hcs
  common shared_var_model3d

  print, 'Adding separators'

  if hcs then colour = [255, 165, 0] else colour = [255, 255, 0]

  seps = read_separators(filename, null_list=nulllist, hcs=hcs)

  if keyword_set(hcs) then to_do = [0] else to_do = nulllist-1

  foreach inull, to_do do begin
    
    for isep = 0, n_elements(seps[inull])-1 do begin
      sep = seps[inull, isep, *, *]
      if csystem_model3d eq 'spherical' then sep = sphr2cart(sep)
      oModel.add, obj_new("IDLgrPolyline", sep[0, *], sep[1, *], sep[2, *], color=colour, thick=4)
    endfor

  endforeach
end

pro model_add_nulls, size=size
  common shared_var_model3d

  boxsize = min([xx[-1] - xx[0], yy[-1] - yy[0], zz[-1] - zz[0]])/40
  if not keyword_set(size) then size = 1
  r = max([boxsize, ds])*size
  radius = ds
  
  foreach inull, nulllist-1 do begin
    mesh_obj, 4, vert, poly, replicate(radius,21,21)
    pos = nulldata[inull].pos
    if csystem_model3d eq 'spherical' then pos = [pos[0]*sin(pos[1])*cos(pos[2]), pos[0]*sin(pos[1])*sin(pos[2]), pos[0]*cos(pos[1])]
    for j = 0, 2 do vert[j, *] = vert[j, *] + pos[j]
    if nulldata[inull].sign gt 0 then colour = [255, 0, 0] else begin
      if nulldata[inull].sign lt 0 then colour = [0, 0, 255] else colour = [0, 255, 0]
    endelse
    oModel.add, obj_new("IDLgrPolygon", data=vert, polygons=poly, color=colour, /shading)
  endforeach
  
end

pro model_add_box
  common shared_var_model3d
  
  print, "Adding box"
  box = dblarr(3,2)
  box[*,0] = [min(xx), min(yy), min(zz)] 
  box[*,1] = [max(xx), max(yy), max(zz)]
  line = double(transpose([[0, 0, 0], [1, 0, 0], [1, 0, 1], $
      [1, 0, 0], [1, 1, 0], [1, 1, 1], $
      [1, 1, 0], [0, 1, 0], [0, 1, 1], $
      [0, 1, 0], [0, 0, 0], [0, 0, 1], $
      [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1]]))
  
  line[*, 0] = line[*, 0]*(box[0, 1] - box[0, 0]) + box[0, 0]
  line[*, 1] = line[*, 1]*(box[1, 1] - box[1, 0]) + box[1, 0]
  line[*, 2] = line[*, 2]*(box[2, 1] - box[2, 0]) + box[2, 0]

  dist = min([n_elements(xx), n_elements(yy), n_elements(zz)])*ds/5
  oModel.add, obj_new('idlgrtext', 'x', locations=[box[0, 0]+dist, box[1, 0], box[2, 0]], /onglass)
  oModel.add, obj_new('idlgrtext', 'y', locations=[box[0, 0], box[1, 0]+dist, box[2, 0]], /onglass)
  oModel.add, obj_new('idlgrtext', 'z', locations=[box[0, 0], box[1, 0], box[2, 0]+dist], /onglass)
  
  oModel.add, obj_new('IDLgrPolyline', line[*, 0], line[*, 1], line[*, 2], color=[0, 0, 0])
end

pro save, dimensions=dimensions

  xobjview_write_image, 'figures/' + rd.prefix(filename) + '-model3d.png', 'png', dimensions=dimensions

end

pro set_null_list, lst
  common shared_var_model3d

  if lst eq !null then begin
    nulllist = nulldata.number
  endif else begin
    if max(lst) gt max(nulldata.number) or min(lst) < 1 then print, 'Invalid list'
    nulllist = lst
  endelse
end

function mk_model, fname, nulls=nulls, separators=separators, sepsurf=sepsurf, hcs=hcs, hcs_sep=hcs_sep, spines=spines, box=box, fanlines=fanlines, nskip=nskip, null_list=null_list, coordsystem=coordsystem
  common shared_var_model3d

  if not keyword_set(coordsystem) then coordsystem = 'cartesian'
  csystem_model3d = coordsystem

  filename = fname

  field = read_field(fname)
  nx = n_elements(field.x)
  ny = n_elements(field.y)
  nz = n_elements(field.z)
  bgrid = dblarr(nx, ny, nz, 3)
  bgrid[*, *, *, 0] = field.bx
  bgrid[*, *, *, 1] = field.by
  bgrid[*, *, *, 2] = field.bz
  xx = field.x
  yy = field.y
  zz = field.z
  field = !null

  ds = min([min(xx[1:-1] - xx[0:-2]), min(yy[1:-1] - yy[0:-2]), min(zz[1:-1] - zz[0:-2])])

  nulldata = read_nulls(filename)
  set_null_list, null_list

  oModel = obj_new("IDLgrModel")
  if keyword_set(box)        then model_add_box
  if keyword_set(nulls)      then model_add_nulls
  if keyword_set(fanlines)   then model_add_fanlines, nskip
  if keyword_set(spines)     then model_add_spines
  if keyword_set(sepsurf)    then model_add_sepsurf, nskip
  if keyword_set(hcs)        then model_add_hcs, nskip
  if keyword_set(hcs_sep)    then model_add_separators, /hcs
  if keyword_set(separators) then model_add_separators
  
  return, oModel
end
