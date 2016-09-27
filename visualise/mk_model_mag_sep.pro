common shared_var, bgrid, xx, yy, zz, oModel, frame, null, ds

pro model_add_sepsurf
  common shared_var

  print, 'Adding Separatrix surface rings'

  for i = 0, n_elements(null)-1 do begin
    if null[i].sign gt 0 then colour = [255,128,128] else begin
      if null[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
    endelse
    
    file = 'output/ringidl' + string(i+1,'(I4.4)') + '.dat'
    print, file
    
    get_lun, fan
    openr, fan, file, /f77_unformatted
    nrings = -1L
    readu, fan, nrings
    print, 'number of rings', nrings
    n = -1L
    readu, fan, n
    
    while n ge 1 do begin
      x = dblarr(n)
      y = dblarr(n)
      z = dblarr(n)
      readu, fan, x, y, z

      break = lonarr(n)
      readu, fan, break

      breaks = [-1, where(break eq 1, /null), n-1]
      for ib = 0, n_elements(breaks)-2 do begin
        if breaks[ib] ne breaks[ib+1] then oModel -> add, obj_new("IDLgrPolyline", x[breaks[ib]+1:breaks[ib+1]], y[breaks[ib]+1:breaks[ib+1]], z[breaks[ib]+1:breaks[ib+1]], color=colour)
      endfor

      readu, fan, n
    endwhile
     
    close,fan
    free_lun,fan
    
  endfor

end

pro model_add_spines
  common shared_var

  print, 'Adding Spines'
  ro = 1d-2*ds

  for i = 0, n_elements(null)-1 do begin
    if null[i].sign gt 0 then col = [250,0,0] else col = [0,0,250]
    foreach dir, [1, -1] do begin
      
      h = 0.01*ds
      hmin = 0.001*ds
      hmax = 0.5*ds
      epsilon = 1.0d-5

      startpt = null[i].pos + dir*ro*null[i].spine
      print, "Starting spine line at", startpt
      line = fieldline3d(startpt,bgrid,xx,yy,zz,null[i].sign*h,hmin,hmax,epsilon,/oneway)

      oModel -> add, obj_new("IDLgrPolyline",line[0,*],line[1,*],line[2,*],color=col,thick=4)
    endforeach
  endfor
end

pro model_add_fanlines
  common shared_var

  print, 'Adding Separatrix surface field lines'
  mxline = 2000
  ro = 0.1*ds

  pp = where(null.sign ge 0, npp)
  nn = where(null.sign le 0, nnn)

  for i = 0, n_elements(null)-1 do begin
    if null[i].sign gt 0 then colour = [255,128,128] else begin
      if null[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
    endelse

    file = 'output/ringidl' + string(i+1,'(I4.4)') + '.dat'
    print, file
    
    get_lun, fan
    openr, fan, file, /f77_unformatted
    nrings = -1L
    readu, fan, nrings

    n = -1L
   
    for iring = 1, nrings/5 do begin
      readu, fan, n
      xf = dblarr(n)
      yf = dblarr(n)
      zf = dblarr(n)
      brk = lonarr(n)
      readu, fan, xf, yf, zf
      readu, fan, brk
    endfor      
    close, fan

    if n le 150 then skip = 2 else skip = fix((n-2)/150.)

    for k = 0, n-2, skip do begin
      startpt = [xf[k],yf[k],zf[k]]

      h = 0.01*ds
      hmin = 0.001*ds
      hmax = 0.5*ds
      epsilon = 1.0d-5

      line = fieldline3d(startpt,bgrid,xx,yy,zz,h,hmin,hmax,epsilon)
     
      dist = (line[0,*]-null[i].pos[0])^2+(line[1,*]-null[i].pos[1])^2+(line[2,*]-null[i].pos[2])^2
      imin = where(dist eq min(dist))
      if null[i].sign gt 0 then begin
        minloc = min(imin)
        line = line[*, 0:minloc]
      endif else begin
        minloc = max(imin)
        line = line[*, minloc:-1]
      endelse
      
      oModel -> add, obj_new("IDLgrPolyline",line[0,*],line[1,*],line[2,*],color=colour,thick=2)
    endfor
  endfor

  free_lun,fan
end

pro model_add_separators
  common shared_var

  print, 'Adding separators'

  for null = 0, n_elements(null)-1 do begin
    seps = read_separators('output/sep' + string(null,'(I4.4)') + '.dat')
    
    if seps ne !null then begin
      for j = 0, (size(seps))[1]-1 do begin
        col = [0,160,60] 
        oModel -> add,obj_new("IDLgrPolyline",seps[j,*,0],seps[j,*,1],seps[j,*,2],color=col,thick=4)
      endfor
    endif
  endfor
end

pro model_add_nulls
  common shared_var

  radius = min([xx[-1]-xx[0],yy[-1]-yy[0],zz[-1]-zz[0]])/200
  
  for i = 0, n_elements(null)-1 do begin
    mesh_obj, 4, vert, poly, replicate(radius,21,21)
    for j = 0, 2 do vert[j,*] = vert[j,*] + null[i].pos[j]
    if null[i].sign gt 0 then colour = [255,0,0] else begin
      if null[i].sign lt 0 then colour = [0,0,255] else colour = [0,255,0]
    endelse
    oModel -> add, obj_new("IDLgrPolygon",data=vert,polygons=poly,color=colour,/shading)
  endfor
  
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
sepsurf=sepsurf,spines=spines,box=box,fanlines=fanlines
  common shared_var

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
  
  null = read_nulls()

  oModel = obj_new("IDLgrModel")
  if keyword_set(box)        then model_add_box
  if keyword_set(nulls)      then model_add_nulls
  if keyword_set(fanlines)   then model_add_fanlines
  if keyword_set(spines)     then model_add_spines
  if keyword_set(sepsurf)    then model_add_sepsurf
  if keyword_set(separators) then model_add_separators
  
  return, oModel
end
