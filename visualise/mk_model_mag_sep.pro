pro model_add_sepsurf,oModel,frame
  
  print, 'Adding Separatrix surface rings'
  nulls = getnulls()

  for i = 0, n_elements(nulls)-1 do begin
    if nulls[i].sign gt 0 then colour = [255,128,128] else begin
      if nulls[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
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
    
  endfor
  
  free_lun,fan
end

pro model_add_spines, oModel, frame, bgrid, xx, yy, zz

  print, 'Adding Spines'
  dx = abs(xx[1]-xx[0])
  ro = 0.1*dx

  nulls = getnulls()
  for i = 0, n_elements(nulls)-1 do begin
    if nulls[i].sign gt 0 then col = [250,0,0] else col = [0,0,250]
    foreach dir, [1, -1] do begin
      
      h = 0.01*dx
      hmin = 0.001*dx
      hmax = 0.5*dx
      epsilon = 1.0d-5

      startpt = nulls[i].gridpos + dir*ro*nulls[i].spine
      print, "Starting spine line at", startpt
      line = fieldline3d(startpt,bgrid,xx,yy,zz,nulls[i].sign*h,hmin,hmax,epsilon,/oneway)

      oModel -> add, obj_new("IDLgrPolyline",line[0,*],line[1,*],line[2,*],color=col,thick=4)
    endforeach
  endfor
end

pro model_add_fanlines, oModel, frame, bgrid, xx, yy, zz

  print, 'Adding Separatrix surface field lines'
  dx = abs(xx[1]-xx[0])
  mxline = 2000
  ro = 0.1*dx   

  nulls = getnulls()

  pp = where(nulls.sign ge 0, npp)
  nn = where(nulls.sign le 0, nnn)

  for i = 0, n_elements(nulls)-1 do begin
    if nulls[i].sign gt 0 then colour = [255,128,128] else begin
      if nulls[i].sign lt 0 then colour = [128,128,255] else colour = [128,255,128]
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

      h = 0.01*dx
      hmin = 0.001*dx
      hmax = 0.5*dx
      epsilon = 1.0d-5

      line = fieldline3d(startpt,bgrid,xx,yy,zz,h,hmin,hmax,epsilon)

      !null = min((line[0,*]-nulls[i].gridpos[0])^2+(line[1,*]-nulls[i].gridpos[1])^2+(line[2,*]-nulls[i].gridpos[2])^2, minloc)
      if nulls[i].sign gt 0 then line = line[*, 0:minloc] else line = line[*, minloc:-1]

      oModel -> add, obj_new("IDLgrPolyline",line[0,*],line[1,*],line[2,*],color=colour,thick=2)
    endfor
  endfor

  free_lun,fan
end

pro model_add_separators,oModel,frame

  print, 'Adding separators'
   
  nulls = getnulls()

  for null = 1, n_elements(nulls) do begin
    seps = read_separators('output/sep' + string(null,'(I4.4)') + '.dat')
    
    if seps ne !null then begin
      for j = 0, (size(seps))[1]-1 do begin
        col = [0,160,60] 
        oModel -> add,obj_new("IDLgrPolyline",seps[j,*,0],seps[j,*,1],seps[j,*,2],color=col,thick=4)
      endfor
    endif
  endfor
end

pro model_add_nulls,oModel,frame,xx,yy,zz

  radius = min([xx[-1]-xx[0],yy[-1]-yy[0],zz[-1]-zz[0]])/40
  nulls = getnulls()
  
  for i = 0, n_elements(nulls)-1 do begin
    mesh_obj, 4, vert, poly, replicate(radius,21,21)
    for j = 0, 2 do vert[j,*] = vert[j,*] + nulls[i].gridpos[j]
    if nulls[i].sign gt 0 then colour = [255,0,0] else begin
      if nulls[i].sign lt 0 then colour = [0,0,255] else colour = [0,255,0]
    endelse
    oModel -> add, obj_new("IDLgrPolygon",data=vert,polygons=poly,color=colour,/shading)
  endfor
  
end

pro model_add_box,oModel,box
  
  print, "Adding box"
  line = transpose([[0,0,0],[1,0,0],[1,0,1],[1,0,0],[1,1,0],[1,1,1],[1,1,0], $
    [0,1,0],[0,1,1],[0,1,0],[0,0,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[0,0,1.]])
  
  line[*,0] = line[*,0]*(box[0,1] - box[0,0]) + box[0,0]
  line[*,1] = line[*,1]*(box[1,1] - box[1,0]) + box[1,0]
  line[*,2] = line[*,2]*(box[2,1] - box[2,0]) + box[2,0]
  
  oModel -> add,obj_new('IDLgrPolyline',line[*,0],line[*,1],line[*,2], $
    color=[0,0,0])
end

function mk_model_mag_sep,fname,nulls=nulls,separators=separators,$
sepsurf=sepsurf,spines=spines,box=box,fanlines=fanlines

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
  
  print,'NOTE: HARDWIRED TO USE GRIDCELL COORDINATES'
  xx = dindgen(nx)+1
  yy = dindgen(ny)+1
  zz = dindgen(nz)+1
  
  box = dblarr(3,2)
  box[*,0] = [min(xx),min(yy),min(zz)] 
  box[*,1] = [max(xx),max(yy),max(zz)]
  
  print, nx, ny, nz

  oModel = obj_new("IDLgrModel")
  if keyword_set(box)        then model_add_box,oModel,box
  if keyword_set(nulls)      then model_add_nulls,oModel,frame,xx,yy,zz
  if keyword_set(fanlines)   then model_add_fanlines,oModel,frame,bgrid,xx,yy,zz
  if keyword_set(spines)     then model_add_spines,oModel,frame,bgrid,xx,yy,zz
  if keyword_set(sepsurf)    then model_add_sepsurf,oModel,frame
  if keyword_set(separators) then model_add_separators,oModel,frame
  
  return, oModel
end
