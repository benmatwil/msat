pro model_add_sepsurf,oModel,frame,rake=rake
  
  print,'Separatrix surface rings'
  get_nulls,nnulls,signs,r,spine,nulls
;nnulls=1
  for i=0,nnulls-1 do begin
    if keyword_set(rake) then colour=[0,255,0] $
    else begin
      if signs[i] gt 0 then colour = [255,128,128] else begin
        if signs[i] lt 0 then colour = [128,128,255] else colour = [128,255,128]
      endelse
    endelse
    
    file='output/ringidl'+string(i+1,'(I4.4)')+'.dat'
    print,file
    
    get_lun,fan
    openr,fan,file,/f77_unformatted
    nrings = -1L
    readu,fan,nrings
    print, 'number of rings', nrings
    n=-1L
    readu,fan,n
    
    while n ge 1 do begin
      x=dblarr(n)
      y=dblarr(n)
      z=dblarr(n)
      readu,fan,x,y,z

      break = lonarr(n)
      readu,fan,break

      breaks = [-1, where(break eq 1, /null), n-1]
      for ib = 0, n_elements(breaks)-2 do begin
        ;print, breaks[ib]+1,breaks[ib+1]
        if breaks[ib] ne breaks[ib+1] then oModel -> add, obj_new("IDLgrPolyline", x[breaks[ib]+1:breaks[ib+1]], y[breaks[ib]+1:breaks[ib+1]], z[breaks[ib]+1:breaks[ib+1]], color=colour)
      endfor

      ;if where(break eq 1, /null) ne !null then stop
          
      readu,fan,n
    endwhile
     
    close,fan
    
  endfor
  
  free_lun,fan
end

pro model_add_spines,oModel,frame,bgrid,xx,yy,zz
  print,'Spines'
  dx = abs(xx[1]-xx[0])
  h = 0.01*dx
  hmin = 0.001*dx
  hmax = 0.5*dx
  epsilon = 1.0d-5
  mxline = 2000
  ro = 0.1*dx

  get_nulls,nnulls,signs,r,spines,nulls
  for i=0,nnulls-1 do begin
    spine = (spines(*,i))
    posi = r(*,i)
    print,spine
    ds = signs(i)/abs(signs(i))
    if (ds gt 0) then col = [250,0,0] else col = [0,0,250]
  startpt = posi+ro*spine
  print,startpt
  linef = tracefield_3D(startpt,bgrid,xx,yy,zz,ds*h,hmin,hmax,epsilon,mxline)
  ii = where(linef[0,*] gt -8000)
  if (ii[0] gt -1) then oModel->add,obj_new("IDLgrPolyline", $
      linef[0,ii],linef[1,ii],linef[2,ii],color=col,thick=4)
  startpt = posi-ro*spine
  print,startpt
  linef = tracefield_3D(startpt,bgrid,xx,yy,zz,ds*h,hmin,hmax,epsilon,mxline)
  ii = where(linef[0,*] gt -8000)
  if (ii[0] gt -1) then oModel->add,obj_new("IDLgrPolyline",linef[0,ii],linef[1,ii],linef[2,ii],color=col,thick=4)
  endfor
end

pro model_add_fanlines,oModel,frame,bgrid,xx,yy,zz

  print,'Separatrix surface field lines'
  dx = abs(xx[1]-xx[0])
  h = 0.01*dx
  hmin = 0.001*dx
  hmax = 0.5*dx
  epsilon = 1.0d-5
  mxline = 2000
  ro = 0.1*dx   

  get_nulls,nnulls,signs,r,spine,nulls

  pp = where(signs ge 0,npp)
  nn = where(signs le 0,nnn)

  for i=0,nnulls-1 do begin
    if keyword_set(rake) then colour=[0,255,0] else begin
      if signs[i] gt 0 then colour = [255,128,128] else begin
        if signs[i] lt 0 then colour = [128,128,255] else colour = [128,255,128]
      endelse
    endelse

    file='output/ringidl'+string(i+1,'(I4.4)')+'.dat'
    print,file
    
    get_lun,fan
    openr,fan,file,/f77_unformatted
    nrings = -1L
    readu,fan,nrings

    n=-1L
    
    fring = nrings/5
    print,i,nrings,fring
    for iring=1,fring do begin
      readu,fan,n
      print, 'n', n
      xf=dblarr(n)
      yf=dblarr(n)
      zf=dblarr(n)
      brk=lonarr(n)
      readu,fan,xf,yf,zf,brk
     endfor      
     close,fan

    if (n le 150) then skip = 1 else skip = fix((n-2)/150.)
    print,i,n, skip
;    help,xf,yf,xf
    for k=0,n-2,skip do begin
       startpt = [xf[k],yf[k],zf[k]]
       linef = tracefield_3D(startpt,bgrid,xx,yy,zz,h,hmin,hmax,epsilon,mxline)
           jj = where(linef[0,*] gt -8000)
           linef = linef[*,jj] 
           if (signs(i) gt 0) then begin
             dist = sqrt((linef[0,*]-r[0,i])^2.+(linef[1,*]-r[1,i])^2.+(linef[2,*]-r[2,i])^2.)
             ll = where(abs(dist) eq min(abs(dist)))
             linef = linef[*,0:ll[0]] 
           endif
           oModel->add,obj_new("IDLgrPolyline",linef[0,*],linef[1,*], $
                linef[2,*],color=colour,thick=2)

           linef = tracefield_3D(startpt,bgrid,xx,yy,zz,-h,hmin,hmax,epsilon,mxline)
           jj = where(linef[0,*] gt -8000)
           linef = linef[*,jj] 
           if (signs(i) lt 0) then begin
             dist = sqrt((linef[0,*]-r[0,i])^2.+(linef[1,*]-r[1,i])^2.+(linef[2,*]-r[2,i])^2.)
             ll = where(abs(dist) eq min(abs(dist)))
             linef = linef[*,0:ll[0]] 
           endif
           oModel->add,obj_new("IDLgrPolyline",linef[0,*],linef[1,*], $
                linef[2,*],color=colour,thick=2)
      endfor
 endfor

  free_lun,fan
end

pro model_add_separators,oModel,frame

  print,'plot separators'
   
    get_nulls,nnulls,signs,r,spine,nulls
;    nnulls=1
    for null=1,nnulls do begin
;    null = 1
    seps=read_separators('output/sep'+string(null,'(I4.4)')+'.dat')
    
    print, null, ' output/sep'+string(null,'(I4.4)')+'.dat', size(seps)
    if keyword_set(seps) then begin
    for j=0,(size(seps))[1]-1 do begin
;    j=4
    col = [0,160,60] 
;    if j eq 0 then col = [80,0,50] else col = [50,160,0]
;    if j eq 1 then col = [0,250,250]
;    if j eq 2 then col = [250,0,250]
;    if j eq 3 then col = [250,250,0]
;    if j eq 4 then col = [0,160,60] 
      oModel->add,obj_new("IDLgrPolyline",seps[j,*,0],seps[j,*,1],seps[j,*,2],$
        color=col,thick=4)  
;      oModel->add,obj_new("IDLgrPolyline",seps[0,*,0],seps[0,*,1],seps[0,*,2],$
;        color=[255,0,0],thick=2)  
;      oModel->add,obj_new("IDLgrPolyline",seps[1,*,0],seps[1,*,1],seps[1,*,2],$
;        color=[0,255,0],thick=2)  
;      oModel->add,obj_new("IDLgrPolyline",seps[2,*,0],seps[2,*,1],seps[2,*,2],$
;        color=[0,0,255],thick=2)  
    endfor
    endif
  endfor
end

pro model_add_nulls,oModel,frame

  radius=0.1
  scl=((8*!dpi)/(3*sqrt(6)))^(1/3.)
  print,"scl=",scl
  
  get_nulls,nnulls,signs,r,spine,nulls
  
  print,nnulls
  print,signs
  print,r
  
  if keyword_set(nulls) eq 0 then return
  for i=0,nnulls-1 do begin
    if abs(signs(i)) eq 2 then begin
      mesh_obj,6,vert,poly,[[0,0,sqrt(3)-sqrt(2)],[1,0,sqrt(3)-sqrt(2)], $
        [0,0,sqrt(3)]]*radius*scl
    endif else begin
      mesh_obj,4,vert,poly,replicate(radius,21,21)
    endelse
    for j=0,2 do vert[j,*]=vert[j,*]+r(j,i)
    colour=(signs(i) gt 0)? [255,0,0] : (signs(i) lt 0)? [0,0,255] $
      : [0,255,0]
    oModel->add,obj_new("IDLgrPolygon",data=vert,polygons=poly,color=colour, $
      /shading)
  endfor
  
end

pro model_add_box,oModel,box
  ; add the box
  line=transpose([[0,0,0],[1,0,0],[1,0,1],[1,0,0],[1,1,0],[1,1,1],[1,1,0], $
    [0,1,0],[0,1,1],[0,1,0],[0,0,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[0,0,1.]])
  
  line[*,0]=line[*,0]*(box[0,1]-box[0,0])+box[0,0]
  line[*,1]=line[*,1]*(box[1,1]-box[1,0])+box[1,0]
  line[*,2]=line[*,2]*(box[2,1]-box[2,0])+box[2,0]
  
  print,'line', max(line), min(line)
  oModel->add,obj_new('IDLgrPolyline',line[*,0],line[*,1],line[*,2], $
    color=[0,0,0])
end

function mk_model_mag_sep,fname,nulls=nulls,separators=separators,$
sepsurf=sepsurf,spines=spines,box=box,fanlines=fanlines

   get_lun,lun
   
   nx=0l
   ny=0l
   nz=0l
   
   openr,lun,fname
   readu,lun,nx,ny,nz
   
   
   xx = dblarr(nx,ny,nz)
   yy = dblarr(nx,ny,nz)
   zz = dblarr(nx,ny,nz)
   
   readu,lun,xx,yy,zz
   bgrid = dblarr(nx,ny,nz,3)
   bgrid[*,*,*,0] = xx
   bgrid[*,*,*,1] = yy
   bgrid[*,*,*,2] = zz
   xx = dblarr(nx)
   yy = dblarr(ny)
   zz = dblarr(nz)
  
   readu,lun,xx,yy,zz
   close,lun
   free_lun,lun
   
   
   print,'NOTE: HARDWIRED TO USE GRIDCELL COORDINATES'
   xx=dindgen(nx)+1
   yy=dindgen(ny)+1
   zz=dindgen(nz)+1
   
   
   box = dblarr(3,2)
   box[*,0] = [min(xx),min(yy),min(zz)] 
   box[*,1] = [max(xx),max(yy),max(zz)]
   
   print,nx,ny,nz
   print,max(xx),min(xx)
   print,max(yy), min(yy)
   print,max(zz),min(zz)

  oModel=obj_new("IDLgrModel")
  if keyword_set(box)        then model_add_box,oModel,box
  if keyword_set(nulls)      then model_add_nulls,oModel,frame
  if keyword_set(fanlines)   then model_add_fanlines,oModel,frame,bgrid,xx,yy,zz
  if keyword_set(spines)     then model_add_spines,oModel,frame,bgrid,xx,yy,zz
  if keyword_set(sepsurf)    then model_add_sepsurf,oModel,frame
  if keyword_set(separators) then model_add_separators,oModel,frame
  
  return,oModel
end
