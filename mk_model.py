pro model_add_sepsurf,oModel,frame,rake=rake
  
  print,'Separatrix surface rings'
  get_nulls,nnulls,signs,r,spine,nulls
  for inull in range(nnulls):
    if (signs[i] > 0):
      colour = [255,128,128]
    elif (signs[i] < 0):
      colour = [128,128,255]
    else:
      colour = [128,255,128]
    
    filename = 'output/ringidl{0:04d}.dat'.format(inull+1)+string(i+1,'(I4.4)')+'.dat'
    print(filename)
    
    with open(filename) as file:#openr,fan,file,/f77_unformatted
      nrings = np.fromfile(file, dtype=np.int32, count=1)
      print('number of rings: {}'.format(nrings))
      n = np.fromfile(file, dtype=np.int32, count=1)
    
      while (n >= 1):
        x = np.fromfile(file, dtype=np.float64, count=n)
        y = np.fromfile(file, dtype=np.float64, count=n)
        z = np.fromfile(file, dtype=np.float64, count=n)

        break = np.fromfile(file, dtype=np.int32, count=n)

        breaks = np.concatenate([[-1],np.where(break == 1),[n-1]])
        for ib in range(len(breaks-1)):
          if (breaks[ib] != breaks[ib+1]):
            #draw rings
          #oModel -> add, obj_new("IDLgrPolyline", x[breaks[ib]+1:breaks[ib+1]], y[breaks[ib]+1:breaks[ib+1]], z[breaks[ib]+1:breaks[ib+1]], color=colour)

        n = np.fromfile(file, dtype=np.int32, count=1)
end

pro model_add_separators,oModel,frame

  print,'plot separators'
   
    get_nulls,nnulls,signs,r,spine,nulls
    nnulls=1
    for null=1,nnulls do begin
    seps=read_separators('output/sep'+string(null,'(I4.4)')+'.dat')
    
    print, null, ' output/sep'+string(null,'(I4.4)')+'.dat', size(seps)
    if keyword_set(seps) then begin
    for j=0,(size(seps))[1]-1 do begin
    col = [0,160,60] 
      oModel->add,obj_new("IDLgrPolyline",seps[j,*,0],seps[j,*,1],seps[j,*,2],$
        color=col,thick=4)  
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
