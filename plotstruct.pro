pro plotstruct, n  
  
  openr,10,'data/newmag.dat'
  nx = 0l
  ny = 0l
  nz = 0l
  readu,10,nx,ny,nz
  bgrid = dblarr(nx,ny,nz,3)
  x=dblarr(nx)
  y=dblarr(ny)
  z=dblarr(nz)
  readu,10,bgrid,x,y,z
  close,10

  nulls = getnullsgordon(/me)
  npos = nulls[n-1].gridpos-1 ; -1 translation from fortran to idl grid spacing
print, npos

  ;r = [2.915757d,534.719921d,1093.833099d] ; null 1787
  ;r = [1.00023d,266.802554d,863.767486d] ; null 1443
  ;r = [2.920712d,471.302378d,1.423052d] ; null 1

  ;x = [floor(r[0]),ceil(r[0])]
  ;y = [floor(r[1]),ceil(r[1])]
  ;z = [floor(r[2]),ceil(r[2])]
  ;bgrid = bgrid[x,y,z,*]

  boxedge = dblarr(2,3)
  rsphere = 1d-3 ; rsphere from sf_converge/params.f90
  boxedge[0,0] = npos[0] - rsphere
  boxedge[0,1] = npos[1] - rsphere
  boxedge[0,2] = npos[2] - rsphere
  boxedge[1,0] = npos[0] + rsphere
  boxedge[1,1] = npos[1] + rsphere
  boxedge[1,2] = npos[2] + rsphere

  startpts = []

  if 1 eq 1 then begin ; null 1
    rad = 2d-4
    nt = 5
    np = 5
    theta = !dpi*((dindgen(nt)+0.5)/nt)
    phi = 2*!dpi*((dindgen(np)+0.5)/np)
    print, theta
    print, phi
    for i = 0, nt-1 do begin
      for j = 0, np-1 do begin
        xp = rad*cos(phi[j])*sin(theta[i])
        yp = rad*sin(phi[j])*sin(theta[i])
        zp = rad*cos(theta[i])
        
        startpts = [[startpts],[xp,yp,zp]]
      endfor
    endfor
  endif

  if 0 eq 1 then begin ; null 1787
    np = 20
    rad = 2d-4
    theta = !dpi/2-0.1
    phi = 2*!dpi*(dindgen(np)/(np-1))
    for j = 0, np-1 do begin
      xp = rad*cos(phi[j])*sin(theta)
      yp = rad*sin(phi[j])*sin(theta)
      zp = rad*cos(theta)
      
      startpts = [[startpts],[xp,yp,zp]]
    endfor
  endif

  help, startpts
  print, npos
print, startpts
  for i = 0, 2 do startpts[i,*] = startpts[i,*] + npos[i]
  print, startpts  
  for j = 0, n_elements(startpts)/3-1 do begin

    h = 1d-5
    line = fieldline3d(startpts[*,j], bgrid, x,y,z, h, 0.1d*h, 10d*h, 0.01d*h, boxedge=boxedge, /oneway)

    line = transpose(line)
    for i = 0, n_elements(line)/3-1 do begin
      line[i,*] = convert(line[i,*],x,y,z)
    endfor
    
    overplot = 1
    plt = plot3d(line[*,0],line[*,1],line[*,2],overplot=overplot)
    
    line = fieldline3d(startpts[*,j], bgrid, x,y,z, -h, 0.1d*h, 10d*h, 0.01d*h, boxedge=boxedge, /oneway)
    line = transpose(line)
    for i = 0, n_elements(line)/3-1 do begin
      line[i,*] = convert(line[i,*],x,y,z)
    endfor

    plt = plot3d(line[*,0],line[*,1],line[*,2],'cyan',overplot=overplot)
    
  endfor



  n=0l
  spine = dblarr(3)
  maxvec = spine
  minvec = spine

  foreach name, ['fandata.dat', 'spinedata.dat'] do begin
    openr,10,name
    readu,10,n
    r = dblarr(3,n)
    readu,10,r
    help,r
    if name eq 'fandata.dat' then begin
      r = transpose(r)
      r = r * rsphere
      for i = 0,2 do r[*,i] = r[*,i] + npos[i]
      for i = 0, n_elements(r)/3-1 do r[i,*] = convert(r[i,*],x,y,z)
      p = plot3d(r[*,0],r[*,1],r[*,2],lines=' ',sym='x', aspect_ratio=1,/overplot)
      
      readu,10,maxvec,minvec
      maxvec = maxvec*rsphere
      minvec = minvec*rsphere
      maxvec = maxvec + npos
      minvec = minvec + npos
      maxvec = convert(maxvec,x,y,z)
      minvec = convert(minvec,x,y,z)
      p = plot3d([maxvec[0]],[maxvec[1]],[maxvec[2]],'cyan',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
      p = plot3d([minvec[0]],[minvec[1]],[minvec[2]],'g',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
      
      if 0 eq 1 then begin
        readu, 10, n & crossfan = dblarr(3,n) & readu, 10, crossfan & help, crossfan
        p = plot3d(reform(crossfan[0,*]),reform(crossfan[1,*]),reform(crossfan[2,*]),'orange',lines=' ',sym='x',/overplot, sym_thick=3, sym_size=2)
        cross = dblarr(3,n_elements(r)-1)
        for i = 1, n_elements(r)/3-1 do begin
          cross[*,i-1] = crossp(r[*,0],r[*,i])
          cross[*,i-1] = cross[*,i-1]/sqrt(cross[0,i-1]^2+cross[1,i-1]^2+cross[2,i-1]^2)
        endfor
        print, cross
        p = plot3d(reform(cross[0,*]),reform(cross[1,*]),reform(cross[2,*]),lines=' ',sym='x',/overplot, sym_thick=4, sym_size=3)
      endif
      
    endif
    if name eq 'spinedata.dat' then begin
      r = transpose(r)
      r = r * rsphere
      for i = 0,2 do r[*,i] = r[*,i] + npos[i]
      for i = 0, n_elements(r)/3-1 do r[i,*] = convert(r[i,*],x,y,z)
      p = plot3d(r[*,0],r[*,1],r[*,2],'r',lines=' ',sym='x',/overplot)
      
      readu, 10, spine
      spine = spine*rsphere
      spine = spine + npos
      spine = convert(spine,x,y,z)
      p = plot3d([spine[0]],[spine[1]],[spine[2]],'blue',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3)
    endif
    close,10
  endforeach

end