function convert, gridcoord, x, y, z, check=check

  r = gridcoord[0]
  t = gridcoord[1]
  p = gridcoord[2]
  
  r = x[floor(r)] + (r-floor(r))(x[ceil(r)]-x[floor(r)])
  t = y[floor(t)] + (t-floor(t))(y[ceil(t)]-y[floor(t)])
  p = z[floor(p)] + (p-floor(p))(z[ceil(p)]-z[floor(p)])
  if keyword_set(check) then print, r,t,p
  return, [r*cos(p)*sin(t), r*sin(p)*sin(t), r*cos(t)]

end

pro plotstruct, n, converge=converge, fan=fan, ball=ball, rsphere=rsphere
  
  openr,10,'data/newmag.dat'
  nx = 0l
  ny = 0l
  nz = 0l
  readu,10,nx,ny,nz
  bgrid = dblarr(nx,ny,nz,3)
  x = dblarr(nx)
  y = dblarr(ny)
  z = dblarr(nz)
  readu,10,bgrid,x,y,z
  close,10
  
  sizebg = size(bgrid)
  xgc = indgen(sizebg[1])
  ygc = indgen(sizebg[2])
  zgc = indgen(sizebg[3])

  nulls = getnullsgordon(/me)
  npos = nulls[n-1].gridpos-1 ; -1 translation from fortran to idl grid spacing

  xyznpos = convert(npos,x,y,z)
  plt = plot3d([xyznpos[0]],[xyznpos[1]],[xyznpos[2]],'red',lines=' ',sym='x', sym_thick=3, sym_size=1.5, dim=[1200,900])
  plt.scale,1.5,1.5,1.5
  
  if not keyword_set(rsphere) then rsphere = 1d-3 ; rsphere from sf_converge/params.f90

  if keyword_set(converge) then begin
    ny=0l
    spine = dblarr(3)
    maxvec = spine
    minvec = spine

    foreach name, ['fandata.dat', 'spinedata.dat'] do begin
      openr,10,name
      readu,10,ny
      r = dblarr(3,ny)
      readu,10,r
      help,r
      if name eq 'fandata.dat' then begin
        r = transpose(r)
        r = r * rsphere
        for i = 0, 2 do r[*,i] = r[*,i] + npos[i]
        for i = 0, n_elements(r)/3-1 do r[i,*] = convert(r[i,*],x,y,z)
        p = plot3d(r[*,0],r[*,1],r[*,2],lines=' ',sym='x', aspect_ratio=1,/overplot)
        
        readu,10,maxvec,minvec
        max1 = maxvec
        min1 = minvec
        maxvec = maxvec*rsphere
        minvec = minvec*rsphere
        maxvec = maxvec + npos
        minvec = minvec + npos
        maxvec = convert(maxvec,x,y,z)
        minvec = convert(minvec,x,y,z)
        p = plot3d([maxvec[0]],[maxvec[1]],[maxvec[2]],'cyan',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3,/aspect_ratio)
        p = plot3d([minvec[0]],[minvec[1]],[minvec[2]],'g',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3,/aspect_ratio)
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
        p = plot3d(r[*,0],r[*,1],r[*,2],'r',lines=' ',sym='x',/overplot,/aspect_ratio)
        
        readu, 10, spine
        spine = spine*rsphere
        spine = spine + npos
        spine = convert(spine,x,y,z)
        p = plot3d([spine[0]],[spine[1]],[spine[2]],'blue',lines=' ',sym='x',/overplot, sym_thick=6, sym_size=3,/aspect_ratio)
      endif
      close,10
    endforeach
  endif
  
  boxedge = dblarr(2,3)
  for i = 0, 2 do boxedge[*,i] = [npos[i] - rsphere, npos[i] + rsphere]
  
  startpts = []
  if keyword_set(fan) then begin
    rad = 0.6*rsphere
    np = 20
    phi = 2*!dpi*((dindgen(np)+0.5)/np)
    for i = 0, np-1 do begin
      startpts = [[startpts],[max1*cos(phi[i])+min1*sin(phi[i])]]
    endfor
    startpts = startpts*rad
  endif
  
  if keyword_set(ball) then begin
    rad = 5*rsphere/2d1
    nt = 6
    np = 6
    theta = !dpi*((dindgen(nt)+0.5)/nt)
    phi = 2*!dpi*((dindgen(np)+0.5)/np)
    for i = 0, nt-1 do begin
      for j = 0, np-1 do begin
        xp = rad*cos(phi[j])*sin(theta[i])
        yp = rad*sin(phi[j])*sin(theta[i])
        zp = rad*cos(theta[i])
        startpts = [[startpts],[xp,yp,zp]]
      endfor
    endfor
  endif

  for i = 0, 2 do startpts[i,*] = startpts[i,*] + npos[i]
  
  for j = 0, n_elements(startpts)/3-1 do begin
    h = rsphere*5d-2
    line = fieldline3d(startpts[*,j], bgrid, xgc,ygc,zgc, h, 0.1d*h, 10d*h, 0.01d*h, boxedge=boxedge, /oneway)
    line = transpose(line)
    for i = 0, n_elements(line)/3-1 do line[i,*] = convert(line[i,*],x,y,z)
    plt = plot3d(line[*,0],line[*,1],line[*,2],'red',/overplot,/aspect_ratio,/perspective)
    
    line = fieldline3d(startpts[*,j], bgrid, xgc,ygc,zgc, -h, 0.1d*h, 10d*h, 0.01d*h, boxedge=boxedge, /oneway)
    line = transpose(line)
    for i = 0, n_elements(line)/3-1 do line[i,*] = convert(line[i,*],x,y,z)
    plt = plot3d(line[*,0],line[*,1],line[*,2],'orange',/overplot,/aspect_ratio)
  endfor
  
  plt.xtitle = "x"
  plt.ytitle = "y"
  plt.ztitle = "z"
end
