function makefield

  field = dblarr(200,360,180,3)
  
  rlim = [-2d,2d]
  plim = [0d,!dpi]
  tlim = [0d,!dpi]
  
  rdiff = rlim[1]-rlim[0]
  pdiff = plim[1]-plim[0]
  tdiff = tlim[1]-tlim[0]
  
  r0 = [1d,1d]
  p0 = [0d,!dpi]
  
  sz = size(field)
  
  nullk = 1
  globalk = 1
  
  for ir = 0, sz[1]-1 do begin
    for ip = 0, sz[2]-1 do begin
      for it = 0, sz[3]-1 do begin
        r = ir*rdiff/(sz[1]-1)
        p = ip*pdiff/(sz[2]-1)
        t = it*tdiff/(sz[3]-1)
        dist = sqrt(r^2 + r0^2 - 2*r*r0*cos(p-p0))
        field[ir,ip,it,0] = total(nullk*(r-r0*cos(p-p0))/dist^5) + globalk*sin(p)
        field[ir,ip,it,1] = total(nullk*r0*sin(p-p0)/dist*5) + globalk*cos(p)
      endfor
    endfor
  endfor

  return, field

end