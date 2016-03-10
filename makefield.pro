function makefield, rs=rs,ps=ps

  field = dblarr(201,361,2)
  
  rlim = [0,2d]
  plim = [0d,!dpi]
  tlim = [0d,!dpi]
  
  rdiff = rlim[1]-rlim[0]
  pdiff = plim[1]-plim[0]
  tdiff = tlim[1]-tlim[0]
  
  r0 = [1d,1d]
  p0 = [0d,!dpi]
  
  sz = (size(field))[1:3]-1
  
  nullk = 1
  globalk = 1
  
  rs = dindgen(sz[0]+1)*rdiff/sz[0]+rlim[0]
  ps = dindgen(sz[1]+1)*pdiff/sz[1]+plim[0]
  ts = dindgen(sz[2]+1)*tdiff/sz[2]+tlim[0]
  
  ;for it = 0, sz[2] do begin
    for ip = 0, sz[1] do begin
      for ir = 0, sz[0] do begin
        r = rs[ir]
        p = ps[ip]
        ;t = ts[it]
        dist = sqrt(r^2 + r0^2 - 2*r*r0*cos(p-p0))
        field[ir,ip,0] = total(nullk*(r-r0*cos(p-p0))/dist^5) + globalk*sin(p)
        field[ir,ip,1] = total(nullk*r0*sin(p-p0)/dist*5) + globalk*cos(p)
      endfor
    endfor
  ;endfor

  return, field

end