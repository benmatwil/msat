pro makefield

  field = dblarr(200,360,180,3)
  
  rlim = [-2d,2d]
  plim = [0d,!dpi]
  tlim = [0d,!dpi]
  
  r0 = [1d,1d]
  p0 = [0d,!dpi]
  
  k = 1
  field[*,*,*,*] = 0d
  foreach r, r0 do begin
    for ir = 0, 199 do begin
      for ip = 0, 359 do begin
        for it = 0, 179 do begin
          r = ir*(rlim[1]-rlim[0])/199
          p = ip*(plim[1]-plim[0])/359
          t = it*(tlim[1]-tlim[0])/179
          dist = r^2 + r0^2 - 2*r*r0*cos(p-p0)
          field[ir,ip,it,0] = field[ir,ip,it,0] + 
          field[ir,ip,it,1] = field[ir,ip,it,1] + 
          field[ir,ip,it,2] = field[ir,ip,it,2] + 
        endfor
      endfor
    endfor
  endforeach

end