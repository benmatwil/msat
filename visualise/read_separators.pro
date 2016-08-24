function read_separators,fname

  get_lun, sep
  openr, sep, fname

  opt = 0l
  ringsmax = 0l

  readu, sep, opt

  seps = !null

  while opt gt 0 do begin
    readu, sep, ringsmax
    x = dblarr(opt)
    y = dblarr(opt)
    z = dblarr(opt)
    readu, sep, x, y, z
    
    slice = dblarr(1,ringsmax,3)
    slice(0,0:opt-1,0) = x
    slice(0,0:opt-1,1) = y
    slice(0,0:opt-1,2) = z
    slice(0,opt:ringsmax-1,0) = x(opt-1)
    slice(0,opt:ringsmax-1,1) = y(opt-1)
    slice(0,opt:ringsmax-1,2) = z(opt-1)
    
    seps = [seps,slice]
    
    readu, sep, opt
  endwhile

  return,seps

end
