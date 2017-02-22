import numpy as np

def prefix(filename):
  return filename[5:-4]

def synmap(filename):

  with open(filename, 'rb') as fieldfile:
    nlon = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
    nlat = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
    bsyn = np.fromfile(fieldfile, count=nlon*nlat, dtype=np.float64).reshape(nlat,nlon)
    lons = np.fromfile(fieldfile, count=nlon, dtype=np.float64)
    lats = np.fromfile(fieldfile, count=nlat, dtype=np.float64)

    return bsyn, lons, lats

def field(filename):
  
  with open(filename, 'rb') as fieldfile:
    nrad = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
    ntheta = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
    nphi = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
    br = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
    bt = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
    bp = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
    rads = np.fromfile(fieldfile, count=nrad, dtype=np.float64)
    thetas = np.fromfile(fieldfile, count=ntheta, dtype=np.float64)
    phis = np.fromfile(fieldfile, count=nphi, dtype=np.float64)
    
  return br, bt, bp, rads, thetas, phis

def nulls(filename, simple=False):
  
  with open('output/'+prefix(filename)+'-nullpos.dat', 'rb') as nullfile:
    nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
    gridpos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
    pos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
  
  if simple == False:
    with open('output/'+prefix(filename)+'-nulldata.dat', 'rb') as nullfile:
      nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
      signs = np.fromfile(nullfile, dtype=np.int32, count=nnulls)
      spines = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
      fans = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
      warning = np.fromfile(nullfile, dtype=np.float64, count=nnulls)

    nulls = np.recarray(nnulls, dtype=[('number',np.int32),
      ('pos',np.float64,3),('gridpos',np.float64,3),
      ('sign',np.int32),('spine',np.float64,3),('fan',np.float64,3)])
  
    nulls.sign = signs
    nulls.spine = spines
    nulls.fan = fans
    nulls.warning = warning
  else:
    nulls = np.recarray(nnulls, dtype=[('number',np.int32),
      ('pos',np.float64,3),('gridpos',np.float64,3)])
  
  nulls.number = np.arange(nnulls, dtype=np.int32)+1
  nulls.pos = pos
  nulls.gridpos = gridpos

  return nulls

def separators(filename, lines=True, connectivity=True):

  nulldata = nulls(filename, simple=True)

  conlist = []
  seplist = []
  coni = []
  sepi = []
  inull = 1

  with open('output/'+prefix(filename)+'-connectivity.dat', 'rb') as sepinfo:
    with open('output/'+prefix(filename)+'-separators.dat', 'rb') as seps:
      start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
      while start > 0 or inull <= 86:
        if (inull == start):
          end = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
          length = np.asscalar(np.fromfile(seps, dtype=np.int32, count=1))
          separator = np.fromfile(seps, dtype=np.float64, count=3*length).reshape(-1,3)
          coni += [end]
          sepi += [separator]
          start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
        else:
          inull += 1
          conlist += [coni]
          seplist += [sepi]
          sepi = []
          coni = []
  
  if lines == False:
    return conlist
  if connectivity == False:
    return seplist
  else:
    return seplist, conlist

def spines(filename):
  
  nulldata = nulls(filename, simple=True)

  spinelist = []

  with open('output/'+prefix(filename)+'-spines.dat', 'rb') as spinefile:
    for inull in nulldata.number:
      spinelisti = []
      for ispine in range(2):
        length = np.asscalar(np.fromfile(spinefile, dtype=np.int32, count=1))
        spinelisti += [np.fromfile(spinefile, dtype=np.float64, count=3*length).reshape(-1,3)]
      spinelist += [spinelisti]
  
  return spinelist

def rings(filename, breakinfo=False):

  nulldata = nulls(filename, simple=True)

  ringlist = []

  with open('output/'+prefix(filename)+'-ringinfo.dat', 'rb') as ringinfo:
    with open('output/'+prefix(filename)+'-rings.dat', 'rb') as ringfile:
      ringsmax = np.asscalar(np.fromfile(ringinfo, dtype=np.int32, count=1))
      for inull in nulldata.number:
        breaklisti = []
        ringlisti = []
        lengths = np.fromfile(ringinfo, dtype=np.int32, count=ringsmax)
        iring = 0
        while iring < ringsmax and lengths[iring] > 0:
          breaklisti += [np.fromfile(ringfile, dtype=np.int32, count=lengths[iring])]
          ringlisti += [np.fromfile(ringfile, dtype=np.float64, count=3*lengths[iring]).reshape(-1,3)]
          iring += 1
        ringlist += [ringlisti]
        breaklist += [breaklisti]

  if breakinfo == True:
    return
  else:
    return ringlist
