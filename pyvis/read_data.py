import numpy as np

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
  
  with open('output/'+filename[5:-4]+'-nullpos.dat', 'rb') as nullfile:
    nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
    gridpos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
    pos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
  
  if simple == False:
    with open('output/'+filename[5:-4]+'-nulldata.dat', 'rb') as nullfile:
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

def separators(filename):

  nulldata = nulls(filename, simple=True)

  seps = []
  for inull in nulldata.number:
    sepsi = []
    with open('output/'+filename[5:-4]+'-sep{:04d}.dat'.format(inull), 'rb') as sepfile:
      length = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))
      while length > 0:
        pts = np.fromfile(sepfile, dtype=np.float64, count=3*length).reshape(length,3)
        # pts = np.empty((length, 3), dtype=np.float64)
        # pts[:,0] = np.fromfile(sepfile, dtype=np.float64, count=length)
        # pts[:,1] = np.fromfile(sepfile, dtype=np.float64, count=length)
        # pts[:,2] = np.fromfile(sepfile, dtype=np.float64, count=length)
        sepsi = sepsi + [pts]
        length = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))
    seps = seps + [sepsi]
  
  return seps

def connectivity(filename):
  
  nulldata = nulls(filename, simple=True)
  
  connectdata = []

  for inull in nulldata.number:
    info = []
    with open('output/'+filename[5:-4]+'-separator{:04d}.dat'.format(inull), 'rb') as connectfile:
      flag = np.asscalar(np.fromfile(connectfile, dtype=np.int32, count=1))
      while flag > 0:
        info = info + [np.asscalar(np.fromfile(connectfile, dtype=np.int32, count=1))]
        np.fromfile(connectfile, dtype=np.int32, count=2)
        flag = np.asscalar(np.fromfile(connectfile, dtype=np.int32, count=1))
    
    connectdata = connectdata + [info]

  return connectdata