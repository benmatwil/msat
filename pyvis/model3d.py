from __future__ import division
import numpy as np
import mayavi.mlab as ml
from . import read as rd
from . import fieldline3d as fl
import sys

def add_sepsurf():
  global nskipglob, nulldata

  print('Adding separatrix surface rings')
  
  rings, breaks = rd.rings(filename, breakinfo=True, nskip=nskipglob)
  
  cols = {-1:(0.5,0.5,1), 0:(0.5,1,0.5), 1:(1,0.5,0.5)}

  for inull in nulllist:
    print('Null {}'.format(inull+1))
    for iring, ring in enumerate(rings[inull]):
      print('Ring {}'.format(iring*nskipglob))
      sys.stdout.write("\033[F")
      brks = np.r_[[-1], np.where(breaks[inull][iring] == 1)[0], [breaks[inull][iring].shape[0]-1]] + 1
      for ib0, ib1 in zip(brks[:-1], brks[1:]):
        if ib0 != ib1:
          ml.plot3d(ring[ib0:ib1, 0], ring[ib0:ib1, 1], ring[ib0:ib1, 2], color=cols[nulldata[inull].sign], line_width=1, tube_radius=None)
    sys.stdout.write("\033[F")

def add_fanlines(nlines, nring):
  global bgrid, xx, yy, zz, ds

  print('Adding separatrix surface field lines')
  
  rings = rd.rings(filename, nskip=nskipglob)

  cols = {-1:(0.5,0.5,1), 0:(0.5,1,0.5), 1:(1,0.5,0.5)}

  for inull in nulllist:
    print('Null {}'.format(inull+1))
    sys.stdout.write("\033[F")
    if nring == None:
      ring = rings[inull][len(rings[inull])//5]
    else:
      ring = rings[inull][nring]
    for ipt in range(0, ring.shape[0], ring.shape[0]//nlines):
      startpt = np.array([ring[ipt, 0], ring[ipt, 1], ring[ipt, 2]])

      h = 0.01*ds
      hmin = 0.001*ds
      hmax = 0.5*ds
      epsilon = 1e-5

      line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon)
      dists = np.sqrt((line[:, 0] - nulldata[inull].pos[0])**2 +
        (line[:, 1] - nulldata[inull].pos[1])**2 + (line[:, 2] - nulldata[inull].pos[2])**2)
      imin = dists.argmin()
      if nulldata[inull].sign < 0:
        line = line[0:imin+1, :]
      else:
        line = line[imin:, :]
      
      ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=cols[nulldata[inull].sign], tube_radius=None)

def add_fieldlines(startpts, col=(0,0,0)):
  global bgrid, xx, yy, zz, ds

  print('Adding separatrix surface field lines')
  
  for ipt in range(startpts.shape[0]):
    startpt = startpts[:,ipt]

    h = 0.01*ds
    hmin = 0.001*ds
    hmax = 0.5*ds
    epsilon = 1e-5

    line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon)
    
    ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=col, tube_radius=None)

def add_spines():
  global nskipglob, nulldata

  print('Adding spines')

  cols = {-1:(0,0,1), 0:(0,1,0), 1:(1,0,0)}

  spines = rd.spines(filename)
  nskip0 = nskipglob//2

  for inull in nulllist:
    for spine in spines[inull]:      
      ml.plot3d(spine[::nskip0, 0], spine[::nskip0, 1], spine[::nskip0, 2], color=cols[nulldata[inull].sign], line_width=4, tube_radius=None)

def add_separators():
  global nulldata, nskipglob

  print('Adding separators')

  seps, conn = rd.separators(filename)
  nskip0 = nskipglob//2

  for inull in range(nulldata.number[-1]):
    for sep in seps[inull]:
      if inull+1 in conn[inull]:
        ml.plot3d(sep[::nskip0, 0], sep[::nskip0, 1], sep[::nskip0, 2], color=(0,1,0), line_width=4, tube_radius=None)

def add_nulls(size):
  global nulldata

  print("Adding nulls")

  cols = {-1:(0,0,1), 0:(0,1,0), 1:(1,0,0)}
  
  r = 5*ds*size
  theta, phi = np.mgrid[0:np.pi:101j, 0:2*np.pi:101j]
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  
  for inull in nulllist:
    ml.mesh(x + nulldata[inull].pos[0], y + nulldata[inull].pos[1], z + nulldata[inull].pos[2], color=cols[nulldata[inull].sign])

def add_box():
  global xx, yy, zz, ds
  
  print("Adding box")
  box = np.zeros((3,2), dtype=np.float64)
  box[:, 0] = np.array([xx.min(), yy.min(), zz.min()])
  box[:, 1] = np.array([xx.max(), yy.max(), zz.max()])
  line = np.array([[0,0,0],[1,0,0],[1,0,1],[1,0,0],[1,1,0],[1,1,1],[1,1,0],
    [0,1,0],[0,1,1],[0,1,0],[0,0,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[0,0,1]], dtype=np.float64)

  line[:,0] = line[:,0]*(box[0,1] - box[0,0]) + box[0,0]
  line[:,1] = line[:,1]*(box[1,1] - box[1,0]) + box[1,0]
  line[:,2] = line[:,2]*(box[2,1] - box[2,0]) + box[2,0]

  # dist = min([n_elements(xx), n_elements(yy), n_elements(zz)])*ds/5
  # oModel -> add, obj_new('idlgrtext', 'x', locations=[box[0,0]+dist,box[1,0],box[2,0]], /onglass)
  # oModel -> add, obj_new('idlgrtext', 'y', locations=[box[0,0],box[1,0]+dist,box[2,0]], /onglass)
  # oModel -> add, obj_new('idlgrtext', 'z', locations=[box[0,0],box[1,0],box[2,0]+dist], /onglass)
  
  ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=(0,0,0), tube_radius=None, line_width=1)

def make(fname, addlist, nulls=None, box=True, fieldlines=None, linecolor=(1,1,1), nskip=20,
  nullrad=1, nfanlines=40, nring=None):

  """Makes a 3D visualisation of the output from Skeleton Codes

    fname: name of the field containing the original magnetic fieldlines
    addlist: list of features to be plotted e.g. ['nulls', 'separators'] will
      only plot the nulls and separators
    
    nulls: if None (default), will plot all nulls, otherwise give a list (starting at 1) of nulls to plot
      e.g. nulls=list(range(45,65)) for all nulls between 45 and 64 inclusive
    nullrad: will scale radius of null spheres by this factor (default 1)
      e.g. nullrad=0.5 will halve the size of the nulls
    box: if True (default), plots a box otherwise set to False
    fieldlines: provide a numpy array of shape (3, n) of start points and it will trace magnetic field lines
    linecolor: color of fieldlines (defaults to black)
    nskip: how many rings to skip in plotting (also skips points in spines and separators for speed)"""
  
  global bgrid, xx, yy, zz, nulldata, ds, filename, nskipglob, nulllist

  ml.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(800,800))

  nskipglob = nskip

  filename = fname

  field = rd.field(filename)
  bgrid = np.zeros((field[0].shape[0], field[0].shape[1], field[0].shape[2], 3), dtype=np.float64)
  for i in range(3):
    bgrid[:, :, :, i] = field[i]

  xx = field[3]
  yy = field[4]
  zz = field[5]

  ds = min([np.diff(xx).min(), np.diff(yy).min(), np.diff(zz).min()])

  nulldata = rd.nulls(filename)

  if nulls == None:
    nulllist = list(range(nulldata.number[-1]))
  else:
    nulllist = []
    for inull in nulls:
      nulllist.append(inull-1)

  if box == True: add_box()
  if 'nulls' in addlist: add_nulls(nullrad)
  if 'fanlines' in addlist: add_fanlines(nfanlines, nring)
  if 'spines' in addlist: add_spines()
  if 'sepsurf' in addlist: add_sepsurf()
  if 'separators' in addlist: add_separators()
  if fieldlines is not None: add_fieldlines(fieldlines, col=linecolor)
