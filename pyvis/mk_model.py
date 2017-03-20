from __future__ import division
import numpy as np
import mayavi.mlab as ml
from . import read as rd
from . import fieldline3d as fl


def model_add_sepsurf():
  global nskipglob, nulldata

  print('Adding separatrix surface rings')
  
  rings, breaks = rd.rings(filename, breakinfo=True, nskip=nskipglob)
  
  cols = {-1:(0.5,0.5,1), 0:(0.5,1,0.5), 1:(1,0.5,0.5)}

  for inull, null in enumerate(nulldata):
    for iring, ring in enumerate(rings[inull]):
      brks = np.r_[[-1], np.where(breaks[inull][iring] == 1)[0], [breaks[inull][iring].shape[0]-1]] + 1
      for ib0, ib1 in zip(brks[:-1], brks[1:]):
        if ib0 != ib1:
          ml.plot3d(ring[ib0:ib1, 0], ring[ib0:ib1, 1], ring[ib0:ib1, 2], color=cols[null.sign], line_width=1, tube_radius=None)

def model_add_fanlines():
  global bgrid, xx, yy, zz, ds

  print('Adding separatrix surface field lines')
  
  mxline = 2000

  rings = rd.rings(filename, nskip=nskipglob)

  cols = {-1:(0.5,0.5,1), 0:(0.5,1,0.5), 1:(1,0.5,0.5)}

  for inull, null in enumerate(nulldata):

    ring = rings[inull][10]
    for ipt in range(0, ring.shape[0], ring.shape[0]//40):
      startpt = np.array([ring[ipt, 0], ring[ipt, 1], ring[ipt, 2]])

      h = 0.01*ds
      hmin = 0.001*ds
      hmax = 0.5*ds
      epsilon = 1e-5

      line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon)
      # dists = (line[:, 0] - null.pos[0])**2 + (line[:, 1] - null.pos[0])**2 + (line[:, 2] - null.pos[0])**2
      # imin = dists.argmin()
      # if null.sign > 0:
      #   line = line[0:imin+1, :]
      # else:
      #   line = line[imin:, :]
      
      ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=cols[null.sign], tube_radius=None)


def model_add_spines():
  global nskipglob, nulldata

  print('Adding spines')

  cols = {-1:(0,0,1), 0:(0,1,0), 1:(1,0,0)}

  spines = rd.spines(filename)
  nskip0 = nskipglob//2

  for inull, null in enumerate(nulldata):
    for spine in spines[inull]:      
      ml.plot3d(spine[::nskip0, 0], spine[::nskip0, 1], spine[::nskip0, 2], color=cols[null.sign], line_width=4, tube_radius=None)

def model_add_separators():
  global nulldata, nskipglob

  print('Adding separators')

  seps = rd.separators(filename, connectivity=False)
  nskip0 = nskipglob//2

  for inull in range(nulldata.shape[0]):
    for sep in seps[inull]:
      ml.plot3d(sep[::nskip0, 0], sep[::nskip0, 1], sep[::nskip0, 2], color=(0,1,0), line_width=4, tube_radius=None)

def model_add_nulls():
  global nulldata

  print("Adding nulls")

  radius = 15*ds

  cols = {-1:(0,0,1), 0:(0,1,0), 1:(1,0,0)}
  
  for sign in [-1, 1]:
    xpos = nulldata[nulldata.sign == sign].pos[:,0]
    ypos = nulldata[nulldata.sign == sign].pos[:,1]
    zpos = nulldata[nulldata.sign == sign].pos[:,2]
    
    ml.points3d(xpos, ypos, zpos, radius, color=cols[sign], scale_factor=1)

def model_add_box():
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
  
  ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=(0,0,0), tube_radius=None)

def mk_model_mag_sep(fname, nulls=False, separators=False, sepsurf=False, spines=False, box=False, fanlines=False, nskip=20):
  
  global bgrid, xx, yy, zz, nulldata, ds, filename, nskipglob

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

  if box == True: model_add_box()
  if nulls == True: model_add_nulls()
  if fanlines == True: model_add_fanlines()
  if spines == True: model_add_spines()
  if sepsurf == True: model_add_sepsurf()
  if separators == True: model_add_separators()
