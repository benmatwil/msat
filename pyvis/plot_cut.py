import numpy as np
import matplotlib.pyplot as plt
import pyvis.read_data as rd
import os

plt.ion()

def start(r, filename):
  global datafile, prefile, nulldata

  datafile = filename
  prefile = rd.prefix(filename)

  nulldata = rd.nulls(datafile, simple=True)

  os.system(f'./make_cut -i {filename} -r {r}')

  plt.figure(figsize=(10,5))
  plt.plot([])
  plt.xlabel('Longitude')
  plt.xlim([0,2*np.pi])
  plt.xticks([0, np.pi/2, np.pi, np.pi*3/2, np.pi*2], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
  plt.ylabel('Latitude')
  plt.ylim([np.pi,0])
  plt.yticks([0, np.pi/4, np.pi/2, np.pi*3/4, np.pi], [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
  plt.tight_layout()

def spines():
  global prefile

  with open('output/'+prefile+'-cut_spines.dat', 'rb') as spinefile:
    null = np.fromfile(spinefile, dtype=np.int32, count=1)
    while null > 0:
      spine = np.fromfile(spinefile, dtype=np.float64, count=3)
      plt.plot(spine[2], spine[1], '.', c='orange')
      null = np.fromfile(spinefile, dtype=np.int32, count=1)

def separators():
  global prefile

  with open('output/'+prefile+'-cut_seps.dat', 'rb') as sepfile:
    null = np.fromfile(sepfile, dtype=np.int32, count=1)
    while null > 0:
      sep = np.fromfile(sepfile, dtype=np.float64, count=3)
      plt.plot(sep[2], sep[1], '.', c='red')
      null = np.fromfile(sepfile, dtype=np.int32, count=1)

def ringssort():
  global prefile, nulldata

  with open('output/'+prefile+'-cut_rings.dat', 'rb') as ringfile:
    for i in range(nulldata.number.max()):
      length = np.asscalar(np.fromfile(ringfile, dtype=np.int32, count=1))
      ring = np.fromfile(ringfile, dtype=np.float64, count=3*length).reshape(-1,3)
      if length > 0:
        ring2 = ring.copy()
        ring2[0,:] = ring2[-1,:]
        ring2[1:,0] = 0
        ring = ring[:-1,:]
        for iring in range(ring2.shape[0]-1):
          ptlist = set(range(ring.shape[0]))
          dists = np.zeros(ring.shape[0], dtype=np.float64)
          for ii, jring in enumerate(ptlist):
            dists[ii] = np.sum((ring2[iring, :] - ring[jring, :])**2)
          ring2[iring+1] = ring[dists.argmin(), :]
          ptlist.remove(dists.argmin())
          ring = ring[list(ptlist),:]
        plt.plot(ring2[:,2], ring2[:,1])

def rings():
  global prefile, nulldata

  with open('output/'+prefile+'-cut_rings.dat', 'rb') as ringfile:
    for i in range(nulldata.number.max()):
      length = np.asscalar(np.fromfile(ringfile, dtype=np.int32, count=1))
      ring = np.fromfile(ringfile, dtype=np.float64, count=3*length).reshape(-1,3)
      plt.plot(ring[:,2], ring[:,1])

def nulls():
  global datafile, nulldata

  plt.plot(nulldata.pos[:,2], nulldata.pos[:,1], '.', c='green')

def field():
  global datafile

  br, _, _, rads, thetas, phis = rd.field(datafile)
  levels = np.linspace(-10,10,101)
  plt.contourf(phis, thetas, br[0,:,:], levels, cmap=plt.cm.RdBu_r, extend='both', vmax=10, vmin=-10)
  cb = plt.colorbar(fraction=0.05, pad=0.025)
  cb.set_ticks([-10,-5,0,5,10])
  cb.set_label('Magnetic Field Strength (G)')
  plt.tight_layout()
