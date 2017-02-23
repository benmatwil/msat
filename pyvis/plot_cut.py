import numpy as np
import matplotlib.pyplot as plt
import pyvis.read_data as rd
import os

plt.ion()

def start(r, filename):
  global datafile, prefile

  datafile = filename
  prefile = rd.prefix(filename)

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

  spines = np.fromfile('output/'+prefile+'-cut_spines.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(spines[:,2], spines[:,1], '.', c='orange')

def separators():
  global prefile

  seps = np.fromfile('output/'+prefile+'-cut_seps.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(seps[:,2], seps[:,1], '.', c='red')

def rings():
  global prefile

  rings = np.fromfile('output/'+prefile+'-cut_rings.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(rings[:,2], rings[:,1], '.', c='blue')

def nulls():
  global datafile

  nulls = rd.nulls(datafile, simple=True)
  plt.plot(nulls.pos[:,2], nulls.pos[:,1], '.', c='green')

def field():
  global datafile

  br, _, _, rads, thetas, phis = rd.field(datafile)
  levels = np.linspace(-10,10,101)
  plt.contourf(phis, thetas, br[0,:,:], levels, cmap=plt.cm.RdBu_r, extend='both', vmax=10, vmin=-10)
  cb = plt.colorbar(fraction=0.05, pad=0.025)
  cb.set_ticks([-10,-5,0,5,10])
  cb.set_label('Magnetic Field Strength (G)')
  plt.tight_layout()
