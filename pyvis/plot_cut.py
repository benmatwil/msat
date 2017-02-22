import numpy as np
import matplotlib.pyplot as plt
import pyvis.read_data as rd
import os

plt.ion()

def start(r, filename):

  global datafile

  datafile = filename

  os.system(f'./make_cut -i {filename} -r {r}')

  plt.figure(figsize=(10,5))
  plt.plot([])
  plt.xlim([0,2*np.pi])
  plt.xticks([0, np.pi/2, np.pi, np.pi*3/2, np.pi*2], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
  plt.ylim([np.pi,0])
  plt.yticks([0, np.pi/4, np.pi/2, np.pi*3/4, np.pi], [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])

def spines():
  spines = np.fromfile('output/field_20100601_0081-fft-cut_spines.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(spines[:,2], spines[:,1], '.', c='orange')

def separators():
  seps = np.fromfile('output/field_20100601_0081-fft-cut_seps.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(seps[:,2], seps[:,1], '.', c='red')

def rings():
  rings = np.fromfile('output/field_20100601_0081-fft-cut_rings.dat', dtype=np.float64).reshape(-1,3)
  plt.plot(rings[:,2], rings[:,1], '.', c='blue')

def nulls():
  global datafile
  nulls = rd.nulls(file, simple=True)
  plt.plot(nulls.pos[:,2], nulls.pos[:,1], '.', c='green')