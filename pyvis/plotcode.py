import pyvis.read as rd
import matplotlib.pyplot as plt
import numpy as np

br, lons, lats = rd.synmap('../sphericalharmonics/hmi/synmap_20100601.dat')
# br = (br[::2, ::2] + br[1::2, ::2] + br[::2, 1::2] + br[1::2, 1::2])/4
# lons = (lons[::2] + lons[1::2])/2
# lats = (lats[::2] + lats[1::2])/2

plt.figure(figsize=(16/2.55, 8/2.55))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim([0, 2*np.pi])
plt.ylim([np.pi, 0])
plt.yticks([0, np.pi/4, np.pi/2, np.pi*3/4, np.pi], [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
plt.xticks([0, np.pi/2, np.pi, np.pi*3/2, np.pi*2], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])

levels = np.linspace(-20,20,101)
ax = plt.gca()
cont = ax.contourf(lons, np.arccos(lats), br, levels, cmap=plt.cm.RdBu_r, extend='both', zorder=-1)
ax.set_rasterization_zorder(0)
# plt.contourf(lons, np.arccos(lats), br, levels, cmap=plt.cm.RdBu_r, extend='both', rasterized=True)
cb = plt.colorbar(cont, ax=ax, fraction=0.05, pad=0.025)

diff = levels[-1] - levels[0]
cb.set_ticks(np.array([0.0, 0.25, 0.5, 0.75, 1.0])*diff + levels[0])
cb.set_label('Magnetic Field Strength (G)')
plt.tight_layout()

plt.savefig('thesis/20100601_baddata.pgf', dpi=1200)