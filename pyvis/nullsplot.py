import read_data as rd
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

filename = 'data/field_20100601_0351.dat'
nulls = rd.nulls(filename, simple=True)

nulls = nulls[np.argsort(nulls.pos[:,0])]
nnulls = nulls.shape[0]
plt.plot(nulls.pos[:,0]-1, np.arange(nnulls)+1)
ax = plt.gca()
ax.set_xscale('log')
plt.xticks([1e-6, 1e-4, 1e-2, 1])
ax.set_xticklabels([r'$1.000001R_{\odot}$', r'$1.0001R_{\odot}$', r'$1.01R_{\odot}$', r'$2R_{\odot}$'])
ax.set_xlabel('Distance above solar surface $R_{\odot}$')
ax.set_ylabel('Null number (ordered in height order)')
plt.title('Null Heights ($l_{max} = 351$)')

plt.savefig('nullheights.png', dpi=400)
plt.close()

###################################################

plt.figure()
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('$r$')
ax.set_ylabel('Number of nulls in the region $>r$')
plt.title('Cumulative number of nulls')

ax.set_xticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
ax.set_xticklabels(['$1.0000001R_{\odot}$', '$1.000001R_{\odot}$',
  '$1.00001R_{\odot}$', '$1.0001R_{\odot}$',
  '$1.001R_{\odot}$', '$1.01R_{\odot}$',
  '$1.1R_{\odot}$', '$2R_{\odot}$'])

date = '20100601'
names = ['', 'cubicr']
linesty = ['-', '--']
lmxs = [81, 251, 351, 701, 1001]
cols = ['blue', 'red', 'green', 'purple', 'orange']

for name, sty in zip(names, linesty):
  for lmax, col in zip(lmxs, cols):
    filename = 'data/field_' + date + '_{:04d}{}.dat'.format(lmax,name)
    nulls = rd.nulls(filename, simple=True)
    nulls = nulls[np.argsort(nulls.pos[:,0])]
    nulls = nulls[::-1]
    plt.plot(nulls.pos[:,0]-1, np.arange(nulls.shape[0]), ls=sty,
      c=col, label=r'$l_{\max} = ' + r'{0}$'.format(lmax) + ' '+name)

# filename = 'data/field_20100601_0351_extrar.dat'
# nulls = rd.nulls(filename, simple=True)
# rads = np.linspace(1,2.5,nrads)
# nulls = nulls[np.argsort(nulls.pos[:,0])]
# cumr = np.empty_like(rads, dtype=np.int32)
# for ir in xrange(rads.shape[0]):
#   cumr[ir] = np.sum(nulls.pos[:,0] > rads[ir])

# plt.plot(rads-1, cumr, label=r'$l_{\max} = ' + r'{0}$'.format(351) + 'newr')

br, bt, bp, rads, thetas, phis = rd.field(filename)
rads = rads-1
for rad in rads:
  plt.axvline(rad, ls='--', c='green')
  #plt.plot([rad,rad],[0,nulls.shape[0]],'--')

cols = ['blue', 'green', 'red', 'teal', 'purple']
for lmax, col in zip(lmxs, cols):
  filename = 'data/field_' + date + '_{:04d}.dat'.format(lmax)
  _,_,_,rads,_,_ = rd.field(filename)
  plt.axvline(rads[1]-1, ls='--', c=col)

plt.savefig('cumulativenullsr.png', dpi=400)
plt.close()

######################################################

plt.figure(figsize=(16,9))
br, bt, bp, rads, thetas, phis = rd.field(filename)
nulls = rd.nulls(filename, simple=True)
plt.contourf(phis, thetas, br[0,:,:], np.linspace(-10,10,21), cmap='RdBu_r', extend='both')

plt.plot(nulls.pos[:,2], nulls.pos[:,1], 'go ')
plt.savefig('nulldist.png', dpi=400)
plt.close()

phistrips = np.linspace(0, 2*np.pi, 121)
nstrip = np.zeros(phistrips.shape[0]-1, dtype=np.int32)
for istrip in xrange(phistrips.shape[0]-1):
  n = 0
  for inull in xrange(nulls.shape[0]):
    if (nulls.pos[inull,2] > phistrips[istrip] and nulls.pos[inull,2] < phistrips[istrip+1]):
      n += 1
  nstrip[istrip] = n

plt.plot(180*(phistrips[:-1]+phistrips[1:])/np.pi/2, nstrip*6)

plt.savefig('nulldistnum.png', dpi=400)
plt.close()