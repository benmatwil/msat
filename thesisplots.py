import pyvis.read as rd
import pyvis.fieldline3d as fl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

xsize = 12
xthick = 3

# fig = plt.figure(figsize=(10/2.55,20/2.55), tight_layout=True)
# axes = []
# for i in range(1,9):
#     axes.append(fig.add_subplot(4,2,i, projection='3d'))
#     axes[i-1].axis('off')

for ifield in range(1,5):
    filename = 'data/nullfield{}.dat'.format(ifield)

    nulls = rd.nulls(filename)
    field = rd.field(filename)
    bgrid = np.zeros((field[0].shape[0], field[0].shape[1], field[0].shape[2], 3), dtype=np.float64)
    for i in range(3):
        bgrid[:, :, :, i] = field[i]

    xs = field[3].copy()
    ys = field[4].copy()
    zs = field[5].copy()

    del field

    n=1
    os.system('./sf -i {0} -n {1:04d}'.format(filename, n) )

    fnames = ['maxvec.dat', 'minvec.dat', 'spine.dat']
    col = ['blue', 'green', 'red']

    for iplot in range(2):
        fig = plt.figure(figsize=(7.5/2.55,7.5/2.55), tight_layout=True)
        ax = fig.add_subplot(111, projection='3d')
        ax.axis('off')
        if iplot == 0:
            rsphere = 0.1
            for j, name in enumerate(fnames):
                with open('output/' + name, 'rb') as file:
                    ny = np.asscalar(np.fromfile(file, count=1, dtype=np.int32))
                    r = np.fromfile(file, count=3*ny, dtype=np.float64).reshape(ny,3).T
                    r = r * rsphere

                    for i in range(3):
                        r[i,:] = r[i,:] + nulls[0].pos[i]

                    plt.plot(r[0,:],r[1,:],r[2,:], '. ', c=col[j])

                    if name == fnames[0]:
                        maxvec = np.fromfile(file, count=3, dtype=np.float64)
                        max1 = maxvec.copy()
                        maxvec = maxvec*rsphere
                        maxvec = maxvec + nulls[0].pos
                        plt.plot([maxvec[0]],[maxvec[1]],[maxvec[2]],'+ ', c=col[j], ms=xsize, mew=xthick)

                    if name == fnames[1]:
                        minvec = np.fromfile(file, count=3, dtype=np.float64)
                        min1 = minvec.copy()
                        minvec = minvec*rsphere
                        minvec = minvec + nulls[0].pos
                        plt.plot([minvec[0]],[minvec[1]],[minvec[2]],'+ ', c=col[j], ms=xsize, mew=xthick)

                    if name == fnames[2]:
                        spine = np.fromfile(file, count=3, dtype=np.float64)
                        spine = spine*rsphere
                        spine = spine + nulls[0].pos
                        plt.plot([spine[0]],[spine[1]],[spine[2]],'+ ', c=col[j], ms=xsize, mew=xthick)

            fanpts = []
            npts = 121
            for iangle in range(npts):
                angle = 2*np.pi*iangle/(npts-1)
                fanpts = fanpts + [(max1*np.cos(angle) + min1*np.sin(angle))*rsphere + nulls[0].pos]
            fanpts = np.array(fanpts)
            ax.plot(fanpts[:,0], fanpts[:,1], fanpts[:,2], c='black', lw=0.8)

            spineline = np.array([spine, 2*nulls[0].pos - spine])
            ax.plot(spineline[:,0], spineline[:,1], spineline[:,2], c='black', lw=0.8)
        else:
            ax = fig.add_subplot(111, projection='3d')
            ax.axis('off')

            startpts = []
            rad = 0.05
            ntheta = 8
            nphi = 8
            theta = np.pi*((np.arange(ntheta)+0.5)/ntheta)
            phi = 2*np.pi*((np.arange(nphi)+0.5)/nphi)
            for it in theta:
                for ip in phi:
                    xp = rad*np.cos(ip)*np.sin(it)
                    yp = rad*np.sin(ip)*np.sin(it)
                    zp = rad*np.cos(it)
                    startpts.append(np.array([xp, yp, zp]))

            h0 = np.diff(xs).min()*0.1
            for startpt in startpts:
                h = h0
                line = fl.fieldline3d(startpt + nulls[0].pos, bgrid, xs, ys, zs, h, 0.1*h, 10*h, 0.01*h)
                ax.plot(line[:, 0], line[:, 1], line[:, 2], c='red', lw=0.2)

            ax.set_xlim((xs.min(),xs.max()))
            ax.set_ylim((ys.min(),ys.max()))
            ax.set_zlim((zs.min(),zs.max()))

            ax.plot(fanpts[:,0], fanpts[:,1], fanpts[:,2], c='black', lw=0.8)
            ax.plot(spineline[:,0], spineline[:,1], spineline[:,2], c='black', lw=0.8)

        ax.view_init(elev=20)
        ax.dist = 5.8

        plt.savefig('sf{}-{}.pgf'.format(ifield, iplot), dpi=1200)
        plt.close()
