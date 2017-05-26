from . import read as rd
import numpy as np
import pyvis.fieldline3d as fl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def plot(n, filename, converge=True, fan=False, ball=True, rsphere=1e-4, h0=3e-2, badtest=False):

    if converge == True:
        os.system('make sf mode=debug')
        os.system('./sf -i {0} -n {1:04d}'.format(filename, n) )

    field = rd.field(filename)
    bgrid = np.zeros((field[0].shape[0], field[0].shape[1], field[0].shape[2], 3), dtype=np.float64)
    for i in range(3):
        bgrid[:, :, :, i] = field[i]
    
    rads = field[3]
    thetas = field[4]
    phis = field[5]

    del field
    
    xgc = np.arange(bgrid.shape[0])
    ygc = np.arange(bgrid.shape[1])
    zgc = np.arange(bgrid.shape[2])

    nulls = rd.nulls(filename, simple=True)
    nullpos = nulls[n-1].gridpos-1 # -1 translation from fortran to idl grid spacing

    xsize = 20
    xthick = 5

    plt.figure(figsize=(12,12), tight_layout=True)
    plt.axes(projection='3d')
    plt.plot([nullpos[0]],[nullpos[1]],[nullpos[2]], 'rx ', ms=xsize, mew=xthick)

    plt.axis('off')

    if converge == True:
        fnames = ['maxvec.dat', 'minvec.dat', 'spine.dat']
        col = ['black', 'green', 'red']

        for j, name in enumerate(fnames):
            with open('output/' + name, 'rb') as file:
                ny = np.asscalar(np.fromfile(file, count=1, dtype=np.int32))
                r = np.fromfile(file, count=3*ny, dtype=np.float64).reshape(ny,3).T
                r = r * rsphere

                for i in range(3):
                    r[i,:] = r[i,:] + nullpos[i]

                plt.plot(r[0,:],r[1,:],r[2,:], 'x ', c=col[j])

                if name == fnames[0]:
                    maxvec = np.fromfile(file, count=3, dtype=np.float64)
                    max1 = maxvec.copy()
                    maxvec = maxvec*rsphere
                    maxvec = maxvec + nullpos
                    plt.plot([maxvec[0]],[maxvec[1]],[maxvec[2]],'x ', c='cyan', ms=xsize, mew=xthick)

                if name == fnames[1]:
                    minvec = np.fromfile(file, count=3, dtype=np.float64)
                    min1 = minvec.copy()
                    minvec = minvec*rsphere
                    minvec = minvec + nullpos
                    plt.plot([minvec[0]],[minvec[1]],[minvec[2]],'x ', c='green', ms=xsize, mew=xthick)

                if name == fnames[2]:
                    spine = np.fromfile(file, count=3, dtype=np.float64)
                    spine = spine*rsphere
                    spine = spine + nullpos
                    plt.plot([spine[0]],[spine[1]],[spine[2]],'x ', c='blue', ms=xsize, mew=xthick)

        fanpts = []
        npts = 121
        for iangle in range(npts):
            angle = 2*np.pi*iangle/(npts-1)
            fanpts = fanpts + [(max1*np.cos(angle) + min1*np.sin(angle))*rsphere + nullpos]
        fanpts = np.array(fanpts)
        plt.plot(fanpts[:,0], fanpts[:,1], fanpts[:,2], c='black')

        spineline = np.array([spine, 2*nullpos - spine])
        plt.plot(spineline[:,0], spineline[:,1], spineline[:,2], c='blue')

    boxedge = np.zeros((2,3), dtype=np.float64)
    for i in range(3):
        boxedge[:,i] = np.array([nullpos[i] - rsphere, nullpos[i] + rsphere])
    
    startpts = []
    if fan == True:
        rad = 0.6*rsphere
        nphi = 20
        phi = 2*np.pi*(np.arange(nphi)+0.5)/nphi
        for i in range(nphi):
            startpts = startpts + [rad*(max1*np.cos(phi[i])+min1*np.sin(phi[i]))]
    
    if ball == True:
        rad = 5*rsphere/2e1
        ntheta = 6
        nphi = 6
        theta = np.pi*((np.arange(ntheta)+0.5)/ntheta)
        phi = 2*np.pi*((np.arange(nphi)+0.5)/nphi)
        for it in theta:
            for ip in phi:
                xp = rad*np.cos(ip)*np.sin(it)
                yp = rad*np.sin(ip)*np.sin(it)
                zp = rad*np.cos(it)
                startpts = startpts + [np.array([xp,yp,zp])]

    if ball == True or fan == True:
        h0 = rsphere*h0
        countf = 0
        countb = 0
        for startpt in startpts:
            h = h0
            line = fl.fieldline3d(startpt + nullpos, bgrid, xgc,ygc,zgc, h, 0.1*h, 10*h, 0.01*h, boxedge=boxedge, oneway=True, gridcoord=True)
            plt.plot(line[:,0],line[:,1],line[:,2],c='red')
            if (line[-1,0] < boxedge[0,0] + rsphere*1e-3 or line[-1,0] > boxedge[1,0] - rsphere*1e-3 or
                line[-1,1] < boxedge[0,1] + rsphere*1e-3 or line[-1,1] > boxedge[1,1] - rsphere*1e-3 or
                line[-1,2] < boxedge[0,2] + rsphere*1e-3 or line[-1,2] > boxedge[1,2] - rsphere*1e-3):
                countf += 1

            h = h0
            line = fl.fieldline3d(startpt + nullpos, bgrid, xgc,ygc,zgc, -h, 0.1*h, 10*h, 0.01*h, boxedge=boxedge, oneway=True, gridcoord=True)
            plt.plot(line[:,0],line[:,1],line[:,2],c='orange')
            if (line[0,0] < boxedge[0,0] + rsphere*1e-3 or line[0,0] > boxedge[1,0] - rsphere*1e-3 or
                line[0,1] < boxedge[0,1] + rsphere*1e-3 or line[0,1] > boxedge[1,1] - rsphere*1e-3 or
                line[0,2] < boxedge[0,2] + rsphere*1e-3 or line[0,2] > boxedge[1,2] - rsphere*1e-3):
                countb += 1
        
        if badtest == True:
            nlines = len(startpts)
            print(n, countf, countb, nlines)
            plt.close()
            if (countf > 0.9*nlines and countb < 0.1*nlines) or (countf < 0.1*nlines and countb > 0.9*nlines):
                print("Actually probs a bad null ({}) at {}".format(n, nullpos))
                return n
            else:
                plt.close()
        else:
            plt.suptitle(r'File: {} Null: {}'.format(filename, n).replace('_', '\_'))

            print("----------------------------------------")
            print("Key:")
            print("----")
            print("Null - Red cross")
            print("Spine - Dark blue cross")
            print("Spine line - Blue line")
            print("Maxvec - Cyan cross")
            print("Minvec - Green cross")
            print("Fan plane - Black circle")
            print("Converged fan - Small black crosses")
            print("Converged spine - Small red crosses")
            print("Lines traced forwards - Red lines")
            print("Lines traced backwards - Orange Lines")
            print("----------------------------------------")
