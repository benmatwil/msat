import numpy as np
from math import *

def trilinear3d(pt, grid, xx, yy, zz):

    ix = np.where(pt[0] > xx)[0].max()
    iy = np.where(pt[1] > yy)[0].max()
    iz = np.where(pt[2] > zz)[0].max()

    x = (pt[0] - xx[ix])/(xx[ix+1] - xx[ix])
    y = (pt[1] - yy[iy])/(yy[iy+1] - yy[iy])
    z = (pt[2] - zz[iz])/(zz[iz+1] - zz[iz])

    btri = np.zeros(3, dtype=np.float64)
    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2, :]

    square = (1 - x)*cube[0, :, :, :] + x*cube[1, :, :, :]
    line = (1 - y)*square[0, :, :] + y*square[1, :, :]
    return (1 - z)*line[0, :] + z*line[1, :]

    # for ib in range(3):
    #   f11 = (1 - x)*cube[0,0,0,ib] + x*cube[1,0,0,ib]
    #   f12 = (1 - x)*cube[0,0,1,ib] + x*cube[1,0,1,ib]
    #   f21 = (1 - x)*cube[0,1,0,ib] + x*cube[1,1,0,ib]
    #   f22 = (1 - x)*cube[0,1,1,ib] + x*cube[1,1,1,ib]

    #   f1 = (1 - y)*f11 + y*f21
    #   f2 = (1 - y)*f12 + y*f22

    #   btri[ib] = (1 - z)*f1 + z*f2

    # return btri
def trilinearscalar3d(pt, grid, xx, yy, zz):
    ix = np.where(pt[0] > xx)[0].max()
    iy = np.where(pt[1] > yy)[0].max()
    iz = np.where(pt[2] > zz)[0].max()
    
    x = (pt[0] - xx[ix])/(xx[ix+1] - xx[ix])
    y = (pt[1] - yy[iy])/(yy[iy+1] - yy[iy])
    z = (pt[2] - zz[iz])/(zz[iz+1] - zz[iz])

    btri = np.zeros(3, dtype=np.float64)
    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2]

    square = (1 - x)*cube[0, :, :] + x*cube[1, :, :]
    line = (1 - y)*square[0, :] + y*square[1, :]
    return (1 - z)*line[0] + z*line[1]

def getdr(r, x, y, z):

    i = np.argwhere(r[0] > x).max()
    j = np.argwhere(r[1] > y).max()
    k = np.argwhere(r[2] > z).max()

    dx = x[i+1] - x[i]
    dy = y[j+1] - y[j]
    dz = z[k+1] - z[k]

    xh = x[i] + dx/2
    yh = y[j] + dy/2
    zh = z[k] + dz/2

    return np.array([dx, dy, dz])


def fieldline3d(startpt, bgrid, x, y, z, h, hmin, hmax, epsilon, mxline=10000, oneway=False, boxedge=None, gridcoord=False):

    # startpt[3,nl] - start point for field line
    # bgrid[nx,ny,nz,3] - magnetic field
    # x[nx],y[ny],z[nz] - grid of points on which magnetic field given

    # h - initial step length
    # hmin - minimum step length
    # hmax - maximum step length
    # epsilon - tolerance to which we require point on field line known

    # define edges of box
    if boxedge is not None:
        xmin = max([boxedge[0,0],x.min()])
        ymin = max([boxedge[0,1],y.min()])
        zmin = max([boxedge[0,2],z.min()])
        xmax = min([boxedge[1,0],x.max()])
        ymax = min([boxedge[1,1],y.max()])
        zmax = min([boxedge[1,2],z.max()])
    else:
        xmin = x.min()
        ymin = y.min()
        zmin = z.min()
        xmax = x.max()
        ymax = y.max()
        zmax = z.max()

    b2 = 0.25
    b3 = 3./32.
    c3 = 9./32.
    b4 = 1932./2197.
    c4 = -7200./2197.
    d4 = 7296./2197.
    b5 = 439./216.
    c5 = -8.
    d5 = 3680./513.
    e5 = -845./4104.
    b6 = -8./27.
    c6 = 2.
    d6 = -3544./2565.
    e6 = 1859./4104.
    f6 = -11./40.

    # used to determine y_i+1 from y_i if using rkf45 (4th order)
    n1 = 25./216.
    n3 = 1408./2565.
    n4 = 2197./4104.
    n5 = -1./5.

    # used to determine y_i+1 from y_i if using rkf54 (5th order)
    nn1 = 16./135.
    nn3 = 6656./12825.
    nn4 = 28561./56430.
    nn5 = -9./50.
    nn6 = 2./55.

    r0 = startpt

    if (r0[0] < xmin or r0[0] > xmax or r0[1] < ymin or r0[1] > ymax or r0[2] < zmin or r0[2] > zmax):
        print("Error: Start point not in range")
        print("Start point is: {} {} {}".format(r0[0],r0[1],r0[2]))
        if (r0[0] < xmin or r0[0] > xmax):
            print("{} (x) is the issue".format(r0[0]))
        if (r0[1] < ymin or r0[1] > ymax):
            print("{} (y) is the issue".format(r0[1]))
        if (r0[2] < zmin or r0[2] > zmax):
            print("{} (z) is the issue".format(r0[2]))
        print("Bounds are: xmin = {}, xmax = {}".format(xmin, xmax))
        print("Bounds are: ymin = {}, ymax = {}".format(ymin, ymax))
        print("Bounds are: zmin = {}, zmax = {}".format(zmin, zmax))

        return np.array([r0[0],r0[1],r0[2]], dtype=np.float64)

    ###################################################

    if oneway == True:
        ih = [h]
    else:
        ih = [h,-h]

    line = [r0]
    
    for h in ih:

        count = 0
        out = False
        bounce = 0
        t = 1.

        while count < mxline:

            if h > 0:
                jl = count
            else:
                jl = 0

            r0 = line[jl]

            if gridcoord == True:
                dr = getdr(r0, x, y, z)
                mindist = min(dr)*h
                hvec = mindist/dr
            else:
                hvec = np.array([h,h,h], dtype=np.float64)

            rt = r0
            b = trilinear3d(rt,bgrid,x,y,z)
            k1 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b2*k1

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k2 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b3*k1 + c3*k2

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k3 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b4*k1 + c4*k2 + d4*k3

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k4 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k5 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b6*k1 + c6*k2 + d6*k3 + e6*k4 + f6*k5

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k6 = hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)

            rtest4 = r0 + n1*k1 + n3*k3 + n4*k4 + n5*k5
            rtest5 = r0 + nn1*k1 + nn3*k3 + nn4*k4 + nn5*k5 + nn6*k6

            #check that line is still in domain
            if rtest4[0] < xmin or rtest4[0] > xmax or rtest4[1] < ymin or rtest4[1] > ymax or rtest4[2] < zmin or rtest4[2] > zmax:
                rout = rtest4
                out = True
                break

            #optimum stepsize
            diff = rtest5 - rtest4
            err = sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
            if err > 0:
                t = (epsilon*abs(h)/(2*err))**0.25 # do we want to update this

            # print (epsilon*abs(h)/(2*err))**0.25, abs(t*h), abs(hmin), hmin, h

            if (abs(t*h) < abs(hmin)):
                t = abs(hmin/h)

            t0 = 1.2
            if t > t0:
                t = t0

            rt = r0
            b = trilinear3d(rt,bgrid,x,y,z)
            k1 = t*hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b2*k1

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k2 = t*hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b3*k1 + c3*k2

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k3 = t*hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b4*k1 + c4*k2 + d4*k3

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k4 = t*hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            b = trilinear3d(rt,bgrid,x,y,z)
            k5 = t*hvec*b/sqrt(b[0]**2 + b[1]**2 + b[2]**2)
            # rf = r0 + b6*k1 + c6*k2 + d6*k3 + e6*k4 + f6*k5

            r4 = r0 + n1*k1 + n3*k3 + n4*k4 + n5*k5
            rt = r4.copy()
            if rt[0] < xmin or rt[0] > xmax or rt[1] < ymin or rt[1] > ymax or rt[2] < zmin or rt[2] > zmax:
                rout = rt
                out = True
                break

            h = t*h
            if abs(h) < hmin:
                h = hmin*h/abs(h)
            if abs(h) > hmax:
                h = hmax*h/abs(h)

            count += 1

            if h > 0:
                line = line + [r4]
            else:
                line = [r4] + line

            #check line is still moving
            if (count >= 2):
                if h > 0:
                    il = -1
                else:
                    il = 1
                dl = line[il] - line[il-1]
                mdl1 = sqrt(dl[0]**2 + dl[1]**2 + dl[2]**2)
                if h > 0:
                    il = -1
                else:
                    il = 2
                dl = line[il] - line[il-2]
                mdl2 = sqrt(dl[0]**2 + dl[1]**2 + dl[2]**2)
                if mdl1 < hmin*0.5:
                    break
                if mdl2 < hmin*0.5:
                    bounce = 1
                    break

        if bounce == 1:
            if il == -1:
                line = line[:-1]
            if il == 2:
                line = line[1:]

        if out:
            if h > 0:
                il = count-1
            else:
                il = 0
            rin = line[il]
            if rout[0] > xmax or rout[0] < xmin:
                if rout[0] > xmax:
                    xedge = xmax
                else:
                    xedge = xmin
                s = (xedge - rin[0])/(rout[0] - rin[0])
                rout = np.array([xedge, s*(rout[1] - rin[1]) + rin[1], s*(rout[2] - rin[2]) + rin[2]], dtype=np.float64)
            if rout[1] > ymax or rout[1] < ymin:
                if rout[1] > ymax:
                    yedge = ymax
                else:
                    yedge = ymin
                s = (yedge - rin[1])/(rout[1] - rin[1])
                rout = np.array([s*(rout[0] - rin[0]) + rin[0], yedge, s*(rout[2] - rin[2]) + rin[2]], dtype=np.float64)
            if rout[2] > zmax or rout[2] < zmin:
                if rout[2] > zmax:
                    zedge = zmax
                else:
                    zedge = zmin
                s = (zedge - rin[2])/(rout[2] - rin[2])
                rout = np.array([s*(rout[0] - rin[0]) + rin[0], s*(rout[1] - rin[1]) + rin[1], zedge], dtype=np.float64)
            if h > 0:
                line = line + [rout]
            else:
                line = [rout] + line
    
    return np.array(line, dtype=np.float64)
