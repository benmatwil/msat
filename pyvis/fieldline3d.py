import numpy as np
from math import *

# rkf45 coefficients
b2 = 0.25
b3, c3 = 3.0/32, 9.0/32
b4, c4, d4 = 1932.0/2197, -7200.0/2197, 7296.0/2197
b5, c5, d5, e5 = 439.0/216, -8.0, 3680.0/513, -845.0/4104
b6, c6, d6, e6, f6 = -8.0/27, 2.0, -3544.0/2565, 1859.0/4104, -11.0/40

# used to determine y_i+1 from y_i if using rkf45 (4th order)
n1, n3, n4, n5 = 25.0/216, 1408.0/2565, 2197.0/4104, -1.0/5

# used to determine y_i+1 from y_i if using rkf54 (5th order)
nn1, nn3, nn4, nn5, nn6 = 16.0/135, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55

def trilinear3d(pt, grid, xx, yy, zz):

    ix = np.where(pt[0] > xx)[0].max()
    iy = np.where(pt[1] > yy)[0].max()
    iz = np.where(pt[2] > zz)[0].max()

    x = (pt[0] - xx[ix])/(xx[ix+1] - xx[ix])
    y = (pt[1] - yy[iy])/(yy[iy+1] - yy[iy])
    z = (pt[2] - zz[iz])/(zz[iz+1] - zz[iz])

    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2, ...]
    square = (1 - x)*cube[0, :, :, ...] + x*cube[1, :, :, ...]
    line = (1 - y)*square[0, :, ...] + y*square[1, :, ...]
    return (1 - z)*line[0, ...] + z*line[1, ...]

def trilinear3d_grid(pt, grid):
    ix, iy, iz = np.floor(pt).astype(np.int)
    x, y, z = pt - np.floor(pt)

    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2, ...]
    square = (1 - x)*cube[0, :, :, ...] + x*cube[1, :, :, ...]
    line = (1 - y)*square[0, :, ...] + y*square[1, :, ...]
    return (1 - z)*line[0, ...] + z*line[1, ...]

def getdr(r, x, y, z):
    i, j, k = np.floor(r).astype(np.int)

    dx = x[i+1] - x[i]
    dy = y[j+1] - y[j]
    dz = z[k+1] - z[k]

    xc = x[i] + dx/2
    yc = y[j] + dy/2

    if csystem == 'spherical':
        return np.array([dx, xc*dy, xc*np.sin(yc)*dz], dtype=np.float64)
    elif csystem == 'cartesian':
        return np.array([dx, dy, dz], dtype=np.float64)
    elif cystem == 'cylindrical':
        return np.array([dx, xc*dy, dz], dtype=np.float64)

def edgecheck(r):
    if csystem == 'spherical':
        if r[1] < ymin or r[1] > ymax:
            if r[1] < ymin: r[1] = 2*ymin - r[1]
            if r[1] > ymax: r[1] = 2*ymax - r[1]
            if r[2] < (zmax + zmin)/2:
                r[2] = r[2] + (zmax - zmin)/2
            else:
                r[2] = r[2] - (zmax - zmin)/2
        if r[2] <= zmin: r[2] = r[2] + (zmax - zmin)
        if r[2] >= zmax: r[2] = r[2] - (zmax - zmin)
    elif csystem == 'cylindrical':
        if r[2] < ymin: r[2] = r[2] + (ymax - ymin)
        if r[2] > ymax: r[2] = r[2] - (ymax - ymin)

def outedge(r):
    outedge = r[0] >= xmax or r[0] <= xmin
    if csystem == 'cartesian':
        outedge = outedge or r[1] >= ymax or r[1] <= ymin or r[2] >= zmax or r[2] <= zmin
    elif csystem == 'cylindrical':
        outedge = outedge or r[2] >= zmax or r[2] <= zmin
    
    return outedge

def fieldline3d(startpt, bgrid, x, y, z, h, hmin, hmax, epsilon, mxline=50000, t_max=1.2, oneway=False, boxedge=None, coordsystem='cartesian', gridcoord=False):
    global xmin, xmax, ymin, ymax, zmin, zmax, csystem
    # need to correct this below
    # startpt[3,nl] - start point for field line
    # bgrid[nx,ny,nz,3] - magnetic field
    # x[nx],y[ny],z[nz] - grid of points on which magnetic field given

    # h - initial step length
    # hmin - minimum step length
    # hmax - maximum step length
    # epsilon - tolerance to which we require point on field line known
    
    csystem = coordsystem

    # define edges of box
    if boxedge is not None:
        xmin = max([boxedge[0,0],x.min()])
        ymin = max([boxedge[0,1],y.min()])
        zmin = max([boxedge[0,2],z.min()])
        xmax = min([boxedge[1,0],x.max()])
        ymax = min([boxedge[1,1],y.max()])
        zmax = min([boxedge[1,2],z.max()])
    else:
        xmin = 0
        ymin = 0
        zmin = 0
        xmax = x.shape[0]-1
        ymax = y.shape[0]-1
        zmax = z.shape[0]-1

    # first convert point into grid coordinates
    ix = np.argwhere(startpt[0] >= x).max()
    iy = np.argwhere(startpt[1] >= y).max()
    iz = np.argwhere(startpt[2] >= z).max()

    r0 = np.empty_like(startpt)
    r0[0] = ix + (startpt[0] - x[ix])/(x[ix+1] - x[ix])
    r0[1] = iy + (startpt[1] - y[iy])/(y[iy+1] - y[iy])
    r0[2] = iz + (startpt[2] - z[iz])/(z[iz+1] - z[iz])

    # Produce an error if the first point isn't in the box
    if startpt[0] < x[0] or startpt[0] > x[-1] or startpt[1] < y[0] or startpt[1] > y[-1] or startpt[2] < z[0] or startpt[2] > z[-1]:
        print("Error: Start point not in range")
        print("Start point is: {} {} {}".format(startpt[0],startpt[1],startpt[2]))
        if (startpt[0] < x[0] or startpt[0] > x[-1]):
            print("{} (x) is the issue".format(startpt[0]))
        if (startpt[1] < y[0] or startpt[1] > y[-1]):
            print("{} (y) is the issue".format(startpt[1]))
        if (startpt[2] < z[0] or startpt[2] > z[-1]):
            print("{} (z) is the issue".format(startpt[2]))
        print("Bounds are: x[0] = {}, x[-1] = {}".format(x[0], x[-1]))
        print("Bounds are: y[0] = {}, y[-1] = {}".format(y[0], y[-1]))
        print("Bounds are: z[0] = {}, z[-1] = {}".format(z[0], z[-1]))

        return np.array([startpt[0],startpt[1],startpt[2]], dtype=np.float64)
    elif not (hmin < h < hmax):
        print('You need to satisfy hmin ({}) < h ({}) < hmax({})'.format(hmin, h, hmax))
        return

    ###################################################

    ih = [h] if oneway else [h, -h]

    line = [r0]
    
    for h in ih:

        count = 0
        out = False
        bounce = 0
        line = line[::-1]

        while count < mxline:

            r0 = line[-1].copy()

            dr = getdr(r0, x, y, z)
            mindist = min(dr)*h
            hvec = mindist/dr

            rt = r0
            b = trilinear3d_grid(rt, bgrid)
            k1 = hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b2*k1
            
            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k2 = hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b3*k1 + c3*k2

            if outedge(rt):
                out = True
                break
            
            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k3 = hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b4*k1 + c4*k2 + d4*k3

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k4 = hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k5 = hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b6*k1 + c6*k2 + d6*k3 + e6*k4 + f6*k5

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k6 = hvec*b/sqrt(np.sum(b**2))

            # 4th order estimate
            rtest4 = r0 + n1*k1 + n3*k3 + n4*k4 + n5*k5
            # 5th order estimate
            rtest5 = r0 + nn1*k1 + nn3*k3 + nn4*k4 + nn5*k5 + nn6*k6

            # optimum stepsize
            diff = rtest5 - rtest4
            err = sqrt(np.sum(diff**2))
            if err > 0:
                t = (epsilon*abs(h)/(2*err))**0.25
            else:
                t = 1

            if (abs(t*h) < abs(hmin)): t = abs(hmin/h)
            if t > t_max: t = t_max

            rt = r0
            b = trilinear3d_grid(rt, bgrid)
            k1 = t*hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b2*k1

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k2 = t*hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b3*k1 + c3*k2

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k3 = t*hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b4*k1 + c4*k2 + d4*k3

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k4 = t*hvec*b/sqrt(np.sum(b**2))
            rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

            if outedge(rt):
                out = True
                break

            edgecheck(rt)
            b = trilinear3d_grid(rt, bgrid)
            k5 = t*hvec*b/sqrt(np.sum(b**2))

            rt = r0 + n1*k1 + n3*k3 + n4*k4 + n5*k5
            edgecheck(rt)
            
            if outedge(rt):
                out = True
                break

            h = t*h
            if abs(h) < hmin: h = hmin*h/abs(h)
            if abs(h) > hmax: h = hmax*h/abs(h)

            count += 1

            line.append(rt)

            # check line is still moving
            if count >= 2:
                dl = line[-1] - line[-2]
                mdl = sqrt(np.sum(dl**2))
                if mdl < hmin*0.5:
                    break

                dl = line[-1] - line[-3]
                mdl = sqrt(np.sum(dl**2))
                if mdl < hmin*0.5:
                    bounce = 1
                    break

        if out:
            rout = rt.copy()
            rin = line[-1].copy()
            if rout[0] > xmax or rout[0] < xmin:
                xedge = xmax if rout[0] > xmax else xmin
                s = (xedge - rin[0])/(rout[0] - rin[0])
                rout = s*(rout - rin) + rin
            if rout[1] > ymax or rout[1] < ymin:
                yedge = ymax if rout[1] > ymax else ymin
                s = (yedge - rin[1])/(rout[1] - rin[1])
                rout = s*(rout - rin) + rin
            if rout[2] > zmax or rout[2] < zmin:
                zedge = zmax if rout[2] > zmax else zmin
                s = (zedge - rin[2])/(rout[2] - rin[2])
                rout = s*(rout - rin) + rin
            line.append(rout)
        elif bounce == 1: line = line[:-1]

    if gridcoord == False:
        for pt in line:
            ix, iy, iz = np.floor(pt).astype(np.int)
            if ix == xmax: ix -= 1
            if iy == ymax: iy -= 1
            if iz == zmax: iz -= 1

            pt[0] = x[ix] + (pt[0] - ix)*(x[ix+1] - x[ix])
            pt[1] = y[iy] + (pt[1] - iy)*(y[iy+1] - y[iy])
            pt[2] = z[iz] + (pt[2] - iz)*(z[iz+1] - z[iz])
    
    return np.array(line[::-1], dtype=np.float64)