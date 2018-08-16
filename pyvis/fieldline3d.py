from __future__ import print_function, division
import numpy as np
from math import sqrt, floor
from numba import njit, float64 as f64, void, boolean

# rkf45 coefficients
b2 = 0.25
b3, c3 = 3/32, 9/32
b4, c4, d4 = 1932/2197, -7200/2197, 7296/2197
b5, c5, d5, e5 = 439/216, -8, 3680/513, -845/4104
b6, c6, d6, e6, f6 = -8/27, 2, -3544/2565, 1859/4104, -11/40

# used to determine y_i+1 from y_i if using rkf45 (4th order)
n1, n3, n4, n5 = 25/216, 1408/2565, 2197/4104, -1/5

# used to determine y_i+1 from y_i if using rkf54 (5th order)
nn1, nn3, nn4, nn5, nn6 = 16/135, 6656/12825, 28561/56430, -9/50, 2/55

# @njit(f64[:](f64[:], f64[:, :, :, :], f64[:], f64[:], f64[:]), cache=True)
@njit
def trilinear3d(pt, grid, xx, yy, zz):
    ix = np.where(pt[0] > xx)[0][-1]
    iy = np.where(pt[1] > yy)[0][-1]
    iz = np.where(pt[2] > zz)[0][-1]

    x = (pt[0] - xx[ix])/(xx[ix+1] - xx[ix])
    y = (pt[1] - yy[iy])/(yy[iy+1] - yy[iy])
    z = (pt[2] - zz[iz])/(zz[iz+1] - zz[iz])

    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2, ...]
    square = (1 - x)*cube[0, :, :, ...] + x*cube[1, :, :, ...]
    line = (1 - y)*square[0, :, ...] + y*square[1, :, ...]
    return (1 - z)*line[0, ...] + z*line[1, ...]

# @njit(f64[:](f64[:], f64[:, :, :, :]), cache=True)
@njit
def trilinear3d_grid(pt, grid):
    ix = floor(pt[0])
    iy = floor(pt[1])
    iz = floor(pt[2])

    x = pt[0] - ix
    y = pt[1] - iy
    z = pt[2] - iz

    cube = grid[ix:ix+2, iy:iy+2, iz:iz+2, ...]
    square = (1 - x)*cube[0, :, :, ...] + x*cube[1, :, :, ...]
    line = (1 - y)*square[0, :, ...] + y*square[1, :, ...]
    return (1 - z)*line[0, ...] + z*line[1, ...]

# @njit(f64[:](f64[:], f64[:], f64[:], f64[:], boolean[:]), cache=True)
@njit
def getdr(r, x, y, z, csystem):
    ix = floor(r[0])
    iy = floor(r[1])
    iz = floor(r[2])

    dx = x[ix+1] - x[ix]
    dy = y[iy+1] - y[iy]
    dz = z[iz+1] - z[iz]

    xp = x[ix] + (r[0] - ix)*dx
    yp = y[iy] + (r[1] - iy)*dy

    if csystem[2]:
        dr = np.array([dx, xp*dy, xp*np.sin(yp)*dz], dtype=np.float64)
    elif csystem[0]:
        dr = np.array([dx, dy, dz], dtype=np.float64)
    elif csystem[1]:
        dr = np.array([dx, xp*dy, dz], dtype=np.float64)
    return dr

# @njit(void(f64[:], f64[:], f64[:], f64[:]), cache=True)
@njit
def gtr(pt, x, y, z):
    ix = floor(pt[0])
    iy = floor(pt[1])
    iz = floor(pt[2])

    if ix == x.shape[0]-1: ix -= 1
    if iy == y.shape[0]-1: iy -= 1
    if iz == z.shape[0]-1: iz -= 1

    pt[0] = x[ix] + (pt[0] - ix)*(x[ix+1] - x[ix])
    pt[1] = y[iy] + (pt[1] - iy)*(y[iy+1] - y[iy])
    pt[2] = z[iz] + (pt[2] - iz)*(z[iz+1] - z[iz])

# @njit(void(f64[:], f64[:], boolean[:], boolean[:]), cache=True)
@njit
def edgecheck(r, minmax, csystem, periodicity):
    # need to implement periodic_t, periodic_p
    if np.any(periodicity):
        if csystem[2]:
            if periodicity[3]:
                if r[1] < minmax[2] or r[1] > minmax[3]:
                    if r[1] < minmax[2]: r[1] = 2*minmax[2] - r[1]
                    if r[1] > minmax[3]: r[1] = 2*minmax[3] - r[1]
                    if r[2] < (minmax[5] + minmax[4])/2:
                        r[2] = r[2] + (minmax[5] - minmax[4])/2
                    else:
                        r[2] = r[2] - (minmax[5] - minmax[4])/2
            if periodicity[4]:
                if r[2] <= minmax[4]: r[2] = r[2] + (minmax[5] - minmax[4])
                if r[2] >= minmax[5]: r[2] = r[2] - (minmax[5] - minmax[4])
        elif csystem[0]:
            if periodicity[0]:
                if r[0] <= minmax[0]: r[0] = r[0] + (minmax[1] - minmax[0])
                if r[0] >= minmax[1]: r[0] = r[0] - (minmax[1] - minmax[0])
            if periodicity[1]:
                if r[1] <= minmax[2]: r[1] = r[1] + (minmax[3] - minmax[2])
                if r[1] >= minmax[3]: r[1] = r[1] - (minmax[3] - minmax[2])
            if periodicity[2]:
                if r[2] <= minmax[4]: r[2] = r[2] + (minmax[5] - minmax[4])
                if r[2] >= minmax[5]: r[2] = r[2] - (minmax[5] - minmax[4])
        elif csystem[1]:
            if periodicity[4]:
                if r[2] < minmax[2]: r[2] = r[2] + (minmax[3] - minmax[2])
                if r[2] > minmax[3]: r[2] = r[2] - (minmax[3] - minmax[2])

# @njit(boolean(f64[:], f64[:], boolean[:], boolean[:]), cache=True)
@njit
def outedge(r, minmax_box, csystem, periodicity):
    # need to implement periodic_t, periodic_p
    if np.any(periodicity):
        if csystem[0]:
            outedge = False
            if not periodicity[0]:
                outedge = outedge or r[0] >= minmax_box[1] or r[0] <= minmax_box[0]
            if not periodicity[1]:
                outedge = outedge or r[1] >= minmax_box[3] or r[1] <= minmax_box[2]
            if not periodicity[2]:
                outedge = outedge or r[2] >= minmax_box[5] or r[2] <= minmax_box[4]
        elif csystem[2]:
            outedge = r[0] >= minmax_box[1] or r[0] <= minmax_box[0]
            if not periodicity[3]:
                outedge = outedge or r[1] >= minmax_box[3] or r[1] <= minmax_box[2]
            if not periodicity[4]:
                outedge = outedge or r[2] >= minmax_box[5] or r[2] <= minmax_box[4]
        elif csystem[1]:
            outedge = r[0] >= minmax_box[1] or r[0] <= minmax_box[0] or r[2] >= minmax_box[5] or r[2] <= minmax_box[4]
            if not periodicity[4]:
                outedge = outedge or r[1] >= minmax_box[3] or r[1] <= minmax_box[2]
    else:
        outedge = r[0] >= minmax_box[1] or r[0] <= minmax_box[0] or r[1] >= minmax_box[3] or r[1] <= minmax_box[2] or r[2] >= minmax_box[5] or r[2] <= minmax_box[4]
        
    return outedge

# @njit(cache=True)
@njit
def rkf45(r0, bgrid, x, y, z, h, hmin, hmax, epsilon, t_max, minmax, minmax_box, csystem, periodicity):
    dr = getdr(r0, x, y, z, csystem)
    mindist = dr.min()*h
    # hvec = mindist/dr
    hvec = np.zeros(3, dtype=np.float64)
    hvec[0] = mindist/dr[0]
    hvec[1] = mindist/dr[1]
    hvec[2] = mindist/dr[2]

    rt = r0
    b = trilinear3d_grid(rt, bgrid)
    k1 = hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b2*k1
    
    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k2 = hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b3*k1 + c3*k2

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h
    
    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k3 = hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b4*k1 + c4*k2 + d4*k3

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k4 = hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k5 = hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b6*k1 + c6*k2 + d6*k3 + e6*k4 + f6*k5

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
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

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k2 = t*hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b3*k1 + c3*k2

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k3 = t*hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b4*k1 + c4*k2 + d4*k3

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k4 = t*hvec*b/sqrt(np.sum(b**2))
    rt = r0 + b5*k1 + c5*k2 + d5*k3 + e5*k4

    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    edgecheck(rt, minmax, csystem, periodicity)
    b = trilinear3d_grid(rt, bgrid)
    k5 = t*hvec*b/sqrt(np.sum(b**2))
    rt = r0 + n1*k1 + n3*k3 + n4*k4 + n5*k5
    edgecheck(rt, minmax, csystem, periodicity)
    
    if outedge(rt, minmax_box, csystem, periodicity): return rt, True, h

    h = t*h
    if abs(h) < hmin: h = hmin*h/abs(h)
    if abs(h) > hmax: h = hmax*h/abs(h)
    
    return rt, False, h

def fieldline3d(startpt, bgrid, x, y, z, h, hmin, hmax, epsilon, mxline=50000, t_max=1.2, oneway=False,
    boxedge=None, coordsystem='cartesian', gridcoord=False, stop_criteria=True, periodicity=''):
    """
    Calculates 3D field line which goes through the point startpt
    startpt - 3 element, 1D array as start point for field line calculation
    bgrid - magnetic field array of shape (nx, ny, nz, 3)
    x, y, z - 1D arrays of grid points on which magnetic field given

    h - initial step length
    hmin - minimum step length
    hmax - maximum step length
    epsilon - tolerance to which we require point on field line known

    maxline - maximum number of points on a fieldline in each direction
    t_max - maximum value of correction factor in RKF45 method
    oneway - whether to only calculate field in one direction (sign of h)
    boxedge - use a smaller domain than edge of the grids x, y and z
    coordsystem - set the coordinate system of the magnetic field grid
    gridcoord - use grid coordinates as input and output to routine
    stop_critieria - use the stopping criteria to detect if field line has stopped inside domain
    periodicity - set the periodicity of the grid e.g. 'xy' for periodicity in both x and y direction
    """
    csystem = np.array([
        True if coordsystem == 'cartesian' else False,
        True if coordsystem == 'cylindrical' else False,
        True if coordsystem == 'spherical' else False
        ], dtype=np.bool_)
    
    periodicity = np.array([
        True if 'x' in periodicity else False,
        True if 'y' in periodicity else False,
        True if 'z' in periodicity else False,
        False if 't' in periodicity else True,
        False if 'p' in periodicity else True
        ], dtype=np.bool_)
    
    # define edges of box
    minmax = np.array([0, x.shape[0]-1, 0, y.shape[0]-1, 0, z.shape[0]-1], dtype=np.float64)
    if boxedge is not None:
        xmin_box = max([boxedge[0,0], x.min()])
        ymin_box = max([boxedge[0,1], y.min()])
        zmin_box = max([boxedge[0,2], z.min()])
        xmax_box = min([boxedge[1,0], x.max()])
        ymax_box = min([boxedge[1,1], y.max()])
        zmax_box = min([boxedge[1,2], z.max()])
        minmax_box = np.array([xmin_box, xmax_box, ymin_box, ymax_box, zmin_box, zmax_box], dtype=np.float64)
    else:
        minmax_box = minmax

    # first convert point into grid coordinates
    if not gridcoord:
        ix = np.argwhere(startpt[0] >= x).max()
        iy = np.argwhere(startpt[1] >= y).max()
        iz = np.argwhere(startpt[2] >= z).max()

        r0 = np.empty_like(startpt)
        r0[0] = ix + (startpt[0] - x[ix])/(x[ix+1] - x[ix])
        r0[1] = iy + (startpt[1] - y[iy])/(y[iy+1] - y[iy])
        r0[2] = iz + (startpt[2] - z[iz])/(z[iz+1] - z[iz])
    else:
        r0 = startpt.copy()

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
        raise ValueError
    elif not (hmin < np.abs(h) < hmax):
        print('You need to satisfy hmin ({}) < h ({}) < hmax({})'.format(hmin, h, hmax))
        raise ValueError
    elif np.all(trilinear3d_grid(r0, bgrid) == np.zeros(3)):
        print('Start point is a null point')
        raise ValueError

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

            rt, out, h = rkf45(r0, bgrid, x, y, z, h, hmin, hmax, epsilon, t_max, minmax, minmax_box, csystem, periodicity)

            if out: break

            count += 1

            line.append(rt)

            if stop_criteria:
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
            if rout[0] > minmax[1] or rout[0] < minmax[0]:
                xedge = minmax[1] if rout[0] > minmax[1] else minmax[0]
                s = (xedge - rin[0])/(rout[0] - rin[0])
                rout = s*(rout - rin) + rin
            if rout[1] > minmax[3] or rout[1] < minmax[2]:
                yedge = minmax[3] if rout[1] > minmax[3] else minmax[2]
                s = (yedge - rin[1])/(rout[1] - rin[1])
                rout = s*(rout - rin) + rin
            if rout[2] > minmax[5] or rout[2] < minmax[4]:
                zedge = minmax[5] if rout[2] > minmax[5] else minmax[4]
                s = (zedge - rin[2])/(rout[2] - rin[2])
                rout = s*(rout - rin) + rin
            line.append(rout)
        elif bounce == 1: line = line[:-1]

    if gridcoord == False:
        for pt in line: gtr(pt, x, y, z)
    
    return np.array(line[::-1], dtype=np.float64)
