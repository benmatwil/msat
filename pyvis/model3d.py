from __future__ import print_function, division
import numpy as np
import mayavi.mlab as ml
from . import read as rd
from . import fieldline3d as fl
import sys
import vtk

try:
    __IPYTHON__
except NameError:
    pass
else:
    from IPython import get_ipython
    ipython = get_ipython()
    ipython.magic(r"%gui qt")

# turn of warnings while vtk/mayavi compatibility is fixed -- still works
vtk.vtkObject.GlobalWarningDisplayOff()

filename = None

sign_names = {-2:'Sink', -1:'Neg', 0:'Zero', 1:'Pos', 2:'Source'}

def make(fname, addlist, null_list=None, box=True, fieldlines=None, linecolor=(0,0,0), nskip=20,
    nullrad=1, nfanlines=40, nring=None, colquant=None, coordsystem='cartesian', no_nulls=False,
    sun=True, outdir=None, periodicity=''):
    """Makes a 3D visualisation of the output from Magnetic Skeleton Analysis Tools

        fname: name of the file containing the original magnetic field
        addlist: list of features to be plotted e.g. ['nulls', 'separators'] will
            only plot the nulls and separators

        nulls: if None (default), will plot all nulls, otherwise give a list (starting at 1) of nulls to plot
            e.g. nulls=list(range(45,65)) for all nulls between 45 and 64 inclusive
        nullrad: will scale radius of null spheres by this factor (default 1)
            e.g. nullrad=0.5 will halve the size of the nulls
        box: if True (default), plots a box otherwise set to False
        fieldlines: provide a numpy array of shape (3, n) of start points and it will trace magnetic field lines
        linecolor: color of fieldlines (defaults to black)
        nskip: how many rings to skip in plotting"""

    global bgrid, xx, yy, zz, nulldata, ds, filename, nskipglob, nulllist, csystem, periodic_check, periodic_dist, periodic_global

    csystem = coordsystem
    print('Using {} coordinates'.format(coordsystem))

    ml.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 800))

    nskipglob = nskip

    if outdir is not None:
        rd.outprefix = outdir

    if fname != filename:
        field = rd.field(fname)
        bgrid = np.zeros(field[0].shape + (3,), dtype=np.float64)
        for i in range(3):
            bgrid[:, :, :, i] = field[i]
        xx = field[3]
        yy = field[4]
        zz = field[5]
        ds = min([np.diff(xx).min(), np.diff(yy).min(), np.diff(zz).min()])/2

    periodic_check = False
    if periodicity != '':
        periodic_dist = 1e100
        if 'x' in periodicity:
            periodic_check = True
            periodic_dist = min((periodic_dist, (xx[-1] - xx[0])**2))
        if 'y' in periodicity:
            periodic_check = True
            periodic_dist = min((periodic_dist, (yy[-1] - yy[0])**2))
        if 'z' in periodicity:
            periodic_check = True
            periodic_dist = min((periodic_dist, (zz[-1] - zz[0])**2))
    periodic_global = periodicity

    filename = fname

    if not no_nulls:
        nulldata = rd.nulls(filename)
        set_null_list(null_list)

    add_structures(*addlist, nullrad=nullrad, nfanlines=nfanlines, nring=nring)

    if coordsystem == 'spherical':
        if sun: add_sun()
        box = False
    if box: add_box()

    if fieldlines is not None: add_fieldlines(fieldlines, col=linecolor)

def set_null_list(lst):
    global nulllist
    if lst is None:
        nulllist = nulldata.number
    else:
        if max(lst) > nulldata.number.max() or min(lst) < 1:
            print('Invalid list')
        nulllist = np.array(lst)

def add_structures(*args, **kwargs):
    if 'nulls' in args: add_nulls(size=kwargs['nullrad'])
    if 'sepsurf_flines' in args: add_sepsurf_flines(kwargs['nfanlines'], kwargs['nring'])
    if 'spines' in args: add_spines()
    if 'sepsurf_rings' in args: add_sepsurf_rings()
    if 'separators' in args: add_separators()
    if 'hcs_rings' in args: add_hcs_rings()
    if 'hcs_flines' in args: add_hcs_flines()
    if 'hcs_sep' in args: add_separators(hcs=True)

def add_sepsurf_rings():
    print('Adding separatrix surface rings')

    rings, breaks = rd.rings(filename, breaks=True, nskip=nskipglob, null_list=nulllist)

    cols = {-2:(218/255, 112/255, 214/255), -1:(0.5, 0.5, 1), 0:(0.5, 1, 0.5), 1:(1, 0.5, 0.5), 2:(1.0, 178/255, 102/255)}

    nulls = nulldata[nulllist-1]

    for isign in np.unique(nulls.sign):
        # new very efficient routine for plotting many lines
        # two lists, one for positive and the other for negative nulls
        x, y, z, s, ptcons = ( [] for _ in range(5) )
        index = 0
        for inull in nulls.number[nulls.sign == isign]-1:
            print('Null {:5d}'.format(inull+1))
            sys.stdout.write("\033[F")
            for iring, ring in enumerate(rings[inull]):
                # convert points if required
                if csystem == 'spherical':
                    ring[:, 0], ring[:, 1], ring[:, 2] = sphr2cart(ring[:, 0], ring[:, 1], ring[:, 2])
                elif csystem == 'cylindrical':
                    ring[:, 0], ring[:, 1], ring[:, 2] = cyl2cart(ring[:, 0], ring[:, 1], ring[:, 2])
                # add ring points to lists
                x.append(ring[:, 0])
                y.append(ring[:, 1])
                z.append(ring[:, 2])
                s.append(np.zeros_like(ring[:, 0]))
                if periodic_check:
                    # use distances between consectutive points to detect extra breaks for periodicity
                    dists = np.r_[np.sum(np.diff(ring, axis=0)**2, axis=1), [0]]
                    # use break data to plot the individual lines in each ring as the break apart
                    brks = np.unique(np.r_[-1, np.where(breaks[inull][iring] == 1)[0], np.where(dists > 0.9*periodic_dist)[0], ring.shape[0]-1])
                else:
                    # use break data to plot the individual lines in each ring as the break apart
                    brks = np.unique(np.r_[-1, np.where(breaks[inull][iring] == 1)[0], ring.shape[0]-1])
                for ib0, ib1 in zip(brks[:-1], brks[1:]):
                    # add the right indicies based on the breaks
                    ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
                index += ring.shape[0]
        # and plot...
        if len(x) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
            src.mlab_source.dataset.lines = np.vstack(ptcons)
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=1, name=sign_names[isign]+'SeparatrixRings')

def add_hcs_rings():
    print('Adding heliospheric current sheet curtain surface rings')

    rings, breaks = rd.rings(filename, breaks=True, nskip=nskipglob, hcs=True)

    x, y, z, s, ptcons = ( [] for _ in range(5) )
    index = 0

    for inull in range(len(rings)):
        print('HCS {:5d}'.format(inull//2+1))
        sys.stdout.write("\033[F")
        for iring, ring in enumerate(rings[inull]):
            # convert points if required
            if csystem == 'spherical':
                ring[:, 0], ring[:, 1], ring[:, 2] = sphr2cart(ring[:, 0], ring[:, 1], ring[:, 2])
            # add ring points to lists
            x.append(ring[:, 0])
            y.append(ring[:, 1])
            z.append(ring[:, 2])
            s.append(np.zeros_like(ring[:, 0]))
            # use break data to plot the individual lines in each ring as the break apart
            brks = np.unique(np.r_[[-1], np.where(breaks[inull][iring] == 1)[0], [ring.shape[0]-1]])
            for ib0, ib1 in zip(brks[:-1], brks[1:]):
                # add the right indicies based on the breaks
                ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
            index += ring.shape[0]
    
    # add points to model
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=(0, 1, 0), line_width=1, name='HCSRings')
    
    for inull in range(0, len(rings), 2):
        ml.plot3d(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2], color=(0, 1, 0), line_width=6, tube_radius=None, name='HCSBase')

def add_sepsurf_flines(nlines, nring=None):
    print('Adding separatrix surface field lines')

    rings = rd.rings(filename, nskip=nskipglob, null_list=nulllist)

    cols = {-2:(218/255, 112/255, 214/255), -1:(0.5, 0.5, 1), 0:(0.5, 1, 0.5), 1:(1, 0.5, 0.5), 2:(1.0, 178/255, 102/255)}

    nulls = nulldata[nulllist-1]

    for isign in np.unique(nulls.sign):
        x, y, z, s, ptcons = ( [] for _ in range(5) )
        index = 0
        for inull in nulls.number[nulls.sign == isign]-1:
            print('Null {}'.format(inull+1))
            sys.stdout.write("\033[F")

            # select the right ring to trace from
            if nring is None:
                ring = rings[inull][len(rings[inull])//5]
            else:
                if type(nring) is int:
                    ring = rings[inull][nring]
                elif type(nring) is float:
                    ring = rings[inull][np.floor(nring*len(rings[inull]))]

            nskip = len(ring[:, 0])//nlines

            for ipt, startpt in enumerate(ring[::nskip, :], start=1):
                # choose some good parameters
                h = 1e-2
                hmin = 1e-3
                hmax = 0.5
                epsilon = 1e-5

                # calculate the fieldline
                line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem, periodicity=periodic_global)

                # cut off the fieldline at the point closest to the null - only want the fan, not the spine
                dists = np.sqrt((line[:, 0] - nulldata[inull].pos[0])**2 +
                    (line[:, 1] - nulldata[inull].pos[1])**2 + (line[:, 2] - nulldata[inull].pos[2])**2)
                imin = dists.argmin()
                line = line[0:imin+1, :] if nulldata[inull].sign == -1 else line[imin:, :]

                if csystem == 'spherical':
                    line[:, 0], line[:, 1], line[:, 2] = sphr2cart(line[:, 0], line[:, 1], line[:, 2])
                elif csystem == 'cylindrical':
                    line[:, 0], line[:, 1], line[:, 2] = cyl2cart(line[:, 0], line[:, 1], line[:, 2])
                
                x.append(line[:, 0])
                y.append(line[:, 1])
                z.append(line[:, 2])
                length = len(line[:, 0])
                s.append(np.zeros(length))
                if csystem == 'spherical' or csystem == 'cylindical' or not periodic_check:
                    ptcons.append(np.vstack([np.arange(index, index+length-1), np.arange(index+1, index+length)]).T)
                else:
                    dists = np.r_[np.sum(np.diff(line, axis=0)**2, axis=1), 0]
                    brks = np.unique(np.r_[-1, np.where(dists > 0.9*periodic_dist)[0], dists.shape[0]-1])
                    for ib0, ib1 in zip(brks[:-1], brks[1:]):
                        ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
                index += length
    
        # add points to model
        if len(x) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
            src.mlab_source.dataset.lines = np.vstack(ptcons)
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=1, name=sign_names[isign]+'SeparatrixFieldlines')

def add_hcs_flines(nlines=100):
    print('Adding heliospheric current sheet curtain surface field lines')

    rings = rd.rings(filename, nskip=nskipglob, hcs=True)

    x, y, z, s, ptcons = ( [] for _ in range(5) )
    index = 0

    for inull in range(0, len(rings), 2):
        print('HCS {:5d}'.format(inull//2+1))
        sys.stdout.write("\033[F")

        iring = 1
        nskip = len(rings[inull][iring][:, 0])//nlines

        for idir in range(2):
            for startpt in rings[inull+idir][iring][::nskip, :]:
                # choose some good parameters
                h = 2e-2
                hmin = h*0.1
                hmax = h*10
                epsilon = h*0.01

                # calculate the fieldline
                line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem)
                imax = np.argmax(line[:, 0])
                if idir == 1:
                    line = line[:imax+1, :]
                else:
                    line = line[imax:, :]

                line[:, 0], line[:, 1], line[:, 2] = sphr2cart(line[:, 0], line[:, 1], line[:, 2])

                x.append(line[:, 0])
                y.append(line[:, 1])
                z.append(line[:, 2])
                length = len(line[:, 0])
                s.append(np.zeros(length))
                ptcons.append(np.vstack([np.arange(index, index+length-1), np.arange(index+1, index+length)]).T)
                index += length
        
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=(0, 1, 0), line_width=1, name='HCSFieldlines')
                
    for inull in range(0, len(rings), 2):
        rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2] = sphr2cart(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2])
        ml.plot3d(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2], color=(0, 1, 0), line_width=6, tube_radius=None, name='HCSBase')

def add_fieldlines(startpts, col=(0, 0, 0), lw=2):
    print('Adding separatrix surface field lines')

    x, y, z, s, ptcons = ( [] for _ in range(5) )
    index = 0

    for i, startpt in enumerate(startpts, start=1):
        print('Calculating fieldline {}'.format(i))
        sys.stdout.write("\033[F")

        # choose fieldline parameters
        h = 1e-2
        hmin = 1e-3
        hmax = 0.2
        epsilon = 1e-5

        # calculate the fieldline
        line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem, periodicity=periodic_global)

        if csystem == 'spherical':
            line[:, 0], line[:, 1], line[:, 2] = sphr2cart(line[:, 0], line[:, 1], line[:, 2])

        x.append(line[:, 0])
        y.append(line[:, 1])
        z.append(line[:, 2])
        length = len(line[:, 0])
        s.append(np.zeros(length))
        if csystem == 'spherical' or csystem == 'cylindical' or not periodic_check:
            ptcons.append(np.vstack([np.arange(index, index+length-1), np.arange(index+1, index+length)]).T)
        else:
            dists = np.r_[np.sum(np.diff(line, axis=0)**2, axis=1), 0]
            brks = np.unique(np.r_[-1, np.where(dists > 0.9*periodic_dist)[0], dists.shape[0]-1])
            for ib0, ib1 in zip(brks[:-1], brks[1:]):
                ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
        index += length
    
    # add points to model
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=col, line_width=lw, name='Fieldlines')

        # ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=col, tube_radius=None, line_width=lw)

def add_spines():
    print('Adding spines')

    cols = {-2:(0.5, 0, 0.5), -1:(0, 0, 1), 0:(0, 1, 0), 1:(1, 0, 0), 2:(1, 165/255, 0)}

    spines = rd.spines(filename, null_list=nulllist)

    nulls = nulldata[nulllist-1]

    # very similar to ring algorithm without breaks
    for isign in np.unique(nulls.sign):
        # set up lists like rings
        x, y, z, s, ptcons = ( [] for _ in range(5) )
        index = 0
        for inull in nulls.number[nulls.sign == isign]:
            print('Null {:5d}'.format(inull))
            sys.stdout.write("\033[F")
            for spine in spines[inull-1]:
                if csystem == 'spherical':
                    spine[:, 0], spine[:, 1], spine[:, 2] = sphr2cart(spine[:, 0], spine[:, 1], spine[:, 2])
                elif csystem == 'cylindical':
                    spine[:, 0], spine[:, 1], spine[:, 2] = cyl2cart(spine[:, 0], spine[:, 1], spine[:, 2])
                x.append(spine[:, 0])
                y.append(spine[:, 1])
                z.append(spine[:, 2])
                s.append(np.zeros_like(spine[:, 0]))
                if csystem == 'spherical' or csystem == 'cylindical' or not periodic_check:
                    ptcons.append(np.vstack([np.arange(index, index+spine.shape[0]-1), np.arange(index+1, index+spine.shape[0])]).T)
                else:
                    dists = np.r_[np.sum(np.diff(spine, axis=0)**2, axis=1), 0]
                    brks = np.unique(np.r_[-1, np.where(dists > 0.9*periodic_dist)[0], dists.shape[0]-1])
                    for ib0, ib1 in zip(brks[:-1], brks[1:]):
                        ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
                index += spine.shape[0]
        # and plot...
        if len(x) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
            src.mlab_source.dataset.lines = np.vstack(ptcons)
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=4, name=sign_names[isign]+'Spines')
    
def add_separators(hcs=False, colour=None):
    """
    Add separators to model
    hcs: controls whether null or hcs separators
    colour: override custom colours
    """
    print('Adding separators')

    seps, conn = rd.separators(filename, null_list=nulllist, hcs=hcs)
    ringsmax, nskip, _, _, ringnums = rd.ringinfo(filename)

    # simpler version of spines and rings - no need for positive and negative
    x, y, z, s, ptcons = ( [] for _ in range(5) )
    index = 0

    if hcs:
        linecol = (1.0, 0.6470588235294118, 0.0) # orange
    else:
        linecol = (1, 1, 0) # yellow
    
    if colour is not None:
        linecol = colour

    if hcs:
        to_do = [0]
    else:
        to_do = nulllist-1

    for inull in to_do:
        print('Null {:5d}'.format(inull+1))
        sys.stdout.write("\033[F")
        num = np.count_nonzero(ringnums[inull] > 0)*nskip
        for con, sep in zip(conn[inull], seps[inull]):
            if con in nulllist:
                if csystem == 'spherical':
                    sep[:, 0], sep[:, 1], sep[:, 2] = sphr2cart(sep[:, 0], sep[:, 1], sep[:, 2])
                elif csystem == 'cylindical':
                    sep[:, 0], sep[:, 1], sep[:, 2] = cyl2cart(sep[:, 0], sep[:, 1], sep[:, 2])
                if sep.shape[0] > num: continue
                x.append(sep[:, 0])
                y.append(sep[:, 1])
                z.append(sep[:, 2])
                s.append(np.zeros_like(sep[:, 0]))
                if csystem == 'spherical' or csystem == 'cylindical' or not periodic_check:
                    ptcons.append(np.vstack([np.arange(index, index+sep.shape[0]-1), np.arange(index+1, index+sep.shape[0])]).T)
                else:
                    dists = np.r_[np.sum(np.diff(sep, axis=0)**2, axis=1), [0]]
                    brks = np.unique(np.r_[[-1], np.where(dists > 0.9*periodic_dist)[0], dists.shape[0]-1])
                    for ib0, ib1 in zip(brks[:-1], brks[1:]):
                        ptcons.append(np.vstack([np.arange(index+ib0+1, index+ib1), np.arange(index+ib0+2, index+ib1+1)]).T)
                index += sep.shape[0]
    
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=linecol, line_width=6, name='Separators')

def add_nulls(size=1):
    print("Adding nulls")

    cols = {-2:(0.5, 0, 0.5), -1:(0, 0, 1), 0:(0, 1, 0), 1:(1, 0, 0), 2:(1, 165/255, 0)}

    boxsize = min([xx[-1] - xx[0], yy[-1] - yy[0], zz[-1] - zz[0]])/40

    # find a nice size for the radius of the null points
    r = max([boxsize, ds])*size

    # pick out only the nulls required
    nulldata1 = nulldata[nulllist-1]

        for sign in np.unique(nulldata1.sign):
        pos = nulldata1.pos[nulldata1.sign == sign]
        if csystem == 'spherical':
            pos[:, 0], pos[:, 1], pos[:, 2] = sphr2cart(pos[:, 0], pos[:, 1], pos[:, 2])
        ml.points3d(pos[:, 0], pos[:, 1], pos[:, 2], color=cols[sign], scale_factor=r, resolution=32, name=sign_names[sign]+'Nulls')

    # r = max([boxsize, ds])
    # r = r*size
    # ntheta, nphi = 16, 32
    # theta, phi = np.mgrid[0:np.pi:ntheta*1j, 0:2*np.pi:nphi*1j]
    # x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)

    # for inull in nulllist:
    #     pos = nulldata[inull].pos
    #     if csystem == 'spherical':
    #         pos[0], pos[1], pos[2] = sphr2cart(*pos)
    #     ml.mesh(x + pos[0], y + pos[1], z + pos[2], color=cols[nulldata[inull].sign])

def change_null_size(size):
    # not working properly... also removes outline
    objs = ml.gcf().children
    count_remove = 0
    for obj in objs:
        if 'Nulls' in obj.name:
            count_remove += 1
    while count_remove > 0:
    objs = ml.gcf().children
    for obj in objs:
        if 'Nulls' in obj.name:
            obj.remove()
                count_remove -= 1
    
    add_nulls(size)

def add_box():
    print("Adding box")
    # create line with all corners to be normalised
    line = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 1],
                     [1, 0, 0], [1, 1, 0], [1, 1, 1],
                     [1, 1, 0], [0, 1, 0], [0, 1, 1],
                     [0, 1, 0], [0, 0, 0], [0, 0, 1],
                     [1, 0, 1], [1, 1, 1], [0, 1, 1],
                     [0, 0, 1]], dtype=np.float64)

    # rescale line to form box around the proper coordinate system
    line[:, 0] = line[:, 0]*(xx[-1] - xx[0]) + xx[0]
    line[:, 1] = line[:, 1]*(yy[-1] - yy[0]) + yy[0]
    line[:, 2] = line[:, 2]*(zz[-1] - zz[0]) + zz[0]

    ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=(0,0,0), tube_radius=None, line_width=1, name='Box')

    # ml.outline(extent=[*xx[[0, -1]], *yy[[0, -1]], *zz[[0, -1]]], color=(0,0,0))
    
def add_axes():
    ml.axes(extent=[*xx[[0, -1]], *yy[[0, -1]], *zz[[0, -1]]])

def add_base(coord, pos, vmin=-10, vmax=10):
    """
    Adds a magnetogram on the side of the box
    coord: 'x', 'y' or 'z'
    pos: 'top' or 'bottom'
    """
    if csystem != 'cartesian':
        print('Note: Designed for cartesian coordinates')

    if pos == 'top':
        ipos = -1
    elif pos == 'bottom':
        ipos = 0
    
    if coord == "x":
        y, z = np.meshgrid(yy, zz, indexing='ij')
        x = np.ones_like(y)*xx[ipos]
        base = ml.mesh(x, y, z, scalars=-bgrid[ipos, :, :, 0], colormap='Greys', vmin=vmin, vmax=vmax)
    elif coord == "y":
        x, z = np.meshgrid(xx, zz, indexing='ij')
        y = np.ones_like(x)*yy[ipos]
        base = ml.mesh(x, y, z, scalars=-bgrid[:, ipos, :, 1], colormap='Greys', vmin=vmin, vmax=vmax)
    elif coord == "z":
        x, y = np.meshgrid(xx, yy, indexing='ij')
        z = np.ones_like(y)*zz[ipos]
        base = ml.mesh(x, y, z, scalars=-bgrid[:, :, ipos, 2], colormap='Greys', vmin=vmin, vmax=vmax)

def add_sun():
    # make theta, phi grids of bgrid coords and create the coords of sphere
    theta, phi = np.meshgrid(yy, zz, indexing='ij')
    x = xx[0] * np.sin(theta) * np.cos(phi)
    y = xx[0] * np.sin(theta) * np.sin(phi)
    z = xx[0] * np.cos(theta)

    # -bgrid because otherwise want reversed colourtable
    ml.mesh(x, y, z, scalars=-bgrid[0, :, :, 0],  colormap='Greys', vmin=-10, vmax=10, name='Solar surface')

    # make z-axis
    ml.plot3d([0, 0], [0, 0], [-xx.max(), xx.max()], color=(0, 0, 0), tube_radius=None, line_width=4, name='Z-axis')

def save():
    ml.savefig('figures/' + rd.prefix(filename) + '-model3d.png')

def add_surface():
    print('Adding separatrix surface rings')

    rings, breaks, assocs = rd.rings(filename, breaks=True, assocs=True)

    cols = {-1:(0.5, 0.5, 1), 0:(0.5, 1, 0.5), 1:(1, 0.5, 0.5)}

    acc = 6

    # ringlist = []
    # trianglelist = []
    # ptnums = [0]
    # for lst in break1: ptnums.append(ptnums[-1] + len(lst))
    # for ring in rings[inull]: ringlist.extend(ring.tolist())

    if True:
        for inull in nulllist-1:
            ringlist = []
            trianglelist = []
            ptnums = [0]
            for lst in breaks[inull]: ptnums.append(ptnums[-1] + len(lst))
            ptnums = ptnums[:-1]
            # print(ptnums)
            for ring in rings[inull]: ringlist.extend(ring.tolist())
            ringlist = np.array(ringlist)
            print('Null {}'.format(inull+1))
            rings1 = rings[inull][::-1]
            breaks1 = breaks[inull][::-1]
            assocs1 = assocs[inull][::-1]
            ptnums = ptnums[::-1]
            for iring, ring in enumerate(rings1[:-1]):
                dists = np.r_[np.sum(np.diff(ring, axis=0)**2, axis=1), [0]]
                breaks1[iring][dists > acc] = 1
                brks = np.r_[[-1], np.where(breaks1[iring] == 1)[0], [breaks1[iring].shape[0]-1]] + 1
                for ipt, pt in enumerate(ring):
                    if ipt not in brks-1 and ipt != brks[-1]:
                        if assocs1[iring][ipt] != assocs1[iring][ipt+1]:
                            tri1 = [ptnums[iring] + ipt, ptnums[iring+1] + assocs1[iring][ipt] - 1, ptnums[iring+1] + assocs1[iring][ipt]]
                            trianglelist.append(tri1)
                        tri2 = [ptnums[iring] + ipt, ptnums[iring+1] + assocs1[iring][ipt], ptnums[iring] + ipt + 1]
                        trianglelist.append(tri2)
            # return trianglelist, ringlist, ptnums
            ml.triangular_mesh(ringlist[:,0], ringlist[:,1], ringlist[:,2], trianglelist, color=cols[nulldata[inull].sign], tube_radius=None, opacity=0.5, mask_points=2)
        
    if False:
        nskp = 20 #nskipglob
        for inull in nulllist-1:
            ringlist = []
            trianglelist = []
            ptnums = [0]
            for lst in breaks[inull][::-nskp]: ptnums.append(ptnums[-1] + len(lst))
            ptnums = ptnums[:-1]
            # print(ptnums)
            for ring in rings[inull][::-nskp]: ringlist.extend(ring.tolist())
            ringlist = np.array(ringlist)
            print('Null {}'.format(inull+1))
            rings1 = rings[inull][::-nskp]
            breaks1 = breaks[inull][::-nskp]
            assocs1 = assocs[inull][::-1]
            # ptnums = ptnums[::-1]
            for iring, ring in enumerate(rings1[:-1]):
                if iring/nskp - iring//nskp > 0: 
                    dists = np.r_[np.sum(np.diff(ring, axis=0)**2, axis=1), [0]]
                    breaks1[iring][dists > acc] = 1
                    brks = np.r_[[-1], np.where(breaks1[iring] == 1)[0], [breaks1[iring].shape[0]-1]] + 1
                    newass = assocs1[iring*nskp]
                    for iskp in range(nskp):
                        newass = assocs1[iring*nskp+iskp+1][newass-1]
                    for ipt, pt in enumerate(ring):
                        if ipt not in brks-1 and ipt != brks[-1]:
                            if newass[ipt] != newass[ipt+1]:
                                tri1 = [ptnums[iring] + ipt, ptnums[iring+1] + newass[ipt] - 1, ptnums[iring+1] + newass[ipt]]
                                trianglelist.append(tri1)
                            tri2 = [ptnums[iring] + ipt, ptnums[iring+1] + newass[ipt], ptnums[iring] + ipt + 1]
                            trianglelist.append(tri2)
            # return trianglelist, ringlist, ptnums
            ml.triangular_mesh(ringlist[:,0], ringlist[:,1], ringlist[:,2], trianglelist, color=cols[nulldata[inull].sign], tube_radius=None, opacity=1)

def sphr2cart(rs, ts, ps):
    # convert (r, theta, phi) to (x, y, z)
    return rs*np.sin(ts)*np.cos(ps), rs*np.sin(ts)*np.sin(ps), rs*np.cos(ts)

def cyl2cart(Rs, ps, zs):
    # convert (R, phi, z) to (x, y, z)
    return Rs*np.cos(ps), Rs*np.sin(ps), zs
