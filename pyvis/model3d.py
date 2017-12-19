from __future__ import print_function, division
import numpy as np
import mayavi.mlab as ml
from . import read as rd
from . import fieldline3d as fl
import sys
import vtk

# turn of warnings while vtk/mayavi compatibility is fixed
vtk.vtkObject.GlobalWarningDisplayOff()

def make(fname, addlist, null_list=None, box=True, fieldlines=None, linecolor=(0,0,0), nskip=20,
    nullrad=1, nfanlines=40, nring=None, colquant=None, coordsystem='cartesian'):
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
        nskip: how many rings to skip in plotting (also skips points in spines and separators for speed)"""

    global bgrid, xx, yy, zz, nulldata, ds, filename, nskipglob, nulllist, csystem

    csystem = coordsystem
    print('Using {} coordinates'.format(coordsystem))

    ml.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 800))

    nskipglob = nskip

    filename = fname

    field = rd.field(filename)
    bgrid = np.zeros(field[0].shape + (3,), dtype=np.float64)
    for i in range(3):
        bgrid[:, :, :, i] = field[i]

    xx = field[3]
    yy = field[4]
    zz = field[5]

    ds = min([np.diff(xx).min(), np.diff(yy).min(), np.diff(zz).min()])/2

    nulldata = rd.nulls(filename)

    set_null_list(null_list)

    if coordsystem == 'spherical':
        add_sun()
        box = False
    if box: add_box()

    add_structures(*addlist, nullrad=nullrad, nfanlines=nfanlines, nring=nring)

    if fieldlines is not None: add_fieldlines(fieldlines, col=linecolor, colquant=colquant)

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

    cols = {-1:(0.5, 0.5, 1), 0:(0.5, 1, 0.5), 1:(1, 0.5, 0.5)}

    # new very efficient routine for plotting many lines
    # two lists, one for positive and the other for negative nulls
    x, y, z, s, ptcons = ( [[],[]] for _ in range(5) )
    index = [0, 0]

    for inull in nulllist-1:
        print('Null {:5d}'.format(inull+1))
        sys.stdout.write("\033[F")
        if nulldata[inull].sign == 1:
            il = 0
        else:
            il = 1
        for iring, ring in enumerate(rings[inull]):
            # convert points if required
            if csystem == 'spherical':
                ring[:, 0], ring[:, 1], ring[:, 2] = sphr2cart(ring[:, 0], ring[:, 1], ring[:, 2])
            # add ring points to lists
            x[il].append(ring[:,0])
            y[il].append(ring[:,1])
            z[il].append(ring[:,2])
            s[il].append(np.zeros_like(ring[:,0]))
            # use break data to plot the individual lines in each ring as the break apart
            brks = np.r_[[-1], np.where(breaks[inull][iring] == 1)[0], [breaks[inull][iring].shape[0]-1]] + 1
            for ib0, ib1 in zip(brks[:-1], brks[1:]):
                if ib0 != ib1:
                    # add the right indicies based on the breaks
                    ptcons[il].append(np.vstack([np.arange(index[il]+ib0, index[il]+ib1-1), np.arange(index[il]+ib0+1, index[il]+ib1)]).T)
            index[il] += ring.shape[0]
    
    # add points to model
    for il, isign in zip(range(2), [1, -1]):
        if len(x[il]) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x[il]), np.hstack(y[il]), np.hstack(z[il]), np.hstack(s[il]))
            src.mlab_source.dataset.lines = np.vstack(ptcons[il])
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=1)

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
            brks = np.r_[[-1], np.where(breaks[inull][iring] == 1)[0], [breaks[inull][iring].shape[0]-1]] + 1
            for ib0, ib1 in zip(brks[:-1], brks[1:]):
                if ib0 != ib1:
                    # add the right indicies based on the breaks
                    ptcons.append(np.vstack([np.arange(index+ib0, index+ib1-1), np.arange(index+ib0+1, index+ib1)]).T)
            index += ring.shape[0]
    
    # add points to model
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=(0, 1, 0), line_width=1)
    
    for inull in range(0, len(rings), 2):
        ml.plot3d(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2], color=(0, 1, 0), line_width=6, tube_radius=None)

def add_sepsurf_flines(nlines, nring=None):
    print('Adding separatrix surface field lines')

    rings = rd.rings(filename, nskip=nskipglob, null_list=nulllist)

    cols = {-1:(0.5, 0.5, 1), 0:(0.5, 1, 0.5), 1:(1, 0.5, 0.5)}

    x, y, z, s, ptcons = ( [[], []] for _ in range(5) )
    index = [0, 0]

    for inull in nulllist-1:
        print('Null {}'.format(inull+1))
        sys.stdout.write("\033[F")

        if nulldata[inull].sign == 1:
            il = 0
        else:
            il = 1

        # select the right ring to trace from
        if nring is None:
            ring = rings[inull][len(rings[inull])//5]
        else:
            ring = rings[inull][nring]
        
        nskip = len(ring[:, 0])//nlines
        
        for startpt in ring[::nskip, :]:
            # choose some good parameters
            h = 1e-2
            hmin = 1e-3
            hmax = 0.5
            epsilon = 1e-5

            # calculate the fieldline
            line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem)
            
            # cut off the fieldline at the point closest to the null - only want the fan, not the spine
            dists = np.sqrt((line[:, 0] - nulldata[inull].pos[0])**2 +
                (line[:, 1] - nulldata[inull].pos[1])**2 + (line[:, 2] - nulldata[inull].pos[2])**2)
            imin = dists.argmin()
            line = line[0:imin+1, :] if nulldata[inull].sign == -1 else line[imin:, :]

            if csystem == 'spherical':
                line[:, 0], line[:, 1], line[:, 2] = sphr2cart(line[:, 0], line[:, 1], line[:, 2])
            
            x[il].append(line[:, 0])
            y[il].append(line[:, 1])
            z[il].append(line[:, 2])
            length = len(line[:, 0])
            s[il].append(np.zeros(length))
            ptcons[il].append(np.vstack([np.arange(index[il], index[il]+length-1), np.arange(index[il]+1, index[il]+length)]).T)
            index[il] += length
    
    # add points to model
    for il, isign in zip(range(2), [1, -1]):
        if len(x[il]) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x[il]), np.hstack(y[il]), np.hstack(z[il]), np.hstack(s[il]))
            src.mlab_source.dataset.lines = np.vstack(ptcons[il])
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=1)

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
        ml.pipeline.surface(lines, color=(0, 1, 0), line_width=1)
                
    for inull in range(0, len(rings), 2):
        rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2] = sphr2cart(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2])
        ml.plot3d(rings[inull][0][:, 0], rings[inull][0][:, 1], rings[inull][0][:, 2], color=(0, 1, 0), line_width=6, tube_radius=None)

def add_fieldlines(startpts, col=(0, 0, 0), colquant=None):
    print('Adding separatrix surface field lines')

    for i, startpt in enumerate(startpts, start=1):
        print('Calculating fieldline {}'.format(i))
        sys.stdout.write("\033[F")

        # choose fieldline parameters
        h = 1e-2
        hmin = 1e-3
        hmax = 0.2
        epsilon = 1e-5

        # calculate the fieldline
        line = fl.fieldline3d(startpt, bgrid, xx, yy, zz, h, hmin, hmax, epsilon, coordsystem=csystem)

        if csystem == 'spherical':
            line[:, 0], line[:, 1], line[:, 2] = sphr2cart(line[:, 0], line[:, 1], line[:, 2])

        if colquant is None:
            ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=col, tube_radius=None)
        else:
            vals = np.zeros(line.shape[0], dtype=np.float64)

            for iline, pt in enumerate(line):
                vals[iline] = fl.trilinearscalar3d(pt, colquant, xx, yy, zz)

            ml.plot3d(line[:, 0], line[:, 1], line[:, 2], vals, tube_radius=None, colormap='jet')

def add_spines():
    print('Adding spines')

    cols = {-1:(0, 0, 1), 0:(0, 1, 0), 1:(1, 0, 0)}

    spines = rd.spines(filename, null_list=nulllist)

    # set up lists like rings
    x, y, z, s, ptcons = ( [[],[]] for _ in range(5) )
    index = [0, 0]

    # very similar to ring algorithm without breaks
    for inull in nulllist-1:
        print('Null {:5d}'.format(inull+1))
        sys.stdout.write("\033[F")
        if nulldata[inull].sign == 1:
            il = 0
        else:
            il = 1
        for spine in spines[inull]:
            if csystem == 'spherical':
                spine[:, 0], spine[:, 1], spine[:, 2] = sphr2cart(spine[:, 0], spine[:, 1], spine[:, 2])
            x[il].append(spine[:,0])
            y[il].append(spine[:,1])
            z[il].append(spine[:,2])
            s[il].append(np.zeros_like(spine[:,0]))
            ptcons[il].append(np.vstack([np.arange(index[il], index[il]+spine.shape[0]-1), np.arange(index[il]+1, index[il]+spine.shape[0])]).T)
            index[il] += spine.shape[0]
    
    for il, isign in zip(range(2), [1, -1]):
        if len(x[il]) > 0:
            src = ml.pipeline.scalar_scatter(np.hstack(x[il]), np.hstack(y[il]), np.hstack(z[il]), np.hstack(s[il]))
            src.mlab_source.dataset.lines = np.vstack(ptcons[il])
            src.update()
            
            lines = ml.pipeline.stripper(src)
            ml.pipeline.surface(lines, color=cols[isign], line_width=4)

def add_separators(hcs=False):
    print('Adding separators')

    seps, conn = rd.separators(filename, null_list=nulllist, hcs=hcs)

    # simpler version of spines and rings - no need for positive and negative
    x, y, z, s, ptcons = ( [] for _ in range(5) )
    index = 0

    if hcs:
        to_do = [0]
    else:
        to_do = nulllist-1

    for inull in to_do:
        print('Null {:5d}'.format(inull+1))
        sys.stdout.write("\033[F")
        for con, sep in zip(conn[inull], seps[inull]):
            if con in nulllist:
                if csystem == 'spherical':
                    sep[:, 0], sep[:, 1], sep[:, 2] = sphr2cart(sep[:, 0], sep[:, 1], sep[:, 2])
                x.append(sep[:,0])
                y.append(sep[:,1])
                z.append(sep[:,2])
                s.append(np.zeros_like(sep[:,0]))
                ptcons.append(np.vstack([np.arange(index, index+sep.shape[0]-1), np.arange(index+1, index+sep.shape[0])]).T)
                index += sep.shape[0]
    
    if len(x) > 0:
        src = ml.pipeline.scalar_scatter(np.hstack(x), np.hstack(y), np.hstack(z), np.hstack(s))
        src.mlab_source.dataset.lines = np.vstack(ptcons)
        src.update()
        
        lines = ml.pipeline.stripper(src)
        ml.pipeline.surface(lines, color=(0, 0.5, 0), line_width=6)

def add_nulls(size=1):
    print("Adding nulls")

    cols = {-1:(0, 0, 1), 0:(0, 1, 0), 1:(1, 0, 0)}

    boxsize = min([xx[-1] - xx[0], yy[-1] - yy[0], zz[-1] - zz[0]])/40

    r = max([boxsize, ds])*size

    nulldata1 = nulldata[nulllist-1]

    for sign in [-1, 0, 1]:
        pos = nulldata1.pos[nulldata1.sign == sign]
        if csystem == 'spherical':
            pos[:, 0], pos[:, 1], pos[:, 2] = sphr2cart(pos[:, 0], pos[:, 1], pos[:, 2])
        ml.points3d(pos[:, 0], pos[:, 1], pos[:, 2], color=cols[sign], scale_factor=r, resolution=32)

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

def add_box():
    print("Adding box")
    box = np.zeros((3, 2), dtype=np.float64)
    box[:, 0] = np.array([xx.min(), yy.min(), zz.min()])
    box[:, 1] = np.array([xx.max(), yy.max(), zz.max()])
    line = np.array([[0,0,0],[1,0,0],[1,0,1],[1,0,0],[1,1,0],[1,1,1],[1,1,0],
        [0,1,0],[0,1,1],[0,1,0],[0,0,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1],[0,0,1]], dtype=np.float64)

    line[:,0] = line[:,0]*(box[0,1] - box[0,0]) + box[0,0]
    line[:,1] = line[:,1]*(box[1,1] - box[1,0]) + box[1,0]
    line[:,2] = line[:,2]*(box[2,1] - box[2,0]) + box[2,0]

    # dist = min([n_elements(xx), n_elements(yy), n_elements(zz)])*ds/5
    # oModel -> add, obj_new('idlgrtext', 'x', locations=[box[0,0]+dist,box[1,0],box[2,0]], /onglass)
    # oModel -> add, obj_new('idlgrtext', 'y', locations=[box[0,0],box[1,0]+dist,box[2,0]], /onglass)
    # oModel -> add, obj_new('idlgrtext', 'z', locations=[box[0,0],box[1,0],box[2,0]+dist], /onglass)

    ml.plot3d(line[:, 0], line[:, 1], line[:, 2], color=(0,0,0), tube_radius=None, line_width=1)

def add_base():
    y, z = np.mgrid[yy[0]:yy[-1]:yy.shape[0]*1j, zz[0]:zz[-1]:zz.shape[0]*1j]
    x = np.ones_like(y)*xx[0]
    base = ml.mesh(x, y, z, scalars=-bgrid[0,:,:,0], colormap='Greys', vmin=-10, vmax=10)

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

    if False:
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
        
    if True:
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


            # for ib0, ib1 in zip(brks[:-1], brks[1:]):
            #     if ib0 != ib1:


            # for ib0, ib1 in zip(brks[:-1], brks[1:]):
            #     if ib0 != ib1:

def add_sun():
    ntheta = bgrid.shape[1]*1j
    nphi = bgrid.shape[2]*1j
    theta, phi = np.mgrid[0:np.pi:ntheta, 0:2*np.pi:nphi]
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    ml.mesh(x, y, z, scalars=-bgrid[0,:,:,0],  colormap='Greys', vmin=-10, vmax=10)

    t = np.linspace(-xx.max(),xx.max(),101)
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    z = t.copy()

    ml.plot3d(x, y, z, color=(0,0,0), tube_radius=None, line_width=4)

def sphr2cart(rs, ts, ps):
    return rs*np.sin(ts)*np.cos(ps), rs*np.sin(ts)*np.sin(ps), rs*np.cos(ts)
