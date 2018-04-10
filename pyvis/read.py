import numpy as np
import glob
import sys
import os.path

files = sorted(glob.glob('data/*.dat'))

outprefix = 'output'

def prefix(filename):
    dot = filename[::-1].find('.')
    slash = filename[::-1].find('/')
    return filename[-slash:-dot-1]

def synmap(filename):
    with open(filename, 'rb') as fieldfile:
        nlon, nlat = np.fromfile(fieldfile, count=2, dtype=np.int32)

        bsyn = np.fromfile(fieldfile, count=nlon*nlat, dtype=np.float64).reshape((nlon, nlat), order='f')

        lons = np.fromfile(fieldfile, count=nlon, dtype=np.float64)
        lats = np.fromfile(fieldfile, count=nlat, dtype=np.float64)

    return bsyn, lons, lats

def field(filename):
    with open(filename, 'rb') as fieldfile:
        nx, ny, nz = np.fromfile(fieldfile, count=3, dtype=np.int32)

        shape = (nx, ny, nz)
        num = np.asscalar(nx)*np.asscalar(ny)*np.asscalar(nz)
        bx = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')
        by = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')
        bz = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')

        x = np.fromfile(fieldfile, count=nx, dtype=np.float64)
        y = np.fromfile(fieldfile, count=ny, dtype=np.float64)
        z = np.fromfile(fieldfile, count=nz, dtype=np.float64)

    return bx, by, bz, x, y, z

def nulls(filename, simple=False):
    with open(outprefix+'/'+prefix(filename)+'-nullpos.dat', 'rb') as nullfile:
        # get three data sets from the null finder file
        nnulls, = np.fromfile(nullfile, dtype=np.int32, count=1)
        gridpos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)
        pos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)

    if simple == False:
        # if you want the sign finder data too...
        with open(outprefix+'/'+prefix(filename)+'-nulldata.dat', 'rb') as nullfile:
            nnulls, = np.fromfile(nullfile, dtype=np.int32, count=1)
            signs = np.fromfile(nullfile, dtype=np.int32, count=nnulls)
            spines = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)
            fans = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)
            warning = np.fromfile(nullfile, dtype=np.int32, count=nnulls)

        # create the right record array
        nulls = np.recarray(nnulls, dtype=[('number',np.int32),
            ('pos',np.float64,3), ('gridpos',np.float64,3),
            ('sign',np.int32), ('spine',np.float64,3), ('fan',np.float64,3), ('warning',np.int32)])

        # fill the array
        nulls.sign = signs
        nulls.spine = spines
        nulls.fan = fans
        nulls.warning = warning
    else:
        # just the null finder data
        nulls = np.recarray(nnulls, dtype=[('number',np.int32),
            ('pos',np.float64,3), ('gridpos',np.float64,3)])

    nulls.number = np.arange(nnulls, dtype=np.int32)+1
    nulls.pos = pos
    nulls.gridpos = gridpos

    return nulls

def separators(filename, null_list=None, lines=True, connectivity=True, hcs=False):
    # read in null data
    nulldata = nulls(filename, simple=True)

    if hcs == False:
        # filenames for the nulls
        connectivityfile = outprefix+'/'+prefix(filename)+'-connectivity.dat'
        separatorsfile = outprefix+'/'+prefix(filename)+'-separators.dat'
        allnulls_list = nulldata.number
    else:
        # filenames for the hcs
        connectivityfile = outprefix+'/'+prefix(filename)+'-hcs-connectivity.dat'
        separatorsfile = outprefix+'/'+prefix(filename)+'-hcs-separators.dat'
        allnulls_list = [1] # need to fix this
    
    # if none set, read in data for all nulls
    if null_list is None:
        null_list = allnulls_list

    # create final lists
    conlist = []
    seplist = []
    
    with open(connectivityfile, 'rb') as sepinfo:
        with open(separatorsfile, 'rb') as seps:
            for inull in allnulls_list:
                coni = []
                sepi = []
                nseps, = np.fromfile(sepinfo, dtype=np.int32, count=1)
                for isep in range(nseps):
                    sepinfo.seek(4, 1)
                    coni.append(np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1)))
                    sepinfo.seek(8, 1)
                    length, = np.fromfile(seps, dtype=np.int32, count=1)
                    sepi.append(np.fromfile(seps, dtype=np.float64, count=3*length).reshape(-1, 3))
                if inull in null_list:
                    # add data to final list, empty if null not in null_list
                    conlist.append(coni)
                    seplist.append(sepi)
                else:
                    conlist.append([])
                    seplist.append([])
    
    # return data based on keyword arguments
    if lines == False:
        return conlist
    elif connectivity == False:
        return seplist
    else:
        return seplist, conlist
                
def spines(filename, null_list=None):
    # read in null data
    nulldata = nulls(filename, simple=True)

    # do all nulls if null_list set to none
    if null_list is None:
        null_list = nulldata.number

    spinelist = []

    with open(outprefix+'/'+prefix(filename)+'-spines.dat', 'rb') as spinefile:
        for inull in nulldata.number:
            spinelisti = []
            for ispine in range(2): # spine in each direction
                length, = np.fromfile(spinefile, dtype=np.int32, count=1)
                if inull in null_list:
                    spinelisti.append(np.fromfile(spinefile, dtype=np.float64, count=3*length).reshape(-1, 3))
                else:
                    # skip spine if not required
                    spinefile.seek(length*3*8, 1)
            spinelist.append(spinelisti)

    return spinelist

def ringinfo(filename):
    with open(outprefix+'/'+prefix(filename)+'-ringinfo.dat', 'rb') as ringinfo:
        ringsmax, writeskip, bytesize = np.fromfile(ringinfo, dtype=np.int32, count=3)
        stepsize, = np.fromfile(ringinfo, dtype=np.float64, count=1)
    return ringsmax, writeskip, bytesize, stepsize

def rings(filename, breaks=False, assocs=False, nskip=1, null_list=None, hcs=False):

    nulldata = nulls(filename, simple=True)

    ringlist = []
    breaklist = []
    if assocs: assoclist = []

    if hcs:
        assocs_filename = outprefix+'/'+prefix(filename)+'-hcs-assocs.dat'
        info_filename = outprefix+'/'+prefix(filename)+'-hcs-ringinfo.dat'
        ring_filename = outprefix+'/'+prefix(filename)+'-hcs-rings.dat'
        break_filename = outprefix+'/'+prefix(filename)+'-hcs-breaks.dat'
    else:
        assocs_filename = outprefix+'/'+prefix(filename)+'-assocs.dat'
        info_filename = outprefix+'/'+prefix(filename)+'-ringinfo.dat'
        ring_filename = outprefix+'/'+prefix(filename)+'-rings.dat'
        break_filename = outprefix+'/'+prefix(filename)+'-breaks.dat'
    
    assoc_exist = os.path.isfile(assocs_filename)
    if assocs and not assoc_exist: print('Not reading in associations: file does not exist')
    do_assocs = assocs and assoc_exist

    ringinfo = open(info_filename, 'rb')
    ringfile = open(ring_filename, 'rb')
    brkfile = open(break_filename, 'rb')

    ringsmax, writeskip, bytesize = np.fromfile(ringinfo, dtype=np.int32, count=3)
    if bytesize == 4:
        floattype = np.float32
    elif bytesize == 8:
        floattype = np.float64
    stepsize, = np.fromfile(ringinfo, dtype=np.float64, count=1)
    if hcs: n_hcs, = np.fromfile(ringinfo, dtype=np.int32, count=1)
    ringsmax = np.ceil(ringsmax/writeskip).astype(np.int)

    if hcs:
        to_do = range(1, n_hcs//2+1)
        ndir = 2
        name = 'hcs'
    else:
        to_do = nulldata.number
        ndir = 1
        name = 'null'

    if null_list is None and not hcs:
        null_list = nulldata.number
    elif hcs:
        null_list = to_do

    for inull in to_do:
        print('Reading rings from {} {:5d}'.format(name, inull))
        sys.stdout.write("\033[F")
        for dir in range(ndir):
            breaklisti = []
            ringlisti = []
            lengths = np.fromfile(ringinfo, dtype=np.int32, count=ringsmax)
            if inull in null_list:
                for iring, length in enumerate(lengths[lengths != 0]):
                    if iring % nskip == 0:
                        breaklisti.append(np.fromfile(brkfile, dtype=np.int32, count=length))
                        ringlisti.append(np.fromfile(ringfile, dtype=floattype, count=3*length).reshape(-1, 3))
                    else:
                        # skip ring if not a multiple of nskip
                        ringfile.seek(length*3*bytesize, 1)
                        brkfile.seek(length*4, 1)
            else:
                # skip all rings from null if not required
                total = lengths.sum(dtype=np.int64)
                ringfile.seek(total*3*bytesize, 1)
                brkfile.seek(total*4, 1)
            ringlist.append(ringlisti)
            breaklist.append(breaklisti)
    
    ringinfo.close()
    ringfile.close()
    brkfile.close()

    if do_assocs:
        assfile = open(assocs_filename, 'rb')
        ringinfo = open(info_filename, 'rb')
        np.fromfile(ringinfo, dtype=np.int32, count=3)
        np.fromfile(ringinfo, dtype=np.float64, count=1)
        for inull in nulldata.number:
            assoclisti = []
            lengths = np.fromfile(ringinfo, dtype=np.int32, count=ringsmax)
            if inull in null_list:
                for iring, length in enumerate(lengths[lengths != 0]):
                    if iring % nskip == 0:
                        assoclisti.append(np.fromfile(assfile, dtype=np.int32, count=length))
                    else:
                        assfile.seek(length*4, 1)
            else:
                assfile.seek(lengths.sum(dtype=np.int64)*4, 1)
            assoclist.append(assoclisti)
        assfile.close()
        ringinfo.close()
    
    if breaks == True and do_assocs == False:
        return ringlist, breaklist
    if breaks == False and do_assocs == True:
        return ringlist, assoclist
    if breaks == True and do_assocs == True:
        return ringlist, breaklist, assoclist
    else:
        return ringlist

def cut_sepsurf(r, filename):
    null_nums = []
    sepsurfs = []
    with open(outprefix+'/'+prefix(filename)+'-rings-cut_'+'{:6.4f}.dat'.format(r), 'rb') as ringfile:
        nlines, = np.fromfile(ringfile, dtype=np.int32, count=1)
        for _ in range(nlines):
            inull, = np.fromfile(ringfile, dtype=np.int32, count=1)
            null_nums.append(inull)
            length, = np.fromfile(ringfile, dtype=np.int32, count=1)
            sepsurf = np.fromfile(ringfile, dtype=np.float64, count=3*length).reshape(-1,3)
            sepsurfs.append(sepsurf)
    return null_nums, sepsurfs

def cut_separators(r, filename, hcs=False):
    null_nums = []
    sep_pts = []
    if hcs:
        with open(outprefix+'/'+prefix(filename)+'-hcs-separators-cut_'+'{:6.4f}.dat'.format(r), 'rb') as sepfile:
            npts, = np.fromfile(sepfile, dtype=np.int32, count=1)
            for _ in range(npts):
                end, = np.fromfile(sepfile, dtype=np.int32, count=1)
                null_nums.append(end)
                sep = np.fromfile(sepfile, dtype=np.float64, count=3)
                sep_pts.append(sep)
    else:
        with open(outprefix+'/'+prefix(filename)+'-separators-cut_'+'{:6.4f}.dat'.format(r), 'rb') as sepfile:
            npts, = np.fromfile(sepfile, dtype=np.int32, count=1)
            for _ in range(npts):
                start, end = np.fromfile(sepfile, dtype=np.int32, count=2)
                null_nums.append(end)
                sep = np.fromfile(sepfile, dtype=np.float64, count=3)
                sep_pts.append(sep)
    return null_nums, sep_pts

def cut_spines(r, filename):
    null_nums = []
    spine_pts = []
    with open(outprefix+'/'+prefix(filename)+'-spines-cut_'+'{:6.4f}.dat'.format(r), 'rb') as spinefile:
        npts, = np.fromfile(spinefile, dtype=np.int32, count=1)
        for _ in range(npts):
            inull, = np.fromfile(spinefile, dtype=np.int32, count=1)
            null_nums.append(inull)
            spine = np.fromfile(spinefile, dtype=np.float64, count=3)
            spine_pts.append(spine)
    return null_nums, spine_pts

def cut_hcs(r, filename):
    lines = []
    with open(outprefix+'/'+prefix(filename)+'-hcs-cut_{:6.4f}.dat'.format(r), 'rb') as hcsfile:
        nlines, = np.fromfile(hcsfile, dtype=np.int32, count=1)
        for _ in range(0, nlines):
            ihcs, = np.fromfile(hcsfile, dtype=np.int32, count=1)
            length, = np.fromfile(hcsfile, dtype=np.int32, count=1)
            lines.append(np.fromfile(hcsfile, dtype=np.float64, count=3*length).reshape(-1,3))
    
    return lines