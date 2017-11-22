import numpy as np
import glob
import sys

files = sorted(glob.glob('data/*.dat'))

def prefix(filename):
    return filename[5:-4]

def synmap(filename):
    with open(filename, 'rb') as fieldfile:
        nlon, nlat = np.fromfile(fieldfile, count=2, dtype=np.int32)

        bsyn = np.fromfile(fieldfile, count=nlon*nlat, dtype=np.float64).reshape((nlon, nlat), order='f')

        lons = np.fromfile(fieldfile, count=nlon, dtype=np.float64)
        lats = np.fromfile(fieldfile, count=nlat, dtype=np.float64)

    return bsyn, lons, lats

def field(filename):
    with open(filename, 'rb') as fieldfile:
        nrad, ntheta, nphi = np.fromfile(fieldfile, count=3, dtype=np.int32)

        shape = (nrad, ntheta, nphi)
        num = nrad*ntheta*nphi
        br = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')
        bt = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')
        bp = np.fromfile(fieldfile, count=num, dtype=np.float64).reshape(shape, order='f')

        rads = np.fromfile(fieldfile, count=nrad, dtype=np.float64)
        thetas = np.fromfile(fieldfile, count=ntheta, dtype=np.float64)
        phis = np.fromfile(fieldfile, count=nphi, dtype=np.float64)

    return br, bt, bp, rads, thetas, phis

def nulls(filename, simple=False):
    with open('output/'+prefix(filename)+'-nullpos.dat', 'rb') as nullfile:
        # get three data sets from the null finder file
        nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
        gridpos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)
        pos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls, 3)

    if simple == False:
        # if you want the sign finder data too...
        with open('output/'+prefix(filename)+'-nulldata.dat', 'rb') as nullfile:
            nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
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
        connectivityfile = 'output/'+prefix(filename)+'-connectivity.dat'
        separatorsfile = 'output/'+prefix(filename)+'-separators.dat'
        allnulls_list = nulldata.number
    else:
        # filenames for the hcs
        connectivityfile = 'output/'+prefix(filename)+'-hcs-connectivity.dat'
        separatorsfile = 'output/'+prefix(filename)+'-hcs-separators.dat'
        allnulls_list = [1, 2] # need to fix this
    
    # if none set, read in data for all nulls
    if null_list is None:
        null_list = allnulls_list

    # create final lists
    conlist = []
    seplist = []
    
    with open(connectivityfile, 'rb') as sepinfo:
        with open(separatorsfile, 'rb') as seps:
            start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
            for inull in allnulls_list:
                coni = []
                sepi = []
                while start == inull:
                    # read length of the separator
                    length = np.asscalar(np.fromfile(seps, dtype=np.int32, count=1))
                    if inull in null_list:
                        # read the connectivity data and separator points
                        coni.append(np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1)))
                        sepinfo.seek(8, 1)
                        sepi.append(np.fromfile(seps, dtype=np.float64, count=3*length).reshape(-1, 3))
                    else:
                        # skip this separator since not connected to this null
                        sepinfo.seek(12, 1)
                        seps.seek(3*length*8, 1)
                    start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
                # add data to final list, empty if null not in null_list
                conlist.append(coni)
                seplist.append(sepi)
    
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

    with open('output/'+prefix(filename)+'-spines.dat', 'rb') as spinefile:
        for inull in nulldata.number:
            spinelisti = []
            for ispine in range(2): # spine in each direction
                length = np.asscalar(np.fromfile(spinefile, dtype=np.int32, count=1))
                if inull in null_list:
                    spinelisti.append(np.fromfile(spinefile, dtype=np.float64, count=3*length).reshape(-1, 3))
                else:
                    # skip spine if not required
                    spinefile.seek(length*3*8, 1)
            spinelist.append(spinelisti)

    return spinelist

def rings(filename, allinfo=False, nskip=1, null_list=None):

    nulldata = nulls(filename, simple=True)

    ringlist = []
    breaklist = []
    assoclist = []

    if null_list is None:
        null_list = nulldata.number

    with open('output/'+prefix(filename)+'-ringinfo.dat', 'rb') as ringinfo:
        with open('output/'+prefix(filename)+'-rings.dat', 'rb') as ringfile:
            ringsmax = np.asscalar(np.fromfile(ringinfo, dtype=np.int32, count=1))
            for inull in nulldata.number:
                print('Reading rings from null {:5d}'.format(inull))
                sys.stdout.write("\033[F")
                assoclisti = []
                breaklisti = []
                ringlisti = []
                lengths = np.fromfile(ringinfo, dtype=np.int32, count=ringsmax)
                if inull in null_list:
                    for iring, length in enumerate(lengths[lengths != 0]):
                        if iring % nskip == 0:
                            assoclisti.append(np.fromfile(ringfile, dtype=np.int32, count=length))
                            breaklisti.append(np.fromfile(ringfile, dtype=np.int32, count=length))
                            ringlisti.append(np.fromfile(ringfile, dtype=np.float64, count=3*length).reshape(-1, 3))
                        else:
                            # skip ring if not a multiple of nskip
                            ringfile.seek(length*32, 1)
                else:
                    # skip all rings from null if not required
                    ringfile.seek(lengths.sum(dtype=np.int64)*32, 1)
                ringlist.append(ringlisti)
                breaklist.append(breaklisti)
                assoclist.append(assoclisti)

    if allinfo == True:
        return ringlist, breaklist, assoclist
    else:
        return ringlist

                    # iring = 0
                    # while iring < ringsmax and lengths[iring] > 0:
                    #     assoclisti.append(np.fromfile(ringfile, dtype=np.int32, count=lengths[iring]))
                    #     breaklisti.append(np.fromfile(ringfile, dtype=np.int32, count=lengths[iring]))
                    #     ringlisti.append(np.fromfile(ringfile, dtype=np.float64, count=3*lengths[iring]).reshape(-1, 3))
                    #     iskip = 1
                    #     while iskip + iring < ringsmax and iskip < nskip:
                    #         ringfile.seek(lengths[iring+iskip]*32, 1)
                    #         iskip += 1
                    #     iring += nskip
