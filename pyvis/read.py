import numpy as np
import io
import glob

files = sorted(glob.glob('data/*.dat'))

def prefix(filename):
    return filename[5:-4]

def synmap(filename):

    with io.open(filename, 'rb') as fieldfile:
        nlon = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
        nlat = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
        bsyn = np.fromfile(fieldfile, count=nlon*nlat, dtype=np.float64).reshape(nlat,nlon)
        lons = np.fromfile(fieldfile, count=nlon, dtype=np.float64)
        lats = np.fromfile(fieldfile, count=nlat, dtype=np.float64)

        return bsyn, lons, lats

def field(filename):

    with io.open(filename, 'rb') as fieldfile:
        nrad = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
        ntheta = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
        nphi = np.asscalar(np.fromfile(fieldfile, count=1, dtype=np.int32))
        br = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
        bt = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
        bp = np.fromfile(fieldfile, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
        rads = np.fromfile(fieldfile, count=nrad, dtype=np.float64)
        thetas = np.fromfile(fieldfile, count=ntheta, dtype=np.float64)
        phis = np.fromfile(fieldfile, count=nphi, dtype=np.float64)

    return br, bt, bp, rads, thetas, phis

def nulls(filename, simple=False):

    with io.open('output/'+prefix(filename)+'-nullpos.dat', 'rb') as nullfile:
        nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
        gridpos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
        pos = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)

    if simple == False:
        with io.open('output/'+prefix(filename)+'-nulldata.dat', 'rb') as nullfile:
            nnulls = np.asscalar(np.fromfile(nullfile, dtype=np.int32, count=1))
            signs = np.fromfile(nullfile, dtype=np.int32, count=nnulls)
            spines = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
            fans = np.fromfile(nullfile, dtype=np.float64, count=3*nnulls).reshape(nnulls,3)
            warning = np.fromfile(nullfile, dtype=np.float64, count=nnulls)

        nulls = np.recarray(nnulls, dtype=[('number',np.int32),
            ('pos',np.float64,3),('gridpos',np.float64,3),
            ('sign',np.int32),('spine',np.float64,3),('fan',np.float64,3)])

        nulls.sign = signs
        nulls.spine = spines
        nulls.fan = fans
        nulls.warning = warning
    else:
        nulls = np.recarray(nnulls, dtype=[('number',np.int32),
            ('pos',np.float64,3),('gridpos',np.float64,3)])

    nulls.number = np.arange(nnulls, dtype=np.int32)+1
    nulls.pos = pos
    nulls.gridpos = gridpos

    return nulls

def separators(filename, lines=True, connectivity=True, hcs=False):

    nulldata = nulls(filename, simple=True)

    if hcs == False:
        connectivityfile = 'output/'+prefix(filename)+'-connectivity.dat'
        separatorsfile = 'output/'+prefix(filename)+'-separators.dat'
    else:
        connectivityfile = 'output/'+prefix(filename)+'-hcs-connectivity.dat'
        separatorsfile = 'output/'+prefix(filename)+'-hcs-separators.dat'

    conlist = []
    seplist = []
    inull = 1
    coni = []
    sepi = []

    with io.open(connectivityfile, 'rb') as sepinfo:
        with io.open(separatorsfile, 'rb') as seps:
            start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
            #need to deal with how big the lists are...
            while start >= 0 or inull <= nulldata.number.max():
                if inull == start:
                    end = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
                    sepinfo.seek(8, 1)
                    length = np.asscalar(np.fromfile(seps, dtype=np.int32, count=1))
                    separator = np.fromfile(seps, dtype=np.float64, count=3*length).reshape(-1,3)
                    coni.append(end)
                    sepi.append(separator)
                    start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
                if inull != start:
                    inull += 1
                    conlist.append(coni)
                    seplist.append(sepi)
                    sepi = []
                    coni = []
    # else:
    #   connectivityfile = 'output/'+prefix(filename)+'-connectivity-hcs.dat'
    #   separatorsfile = 'output/'+prefix(filename)+'-separators-hcs.dat'

    #   with io.open(connectivityfile, 'rb') as sepinfo:
    #     with io.open(separatorsfile, 'rb') as seps:
    #       start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
    #       while start == 0:
    #         end = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))
    #         sepinfo.seek(8, 1)
    #         length = np.asscalar(np.fromfile(seps, dtype=np.int32, count=1))
    #         separator = np.fromfile(seps, dtype=np.float64, count=3*length).reshape(-1,3)
    #         conlist += [end]
    #         seplist += [separator]
    #         start = np.asscalar(np.fromfile(sepinfo, dtype=np.int32, count=1))

    if lines == False:
        return conlist
    elif connectivity == False:
        return seplist
    else:
        return seplist, conlist

def spines(filename):

    nulldata = nulls(filename, simple=True)

    spinelist = []

    with io.open('output/'+prefix(filename)+'-spines.dat', 'rb') as spinefile:
        for inull in nulldata.number:
            spinelisti = []
            for ispine in range(2):
                length = np.asscalar(np.fromfile(spinefile, dtype=np.int32, count=1))
                spinelisti.append(np.fromfile(spinefile, dtype=np.float64, count=3*length).reshape(-1,3))
            spinelist.append(spinelisti)

    return spinelist

def rings(filename, allinfo=False, nskip=1):

    nulldata = nulls(filename, simple=True)

    ringlist = []
    breaklist = []
    assoclist = []

    with io.open('output/'+prefix(filename)+'-ringinfo.dat', 'rb') as ringinfo:
        with io.open('output/'+prefix(filename)+'-rings.dat', 'rb') as ringfile:
            ringsmax = np.asscalar(np.fromfile(ringinfo, dtype=np.int32, count=1))
            for inull in nulldata.number:
                assoclisti = []
                breaklisti = []
                ringlisti = []
                lengths = np.fromfile(ringinfo, dtype=np.int32, count=ringsmax)
                iring = 0
                while iring < ringsmax and lengths[iring] > 0:
                    assoclisti.append(np.fromfile(ringfile, dtype=np.int32, count=lengths[iring]))
                    breaklisti.append(np.fromfile(ringfile, dtype=np.int32, count=lengths[iring]))
                    ringlisti.append(np.fromfile(ringfile, dtype=np.float64, count=3*lengths[iring]).reshape(-1,3))
                    iskip = 1
                    while iskip + iring < ringsmax and iskip < nskip:
                        ringfile.seek(lengths[iring+iskip]*32, 1)
                        iskip += 1
                    iring += nskip
                ringlist.append(ringlisti)
                breaklist.append(breaklisti)
                assoclist.append(assoclisti)

    if allinfo == True:
        return ringlist, breaklist, assoclist
    else:
        return ringlist
