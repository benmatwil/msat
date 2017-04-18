import numpy as np
import matplotlib.pyplot as plt
from . import read as rd
import os, io
from random import shuffle

plt.ion()

# cols = {-1:'blue', 0:'green', 1:'red'}

def get_colours():
    negcolors, poscolors = plt.get_cmap('Blues'), plt.get_cmap('Reds')
    negcolors = [negcolors(i) for i in np.linspace(0.5, 1, np.count_nonzero(nulldata.sign == -1))]
    poscolors = [poscolors(i) for i in np.linspace(0.5, 1, np.count_nonzero(nulldata.sign == 1))]
    shuffle(poscolors)
    shuffle(negcolors)
    colours, ipos, ineg = [], 0, 0
    for inull in range(nulldata.number[-1]):
        if nulldata[inull].sign == -1:
            col = negcolors[ineg]
            ineg += 1
        elif nulldata[inull].sign == 1:
            col = poscolors[ipos]
            ipos += 1
        colours.append(col)
    
    return colours

#################################################################

def start(r, filename, title=False, levels=np.linspace(-20,20,101), colmap=plt.cm.Greys_r):
    global datafile, prefile, nulldata, rstr, colours

    datafile = filename
    prefile = rd.prefix(filename)

    nulldata = rd.nulls(datafile)

    rstr = '{:6.4f}'.format(r)

    def ticks(nsplit, npi):
        from fractions import Fraction as fr
        list = []
        for i in range(nsplit+1):
            string = r'$'
            fract = fr(npi*i,nsplit)
            if fract.numerator == 0:
                string = r'$0$'
            elif fract.denominator != 1:
            string = string + r'\frac{'
            if fract.numerator != 1:
                string = string + str(fract.numerator) 
            string = string + r'\pi$' + r'}{' + str(fract.denominator) + r'}$'
            else:
            if fract.numerator != 1:
                string = string + str(fract.numerator)
            string = string + r'\pi$'
            list.append(string)
        return list

    # os.system(f'./make_cut -i {filename} -r {r}')

    plt.figure(figsize=(16/2.55, 8/2.55)) # set-up for A4
    ax = plt.gca()
    plt.xlabel('Longitude')
    plt.xlim([0, 2*np.pi])
    plt.xticks([0, np.pi/2, np.pi, np.pi*3/2, np.pi*2], [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
    plt.ylabel('Latitude')
    plt.ylim([np.pi, 0])
    plt.yticks([0, np.pi/4, np.pi/2, np.pi*3/4, np.pi], [r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])

    ax.set_rasterization_zorder(0)
    br, _, _, rads, thetas, phis = rd.field(datafile)
    ir = np.where(r >= rads)[0].max()
    plotfield = br[ir, :, :] + (r - rads[ir])/(rads[ir+1] - rads[ir])*(br[ir+1, :, :] - br[ir, :, :])
    plt.contourf(phis, thetas, plotfield, levels, cmap=colmap, extend='both', zorder=-10)
    cb = plt.colorbar(fraction=0.05, pad=0.025)
    diff = levels[-1] - levels[0]
    cb.set_ticks(np.array([0.0, 0.25, 0.5, 0.75, 1.0])*diff + levels[0])
    # print(np.array([0.0, 0.25, 0.5, 0.75, 1.0])*diff + levels[0])
    cb.set_label('Magnetic Field Strength (G)')
    plt.tight_layout()
    if title == True: plt.title('Cut: ' + datafile + ' r = ' + rstr)

    colours = get_colours()

#################################################################

def all(labels=False, ms=3):
    lw = 1
    rings(labels=labels, lw=lw)
    hcs(lw=lw)
    spines(labels=labels, ms=ms)
    separators(labels=labels, ms=ms)

#################################################################

def spines(labels=False, ms=3):
    with io.open('output/'+prefile+'-spines-cut_'+rstr+'.dat', 'rb') as spinefile:
        inull = np.asscalar(np.fromfile(spinefile, dtype=np.int32, count=1))
        while inull > 0:
            spine = np.fromfile(spinefile, dtype=np.float64, count=3)
            plt.plot(spine[2], spine[1], '.', c=colours[inull-1], ms=ms)
            if labels == True:
                plt.text(spine[2], spine[1], '{}'.format(inull), color='green')
            inull = np.asscalar(np.fromfile(spinefile, dtype=np.int32, count=1))

#################################################################

def separators(labels=False, ms=3):
    with io.open('output/'+prefile+'-separators-cut_'+rstr+'.dat', 'rb') as sepfile:
        null = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))
        while null > 0:
            sep = np.fromfile(sepfile, dtype=np.float64, count=3)
            plt.plot(sep[2], sep[1], '*', c='yellow', ms=ms)
            if labels == True:
                plt.text(sep[2], sep[1], '{}'.format(null))
            null = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))

    with io.open('output/'+prefile+'-hcs-separators-cut_'+rstr+'.dat', 'rb') as sepfile:
        null = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))
        while null > 0:
            sep = np.fromfile(sepfile, dtype=np.float64, count=3)
            plt.plot(sep[2], sep[1], '*', c='orange', ms=ms)
            null = np.asscalar(np.fromfile(sepfile, dtype=np.int32, count=1))

#################################################################

def rings(dots=False, labels=False, lw=1):
    with io.open('output/'+prefile+'-rings-cut_'+rstr+'.dat', 'rb') as ringfile:
        inull = np.asscalar(np.fromfile(ringfile, dtype=np.int32, count=1))
        while inull >= 0:
            length = np.asscalar(np.fromfile(ringfile, dtype=np.int32, count=1))
            ring = np.fromfile(ringfile, dtype=np.float64, count=3*length).reshape(-1,3)
            plt.plot(ring[:,2], ring[:,1], c=colours[inull-1], lw=lw)
            if labels == True:
                plt.text(ring[length//2,2], ring[length//2,1], '{}'.format(inull))
            if dots == True:
                plt.plot(ring[:,2], ring[:,1], '.', c='green')
            inull = np.asscalar(np.fromfile(ringfile, dtype=np.int32, count=1))

#################################################################

def hcs(dots=False, lw=1):
    with io.open('output/'+prefile+'-hcs-cut_'+rstr+'.dat', 'rb') as hcsfile:
        ihcs = np.asscalar(np.fromfile(hcsfile, dtype=np.int32, count=1))
        while ihcs >= 0:
            length = np.asscalar(np.fromfile(hcsfile, dtype=np.int32, count=1))
            # print(length)
            line = np.fromfile(hcsfile, dtype=np.float64, count=3*length).reshape(-1,3)
            plt.plot(line[:, 2], line[:, 1], ls=':', c='lime', lw=lw)
            if dots == True:
                plt.plot(line[:, 2], line[:, 1], '.', c='blue')
            ihcs = np.asscalar(np.fromfile(hcsfile, dtype=np.int32, count=1))

#################################################################

def nulls(labels=False, size=1):
    for i in range(nulldata.shape[0]):
        plt.plot(nulldata.pos[i, 2], nulldata.pos[i, 1], '.', color=colours[i], ms=size)
        if labels == True:
            plt.text(nulldata.pos[i,2], nulldata.pos[i,1], '{}, {:04f}'.format(i, nulldata.pos[i,0]))

#################################################################

def field():
    br, _, _, rads, thetas, phis = rd.field(datafile)
    levels = np.linspace(-10,10,101)
    plt.contourf(phis, thetas, br[0,:,:], levels, cmap=plt.cm.RdBu_r, extend='both', vmax=10, vmin=-10)
    cb = plt.colorbar(fraction=0.05, pad=0.025)
    cb.set_ticks([-10,-5,0,5,10])
    cb.set_label('Magnetic Field Strength (G)')
    plt.tight_layout()

#################################################################
def save(type, dpi=1200, figdir='figures', extra=None):
    if extra is not None:
        addin = '-' + extra
    else:
        addin = ''
    savename = figdir+'/'+prefile+'-cut_'+rstr+addin+'.{}'.format(type)
    plt.savefig(savename, dpi=dpi)
    print('Saved to {}'.format(savename))
