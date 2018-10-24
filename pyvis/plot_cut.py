# example file to plot radial cuts through global magnetic field models of the solar corona
# can produce the majority of plots found in the Williams (2018)

import numpy as np
import matplotlib.pyplot as plt
from . import read as rd
import os
from random import shuffle
from fractions import Fraction as fr

plt.ion()

def get_colours():
    """
    Pick random colours for each null point in shades of red and blue
    """
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

def latex_ticks(nsplit, npi):
    """
    Produces strings for the ticks used for some spherical coordinate plots in
    latex for multiples of pi
    """
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
            string = string + r'\pi' + r'}{' + str(fract.denominator) + r'}$'
        else:
            if fract.numerator != 1:
                string = string + str(fract.numerator)
            string = string + r'\pi$'
        list.append(string)
    return list

#################################################################

def start(filename, norm, d, title=False, levels=np.linspace(-20,20,101), colmap=plt.cm.Greys_r, nticks=4, newcolours=True, linerast=False, output_dir=None):
    """
    Function to create the plot and prepare everything plotting using other functions
    """
    global datafile, prefile, nulldata, plane_norm, plane_d, colours

    datafile = filename
    prefile = rd.prefix(filename)

    if output_dir is not None:
        rd.outprefix = output_dir

    nulldata = rd.nulls(datafile, simple=True)

    plane_norm = norm
    try:
        len(d)
    except TypeError:
        plane_d = d
    else:
        pt = np.array(d, dtype=np.float64)
        plane_d = np.sum(pt*norm)

    # os.system(f'./make_cut -i {filename} -r {r}')

    plt.figure(figsize=(15/2.55, 6.6/2.55)) # set-up for A4
    ax = plt.gca()
    plt.xlabel(r'Longitude, \(\phi\)')
    plt.xlim([0, 2*np.pi])
    plt.xticks(np.linspace(0, 2*np.pi, nticks+1), latex_ticks(nticks, 2))
    plt.ylabel(r'Colatitude, \(\theta\)')
    plt.ylim([np.pi, 0])
    plt.yticks(np.linspace(0, np.pi, nticks+1), latex_ticks(nticks, 1))

    if linerast == False:
        ax.set_rasterization_zorder(-6)
    else:
        ax.set_rasterization_zorder(0)
        
    br, _, _, rads, thetas, phis = rd.field(datafile)
    ir = np.where(plane_d >= rads)[0].max()
    plotfield = br[ir, :, :] + (plane_d - rads[ir])/(rads[ir+1] - rads[ir])*(br[ir+1, :, :] - br[ir, :, :])
    plt.contourf(phis, thetas, plotfield, levels, cmap=colmap, extend='both', zorder=-10)
    cb = plt.colorbar(fraction=0.05, pad=0.025)
    diff = levels[-1] - levels[0]
    cb.set_ticks(np.linspace(0, 1.0, 5)*diff + levels[0])
    cb.set_label('Radial Field Strength (G)')
    plt.tight_layout()
    if title: plt.title('Cut: ' + datafile + ' r = ' + rstr)

    if newcolours: colours = get_colours()

#################################################################

def all(labels=False, ms=3):
    """
    Just calls all other functions for convenience
    """
    lw = 1
    rings(labels=labels, lw=lw)
    hcs(lw=lw)
    spines(labels=labels, ms=ms)
    separators(labels=labels, ms=ms)
    separators(labels=labels, ms=ms, hcs=True)

#################################################################

def spines(labels=False, ms=3):
    """
    Plots all the spine points
    """
    null_nums, spines = rd.cut_spines(datafile, plane_norm, plane_d)
    for inull, spine in zip(null_nums, spines):
        plt.plot(spine[2], spine[1], '.', c=colours[inull-1], ms=ms, zorder=-5)
        if labels:
            plt.text(spine[2], spine[1], '{}'.format(inull), color='green')

#################################################################

def separators(labels=False, ms=3, hcs=False):
    """
    Plots all the separator points
    """
    null_nums, separators = rd.cut_separators(datafile, plane_norm, plane_d, hcs=hcs)
    col = 'orange' if hcs else 'yellow'
    for inull, sep in zip(null_nums, separators):
        plt.plot(sep[2], sep[1], '*', c=col, ms=ms, zorder=-5)
        if labels:
            plt.text(sep[2], sep[1], '{}'.format(inull))

#################################################################

def rings(dots=False, labels=False, lw=1):
    """
    Plots all the lines denoting the ring crossings
    """
    null_nums, rings = rd.cut_sepsurf(datafile, plane_norm, plane_d)
    for inull, ring in zip(null_nums, rings):
        plt.plot(ring[:, 2], ring[:, 1], c=colours[inull-1], lw=lw, zorder=-5)
        if labels:
            plt.text(ring[ring.length[0]//2, 2], ring[ring.length[0]//2, 1], '{}'.format(inull))
        if dots:
            plt.plot(ring[:, 2], ring[:, 1], '.', c='green')

#################################################################

def hcs(dots=False, lw=1):
    """
    Plots all the hcs crossings
    """
    lines = rd.cut_hcs(datafile, plane_norm, plane_d)
    for line in lines:
        plt.plot(line[:, 2], line[:, 1], c='lime', ls=':', lw=lw, zorder=-0.5)
        if dots:
            plt.plot(line[:, 2], line[:, 1], '.', c='blue')

    # with open(outprefix+'/'+prefile+'-hcs-cut_'+rstr+'.dat', 'rb') as hcsfile:
    #     nlines, = np.fromfile(hcsfile, dtype=np.int32, count=1)
    #     print(nlines)
    #     if rad > 2.49:
    #         skip = 2
    #     else:
    #         skip = 1
    #     for _ in range(0, nlines, skip):
    #         ihcs, = np.fromfile(hcsfile, dtype=np.int32, count=1)
    #         length, = np.fromfile(hcsfile, dtype=np.int32, count=1)
    #         line = np.fromfile(hcsfile, dtype=np.float64, count=3*length).reshape(-1,3)

#################################################################

def nulls(labels=False, size=1, zorder=-1):
    """
    Plots all the null point locations in theta/phi
    """
    for i in range(nulldata.shape[0]):
        plt.plot(nulldata.pos[i, 2], nulldata.pos[i, 1], '.', color=colours[i], ms=size, zorder=zorder)
        if labels == True:
            plt.text(nulldata.pos[i,2], nulldata.pos[i,1], '{}, {:04f}'.format(i, nulldata.pos[i,0]))

#################################################################

def field():
    """
    Plots the background radial field on the solar surface
    """
    br, _, _, rads, thetas, phis = rd.field(datafile)
    levels = np.linspace(-10, 10, 101)
    plt.contourf(phis, thetas, br[0, :, :], levels, cmap=plt.cm.RdBu_r, extend='both', vmax=10, vmin=-10, zorder=-10)
    cb = plt.colorbar(fraction=0.05, pad=0.025)
    cb.set_ticks([-10, -5, 0, 5, 10])
    cb.set_label('Magnetic Field Strength (G)')
    plt.tight_layout()
