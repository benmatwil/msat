from . import read as rd
import igraph as ig
import os
from math import *

def plot(filename, layout=None, labels=False, save=False):
    print('Filename is {}'.format(filename))
    nulls = rd.nulls(filename)
    con = rd.separators(filename, lines=False)
    prefile = rd.prefix(filename)

    g = ig.Graph()
    g.add_vertices(nulls.number[-1])
    for inull in range(nulls.number[-1]):
      if nulls[inull].sign != 0:
        if nulls[inull].sign == 1:
          col = 'red'
        elif nulls[inull].sign == -1:
          col = 'blue'
        g.vs[inull]['color'] = col
        if labels == True:
          g.vs[inull]['label'] = '{}'.format(inull+1)

    for inull, cnulls in enumerate(con):
        toremove = []
        for jnull in cnulls:
            if inull+1 in con[jnull-1]:
                g.add_edge(inull, jnull-1)
                con[jnull-1].remove(inull+1)
                toremove.append(jnull)
                # con[inull].remove(jnull)
            else:
                g.add_edge(inull, jnull-1, color=g.vs[inull]['color'])
        for iremove in toremove:
            con[inull].remove(iremove)
    # for inull, cnulls in enumerate(con):
    #   print(inull+1, cnulls)
    #   for jnull in cnulls:
    #     g.add_edge(inull, jnull-1, color=g.vs[inull]['color'])

    if os.path.isfile('output/'+prefile+'-hcs-connectivity.dat'):
        conhcs = rd.separators(filename, lines=False, hcs=True)
        for ihcs, connect in enumerate(conhcs):
            if len(connect) > 0:
                g.add_vertex(1)
                g.vs[ihcs+nulls.number[-1]]['color'] = 'green'
                for inull in connect:
                    g.add_edge(ihcs+nulls.number[-1], inull-1)

    # layout = g.layout('kk')
    if layout is None:
        layout = g.layout('fr')
    else:
        layout = g.layout(layout)

    vsize = 10

    if save == False:
        ig.plot(g, layout=layout, bbox=(1000,1000), margin=50, vertex_size=vsize)
    else:
        ig.plot(g, 'figures/'+prefile+'-graph.png', layout=layout, bbox=(1000,1000), margin=50, vertex_size=vsize)

    # return g
