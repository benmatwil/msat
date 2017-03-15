from . import read as rd
import igraph as ig
import os

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
    for jnull in cnulls:
      if inull+1 in con[jnull-1]:
        g.add_edge(inull, jnull-1)
        con[jnull-1].remove(inull+1)
      else:
        g.add_edge(inull, jnull-1, color=g.vs[inull]['color'])
  
  if os.path.isfile('output/'+prefile+'-connectivity-hcs.dat'):
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

  if save == False:
    ig.plot(g, layout=layout, bbox=(1000,1000), margin=50, vertex_size=20)
  else:
    ig.plot(g, 'pngs/'+prefile+'-graph.png', layout=layout, bbox=(1000,1000), margin=50, vertex_size=20)

  return g
