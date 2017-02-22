# import matplotlib.pyplot as plt
import read_data as rd
# import numpy as np
import igraph as ig

filename = 'data/field_20100601_0081-fft.dat'
# filename = 'data/field_test.dat'

def plot(filename, layout=None):
  nulls = rd.nulls(filename)
  con = rd.connectivity(filename)

  g = ig.Graph()
  for inull in range(nulls.shape[0]):
    if nulls[inull].sign == 1:
      col = 'red'
    elif nulls[inull].sign == -1:
      col = 'blue'
    else:
      col = 'green'
    g.add_vertex(color=col)

  for inull, cnulls in enumerate(con):
    for jnull in cnulls:
      g.add_edge(inull, jnull-1, color=g.vs[inull]['color'])


  # layout = g.layout('kk')
  if layout is None:
    layout = g.layout('fr')
  else:
    layout = g.layout(layout)

  ig.plot(g, layout=layout, bbox=(1000,1000), margin=50, vertex_size=20)

# g = nx.Graph()
# cols = {}
# for inull in nulls.number:
#   g.add_node(inull)
#   # cols[inull] = 'blue'

# for inull, cnulls in enumerate(con):
#   for jnull in cnulls:
#     g.add_edge(inull+1, jnull)

# pos = nx.spring_layout(g, iterations=200, k=0.5)
# # pos = nx.pydotplus.graphviz_layout(g, prog='neato')

# nx.draw_networkx_nodes(g, pos=pos, nodelist=(np.where(nulls.sign==1)[0]+1).tolist(), node_color='r')
# nx.draw_networkx_nodes(g, pos=pos, nodelist=(np.where(nulls.sign==-1)[0]+1).tolist(), node_color='b')
# nx.draw_networkx_nodes(g, pos=pos, nodelist=(np.where(nulls.sign==0)[0]+1).tolist(), node_color='g')

# nx.draw_networkx_edges(g, pos=pos)

# nx.draw_networkx_labels(g, pos=pos, font_size=10)

# adjmat = nx.to_numpy_matrix(g).astype(np.int32)

# na = rd.nulls('data/field_20100601_0081-anal.dat')
# nf = rd.nulls('data/field_20100601_0081-fft.dat')

# xa = na.pos[:,0]*np.sin(na.pos[:,1])*np.cos(na.pos[:,2])
# ya = na.pos[:,0]*np.sin(na.pos[:,1])*np.sin(na.pos[:,2])
# za = na.pos[:,0]*np.cos(na.pos[:,1])

# xf = nf.pos[:,0]*np.sin(nf.pos[:,1])*np.cos(nf.pos[:,2])
# yf = nf.pos[:,0]*np.sin(nf.pos[:,1])*np.sin(nf.pos[:,2])
# zf = nf.pos[:,0]*np.cos(nf.pos[:,1])

# dists = np.zeros((xa.shape[0], xf.shape[0]), dtype=np.float64)

# for i in xrange(xa.shape[0]):
#   for j in xrange(xf.shape[0]):
#     dists[i,j] = np.sqrt((xa[i] - xf[j])**2 + (ya[i] - yf[j])**2 + (za[i] - zf[j])**2)

# list = np.zeros(xa.shape[0], dtype=np.int32)-1
# for i in xrange(xa.shape[0]):
#   if dists[i,:].min() < 0.01:
#     list[i] = dists[i,:].argmin()




# br, bt, bp, rg, tg, pg = rd.field('data/field_20100601_0081-fft.dat')

# shp = (2*rg.shape[0]-1, 2*tg.shape[0]-1, 2*pg.shape[0]-1, 3)

# bgrid = np.zeros(shp, dtype=np.float64)

# bgrid[::2,::2,::2,0] = br
# bgrid[::2,::2,::2,1] = bt
# bgrid[::2,::2,::2,2] = bp

# bgrid[1::2,:,:,:] = (bgrid[:-2:2,:,:,:] + bgrid[2::2,:,:,:])/2
# bgrid[:,1::2,:,:] = (bgrid[:,:-2:2,:,:] + bgrid[:,2::2,:,:])/2
# bgrid[:,:,1::2,:] = (bgrid[:,:,:-2:2,:] + bgrid[:,:,2::2,:])/2

# rgrid = np.zeros(2*rg.shape[0]-1, dtype=np.float64)
# tgrid = np.zeros(2*tg.shape[0]-1, dtype=np.float64)
# pgrid = np.zeros(2*pg.shape[0]-1, dtype=np.float64)

# rgrid[::2] = rg
# rgrid[1::2] = (rg[:-1] + rg[1:])/2
# tgrid[::2] = tg
# tgrid[1::2] = (tg[:-1] + tg[1:])/2
# pgrid[::2] = pg
# pgrid[1::2] = (pg[:-1] + pg[1:])/2

# with open('data/field_20100601_0081-tri.dat', 'wb') as file:
#   np.array([bgrid.shape[0], bgrid.shape[1], bgrid.shape[2]], np.int32).tofile(file)
#   bgrid.T.tofile(file)
#   rgrid.tofile(file)
#   tgrid.tofile(file)
#   pgrid.tofile(file)