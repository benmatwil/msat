import matplotlib.pyplot as plt
import numpy as np
import read_data as rd

filename = 'data/field_20100601_0081-fft.dat'

nulls = rd.nulls(filename)
seps = rd.separators(filename)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make data
u = np.linspace(0, 2 * np.pi, 200)
v = np.linspace(0, np.pi, 200)
rad = 1
x = rad * np.outer(np.cos(u), np.sin(v))
y = rad * np.outer(np.sin(u), np.sin(v))
z = rad * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='yellow')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
rad = 0.01
x = rad * np.outer(np.cos(u), np.sin(v))
y = rad * np.outer(np.sin(u), np.sin(v))
z = rad * np.outer(np.ones(np.size(u)), np.cos(v))

for null in nulls:
  if null.sign == 1:
    col = 'red'
  elif null.sign == -1:
    col = 'blue'
  else:
    col = 'green'
  pos = null.pos
  pos = np.array([pos[0]*np.sin(pos[1])*np.cos(pos[2]), 
    pos[0]*np.sin(pos[1])*np.sin(pos[2]),
    pos[0]*np.cos(pos[1])])
  # ax.plot_surface(x+pos[0], y+pos[1], z+pos[2], color=col, shade=False)
  for isep in range(len(seps[null.number-1])):
    pts = seps[null.number-1][isep]
    xpts = pts[:,0]*np.sin(pts[:,1])*np.cos(pts[:,2])
    ypts = pts[:,0]*np.sin(pts[:,1])*np.sin(pts[:,2])
    zpts = pts[:,0]*np.cos(pts[:,1])
    ax.plot(xpts, ypts, zpts, c='green')

plt.draw()

filename = 'data/field_20100601_0081-fft.dat'

nulls = rd.nulls(filename)
seps = rd.separators(filename)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for null in nulls:
  if null.sign == 1:
    col = 'red'
  elif null.sign == -1:
    col = 'blue'
  else:
    col = 'green'
  pos = null.pos
  # pos = np.array([pos[0]*np.sin(pos[1])*np.cos(pos[2]), 
  #  pos[0]*np.sin(pos[1])*np.sin(pos[2]),
  #  pos[0]*np.cos(pos[1])])
  # ax.plot_surface(x+pos[0], y+pos[1], z+pos[2], color=col, shade=False)
  for isep in range(len(seps[null.number-1])):
    pts = seps[null.number-1][isep]
    xpts = pts[:,0]*np.sin(pts[:,1])*np.cos(pts[:,2])
    ypts = pts[:,0]*np.sin(pts[:,1])*np.sin(pts[:,2])
    zpts = pts[:,0]*np.cos(pts[:,1])
    ax.plot(xpts, ypts, zpts, c='green')
    # ax.plot(pts[:,0], pts[:,1], pts[:,2], c='green')

plt.draw()




