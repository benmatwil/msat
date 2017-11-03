import numpy as np

nx = 200
ny = 200
nz = 200

grid = np.linspace(-0.1, 0.1, nx, dtype=np.float64)
x, y, z = np.meshgrid(grid, grid, grid, indexing='ij')

params = [(1, 0, 0, 0), (1.5, 0, 0, 0), (4, 0, 0, 0), (1, 0, 1, 0)]

for param in params:
    p, q, jpar, jperp = param
    bx = x + (q - jpar)*y/2
    by = (q + jpar)*x/2 + p*y
    bz = jperp*y - (p + 1)*z

    with open('data/nullfield{}.dat'.format(i+1), 'wb') as nf:
        np.array([nx, ny, nz], dtype=np.int32).tofile(nf)
        bx.T.tofile(nf)
        by.T.tofile(nf)
        bz.T.tofile(nf)
        grid.tofile(nf)
        grid.tofile(nf)
        grid.tofile(nf)

# bx = np.empty((nx,ny,nz), dtype=np.float64)
# by = bx.copy()
# bz = bx.copy()
# for ix, x in enumerate(grid):
#     for iy, y in enumerate(grid):
#         for iz, z in enumerate(grid):
#             bx[ix, iy, iz] = x + (q - jpar)*y/2
#             by[ix, iy, iz] = (q + jpar)*x/2 + p*y
#             bz[ix, iy, iz] = jperp*y - (p + 1)*z