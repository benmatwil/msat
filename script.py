import pyvis.model3d as m3d

m3d.make('data/mhs_field_20100601_0161_fft_alpha_0.10_d_0.50-finite.dat', [], coordsystem='spherical')

ntheta = 6
nphi = 12
phis = np.linspace(0, 2*np.pi, nphi, endpoint=False)
thetas = np.linspace(0, np.pi, ntheta, endpoint=False)
thetas = thetas[1:]
rl(m3d.fl)
# for r, col in zip([1.2],[(1,0,0)]):
for r, col in zip([1.05,1.2,2.4],[(0,0,1),(1,0,0),(0,1,0)]):
    pts = []
    for it, t in enumerate(thetas):
        for ip, p in enumerate(phis):
            pts.append(np.array([r, t, p]))
    m3d.add_fieldlines(pts, col=col)

x, z = np.mgrid[0:2.5:10j, -2.5:2.5:10*1j]
y = np.ones_like(x)
y = np.zeros_like(x)
m3d.ml.mesh(x, y, z, color=(0,0,0))

m3d.make('data/mhs_field_20100601_0081_fft_alpha_1.00_d_0.00-infinite.dat', [], coordsystem='spherical')

ntheta = 5
nphi = 8
phis = np.linspace(0, 2*np.pi, nphi, endpoint=False)
thetas = np.linspace(0, np.pi, ntheta, endpoint=False)
thetas = thetas[1:-1]
for r, col in zip([1.2,2.4],[(1,0,0),(0,1,0)]):
    pts = []
    for it, t in enumerate(thetas):
        for ip, p in enumerate(phis):
            pts.append(np.array([r, t, p]))
    m3d.add_fieldlines(pts, col=col)

x, z = np.mgrid[0:2.5:10j, -2.5:2.5:10*1j]
y = np.ones_like(x)
y = np.zeros_like(x)
m3d.ml.mesh(x, y, z, color=(0,0,0))
