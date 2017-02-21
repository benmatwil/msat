map = rd.synmap('hmi/synmap_20100601.dat')

plt.figure(figsize=(11,5))
plt.contourf(map[1], map[2], map[0], np.linspace(-10,10,21), cmap=plt.cm.RdBu_r, extend='both')
#plt.xticks([0, 2*np.pi/3, 4*np.pi/3, 2*np.pi], ['0', r'$\displaystyle \frac{2\pi}{3}$', r'$\displaystyle \frac{4\pi}{3}$', r'$2\pi$'])
plt.xticks([0, 2*np.pi/3, 4*np.pi/3, 2*np.pi], ['0', r'$\frac{2\pi}{3}$', r'$\frac{4\pi}{3}$', r'$2\pi$'])
plt.yticks([-1,-0.5,0,0.5,1])
plt.xlabel('Longitude')
plt.ylabel('Sine Latitude')

cb = plt.colorbar(fraction=0.05, pad=0.025)

cb.set_ticks([-10,-5,0,5,10])
cb.set_label('Magnetic Field Strength (G)')

plt.tight_layout()
plt.savefig('/user/bmw/Pictures/baddata.png')

