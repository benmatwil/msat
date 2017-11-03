from fractions import Fraction as fr
import pyvis.read as rd

br, ps, ts = rd.synmap('../sphericalharmonics/hmi/synmap_20100601.dat')

ts = np.arccos(ts)

def ticks(nsplit, npi):
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

levels=np.linspace(-20,20,101)
colmap=plt.cm.RdBu_r
nticks=4

plt.figure(figsize=(15/2.55, 7.2/2.55)) # set-up for A4
ax = plt.gca()
plt.xlabel('Longitude')
plt.xlim([0, 2*np.pi])
plt.xticks(np.linspace(0, 2*np.pi, nticks+1), ticks(nticks, 2))
plt.ylabel('Latitude')
plt.ylim([np.pi, 0])
plt.yticks(np.linspace(0, np.pi, nticks+1), ticks(nticks, 1))

ax.set_rasterization_zorder(0)
plt.contourf(ps, ts, br, levels, cmap=colmap, extend='both', zorder=-10)
cb = plt.colorbar(fraction=0.05, pad=0.025)
diff = levels[-1] - levels[0]
cb.set_ticks(np.array([0.0, 0.25, 0.5, 0.75, 1.0])*diff + levels[0])
# print(np.array([0.0, 0.25, 0.5, 0.75, 1.0])*diff + levels[0])
cb.set_label('Magnetic Field Strength (G)')
plt.tight_layout()
