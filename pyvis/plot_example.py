import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

X = np.linspace(-np.pi, np.pi, 256, endpoint=True)
C, S = np.cos(X), np.sin(X)

plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})

plt.figure(figsize=(10, 6))

plt.plot(X, C, color="blue", linewidth=2.5, linestyle="-", label=r"$\cos(x)$")
plt.plot(X, S, color="red",  linewidth=2.5, linestyle="-", label=r"$\sin(x)$")

plt.legend(loc='upper left')

plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
         [r'$-\pi$', r'$-\frac{\pi}{2}$', r'$0$', r'$+\frac{\pi}{2}$', r'$+\pi$'])

plt.yticks([-1, 0, +1],
         [r'$-1$', r'$0$', r'$+1$'])

ax = plt.gca()  # gca stands for 'get current axis'
ax.spines['right'].set_color('none') # should this really be hide rather than no colour
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position('zero')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position('zero')

t = 2 * np.pi / 3
plt.plot([t, t], [0, np.cos(t)], color='blue', linewidth=2.5, linestyle="--")
plt.scatter([t, ], [np.cos(t), ], 50, color='blue')

plt.annotate(r'$\sin\left(\frac{2\pi}{3}\right)=\frac{\sqrt{3}}{2}$',
             xy=(t, np.sin(t)), xycoords='data',
             xytext=(+10, +30), textcoords='offset points', fontsize=16,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

plt.plot([t, t],[0, np.sin(t)], color='red', linewidth=2.5, linestyle="--")
plt.scatter([t, ],[np.sin(t), ], 50, color='red')

plt.annotate(r'$\cos\left(\frac{2\pi}{3}\right)=-\frac{1}{2}$',
             xy=(t, np.cos(t)), xycoords='data',
             xytext=(-90, -50), textcoords='offset points', fontsize=16,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
    label.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0))

plt.savefig('figure_1.png',dpi=240)