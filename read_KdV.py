import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
def format(fig):
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['xtick.labelsize'] = 19
    plt.rcParams['ytick.labelsize'] = 19
    plt.rcParams['font.size'] = 19
    plt.rcParams['figure.figsize'] = [5.6*6, 4*3]
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams['legend.fontsize'] = 15
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.linewidth'] = 1.5
    # plt.style.use('dark_background')


def ax_format(ax, xmaj, xmin, ymaj, ymin):
    ax.xaxis.set_tick_params(which='major', size=5, width=1,
                            direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                            direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=5, width=1,
                            direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                            direction='in', right='on')
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xmin))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ymaj))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ymin))

path = 'x_explicit.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
x_explicit = np.array(list)  

path = 'rho_explicit_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_explicit = np.array(list)

path = 'rho_theo_explicit_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_theo_explicit = np.array(list)

path = 'x_implicit.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
x_implicit = np.array(list)

path = 'rho_implicit_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_implicit = np.array(list)

path = 'rho_theo_implicit_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_theo_implicit = np.array(list)

path = 'x_explicit_matrix.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
x_explicit_matrix = np.array(list)

path = 'rho_explicit_matrix_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_explicit_matrix = np.array(list)

path = 'rho_theo_explicit_matrix_plot.txt'
data_list = []
list = []
fo = open(path, 'r')
for line in fo:
    s = line.split(',')
    for i in range(len(s)):
        data_list.append(float(s[i]))
    list.append(data_list)
    data_list = []
fo.close()
rho_theo_explicit_matrix = np.array(list)

fig1 = plt.figure(figsize=(16,8))
format(fig1)
fig1.suptitle('Solve Korteweg-de Vries (KdV) equation')

ax1 = fig1.add_subplot(1, 3, 1)
ax1.plot(x_explicit, rho_theo_explicit[:,0], '--', label = 'Initial')
ax1.plot(x_explicit, rho_theo_explicit[:,-1], label = 'Exact')
ax1.plot(x_explicit, rho_explicit[:,-1], '+', label = 'Computed')
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$\rho (x,t)$')
ax1.set_xlim(x_explicit[0], x_explicit[-1])
ax1.legend()
ax1.set_title('Use iteration explicit scheme')

ax1 = fig1.add_subplot(1, 3, 2)
ax1.plot(x_explicit_matrix, rho_theo_explicit_matrix[:,0], '--', label = 'Initial')
ax1.plot(x_explicit_matrix, rho_theo_explicit_matrix[:,-1], label = 'Exact')
ax1.plot(x_explicit_matrix, rho_explicit_matrix[:,-1], '+', label = 'Computed')
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$\rho (x,t)$')
ax1.set_xlim(x_explicit_matrix[0], x_explicit_matrix[-1])
ax1.legend()
ax1.set_title('Use matrix explicit scheme')

ax1 = fig1.add_subplot(1, 3, 3)
ax1.plot(x_implicit, rho_theo_implicit[:,0], '--', label = 'Initial')
ax1.plot(x_implicit, rho_theo_implicit[:,-1], label = 'Exact')
ax1.plot(x_implicit, rho_implicit[:,-1], '+', label = 'Computed')
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$\rho (x,t)$')
ax1.set_xlim(x_implicit[0], x_implicit[-1])
ax1.legend()
ax1.set_title('Use matrix implicit scheme (Crank-Nicolson)')

plt.tight_layout()
#plt.show()
plt.savefig('KdV.png')
  

