import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Set matplotlib
plt.rcParams['font.family'] ='sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['figure.subplot.bottom'] = 0.15
plt.rcParams['figure.subplot.right'] = 0.85

# input file 
from data.indata import *

# analytical solution
data = np.loadtxt(fl_ana.strip(), skiprows=1, dtype="float")
x = data[:,0].reshape(ny,nx)
y = data[:,1].reshape(ny,nx)
u = data[:,2].reshape(ny,nx)
fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Analytical Solution")
c = ax.contourf(x, y, u, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(c, cax=cax).set_label('u(x,y)')
plt.show()

# numerical solution
data = np.loadtxt(fl_num.strip(), skiprows=1, dtype="float")
x = data[:,0].reshape(ny,nx)
y = data[:,1].reshape(ny,nx)
u = data[:,2].reshape(ny,nx)
fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Numerical Solution")
c = ax.contourf(x, y, u, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(c, cax=cax).set_label('u(x,y)')
plt.show()

# relative error
data = np.loadtxt(fl_err.strip(), skiprows=1, dtype="float")
x = data[:,0].reshape(ny,nx)
y = data[:,1].reshape(ny,nx)
err = data[:,2].reshape(ny,nx)
data = np.loadtxt(fl_num.strip(), skiprows=1, dtype="float")
fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Numerical Solution")
c = ax.contourf(x, y, err, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(c, cax=cax).set_label('u(x,y)')
plt.show()

err = err.reshape(n_n)
fig, ax = plt.subplots()
plt.hist(err)
ax.set_xscale('log')
ax.set_xlabel('Relative Error')
ax.set_ylabel('Frequency')
plt.show()
