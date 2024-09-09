import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

fontsize=15
tick_dir='out'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.prop_cycle'] = cycler(color=['k','r','b'])
plt.rcParams['axes.prop_cycle'] = cycler(ls=['-'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size

Ngas=400
nb_depth = 30 #mpc/h
sim_depth = 40 #mpc/h
index = int(nb_depth/sim_depth*Ngas+0.5)
index_half = int(index/2)

# slice stack is a file created on expanse: FlexRT/postprocessing/extract_results.py
Delta_slice, xHII_slice, F_lyC_inc_slice, vIF_slice, F_lya_slice = np.load("data/march/APR_slice_stack.npy")
xHI_slice  = 1-xHII_slice
z=5.7
omega_b = 0.048
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rhoBaryons_0 = rho_crit_z*omega_b
Y=0.24
m_H = 1.672622e-24 #g

n=1
m=1
bottom = 0.1; left=0.06
top=1.-bottom; right = 1.-left
figheight = 4.5 #inch
figwidth = figheight*1.5

fig, ax = plt.subplots(nrows=n, ncols=m, figsize=(figwidth, figheight),sharex=True,sharey=True,dpi=800)
plt.subplots_adjust(bottom=bottom,top=top, left=left, right=right)

ticks = np.linspace(0,Ngas,5)
labels = np.asarray(np.linspace(0,40,5),int)
ax.set_xticks(ticks=ticks)
ax.set_xticklabels(labels)
ax.set_yticks(ticks=ticks)
ax.set_yticklabels(labels)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,Ngas)
ax.set_ylim(0,Ngas)
ax.set_ylabel("cMpc/h")

#plotting density field. Note that density is in log-space
im = ax.imshow(np.log10(Delta_slice),cmap=cm.inferno,vmin=-0.5,vmax=1.5,origin='lower')
ax_divider = make_axes_locatable(ax)
cax = ax_divider.append_axes("top", size="7%", pad="2%")
cb = fig.colorbar(im,cax=cax, orientation="horizontal",pad=-0.02)
cax.xaxis.set_ticks_position("top")

#neutral region is chosen to be where gas is less than 10% ionized (more than 90% neutral)
xHII_slice = 1-xHI_slice
xHII_slice[xHII_slice>0.1]=np.nan #pixels where gas is more than 10% ionized are set to be nan (unshaded)
ax.imshow(xHII_slice,cmap=cm.binary,vmin=0,vmax=1,origin='lower',alpha=0.5)
pad_title = 50
ax.set_title(r"Overdensity, log$_{10}\Delta$", pad=pad_title)

plt.savefig("figures/gas_slice_shaded.png",bbox_inches='tight')
