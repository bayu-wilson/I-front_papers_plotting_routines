import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

#constants
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
sr_to_arcsec2 = 4.25e10
E_lya = h*c/lambda_lya_cm

#plotting stuff
fontsize=13
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
z=5.7

slice_stack = np.load("data/march/APR_slice_stack.npy")
xchunk = [5/40,15/40]
ychunk = [5/40,15/40]
xlim = [int(xchunk[0]*Ngas),int(xchunk[1]*Ngas)]
ylim = [int(ychunk[0]*Ngas),int(ychunk[1]*Ngas)]
Ngas = xlim[1]-xlim[0]
slice_stack = slice_stack[:,ylim[0]:ylim[1],xlim[0]:xlim[1]]

Delta_slice, xHII_slice, F_lyC_inc_slice, vIF_slice, F_lya_slice = slice_stack
xHI_slice = 1-xHII_slice

slice_list = [np.log10(Delta_slice),F_lyC_inc_slice, vIF_slice,F_lya_slice]
cmap_list = [cm.inferno, cm.gnuplot, cm.gnuplot, cm.gnuplot]

n=1 #rows
m=len(slice_list) #columns
bottom = 0.1; left=0.06
top=1.-bottom; right = 1.-left
figheight = 5
figwidth = figheight*2.5

fig, ax = plt.subplots(nrows=n, ncols=m, figsize=(figwidth, figheight),sharex=True,sharey=True,dpi=500)
ax = ax.flatten()
plt.subplots_adjust(bottom=bottom,top=top, left=left, right=right)


vrange = [[-0.5,1.5],[0,105_000],[0,25000],[0,105_000/5]]
for i in range(m): #looping over columns
	im = ax[i].imshow(slice_list[i],cmap=cmap_list[i],vmin=vrange[i][0],vmax=vrange[i][1],origin='lower')
	ax_divider = make_axes_locatable(ax[i])
	cax = ax_divider.append_axes("top", size="7%", pad="2%")
	cb = fig.colorbar(im,cax=cax, orientation="horizontal",pad=-0.02)
	cax.xaxis.set_ticks_position("top")

#overplotting the neutral gas as the shaded region in the overdensity panel
xHII_slice[xHII_slice>0.1]=np.nan #ignoring all regions that are ionized
ax[0].imshow(xHII_slice,cmap=cm.binary,vmin=0,vmax=1,origin='lower',alpha=0.5)

Nticks = 5
ticks= np.linspace(0,Ngas,Nticks)
xlabels = np.asarray(np.linspace(xchunk[0]*40,xchunk[1]*40,Nticks)+0.5,int)
ylabels =  np.asarray(np.linspace(ychunk[0]*40,ychunk[1]*40,Nticks)+0.5,int)

for j in range(m):
    axis = ax[j]
	#option to add a red dot in the center
    # axis.scatter(x=[Ngas/2],y=[Ngas/2],marker="+",color="red"))
    axis.set_xticks(ticks=ticks)
    axis.set_xticklabels(xlabels)
    axis.set_yticks(ticks=ticks)
    axis.set_yticklabels(ylabels)

    axis.xaxis.set_ticks_position('both')
    axis.yaxis.set_ticks_position('both')
    axis.xaxis.set_ticks_position('both')
    axis.yaxis.set_ticks_position('both')

    axis.set_xlim(0,Ngas)
    axis.set_ylim(0,Ngas)
ax[0].set_ylabel("cMpc/h")
ax[3].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True,
                 bottom=True, top=True, left=True, right=True,direction='out')

pad_title = 55
ax[0].set_title("Overdensity\n"+r"log$_{10}\Delta$", pad=pad_title)
ax[1].set_title("Inc. LyC Flux\n"+r"[s$^{-1}$ cm$^{-2}$]", pad=pad_title)
ax[2].set_title("IF speed\n"+r"[km s$^{-1}$]", pad=pad_title)
ax[3].set_title(r"Emitted Ly$\alpha$ Flux "+"\n"+r"[s$^{-1}$ cm$^{-2}$]", pad=pad_title)

plt.savefig("figures/zoom_slice.png",bbox_inches='tight')
