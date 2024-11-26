import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from functions import recenter_image

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
sr_to_arcsec2 = 4.25e+10
xx_recenter, yy_recenter = 8/40,28/40

recomb_path = "data/march/RR_a1.5_f1.5.npy"

start_idx = 15
SB_Recomb_intrinsic = np.load(recomb_path)/sr_to_arcsec2 #need to convert from 1/sr to 1/arcsec2
SB_Recomb_intrinsic = np.sum(SB_Recomb_intrinsic[:,start_idx:start_idx+index,:],axis=1)
SB_Recomb_intrinsic = recenter_image(SB_Recomb_intrinsic,xx_recenter, yy_recenter,Ngas)


vmax = 8
multiplier=1e21
fig,ax = plt.subplots(1,figsize=(6,6))
im = ax.imshow((SB_Recomb_intrinsic)*multiplier,origin='lower', vmin=0,vmax=vmax,cmap=cm.inferno)
ax_divider = make_axes_locatable(ax)
cax = ax_divider.append_axes("top", size="7%", pad="2%")
cb = fig.colorbar(im,cax=cax,ticklocation = 'top',orientation="horizontal",pad=-0.02)
cb.set_label("Recomb. Surface Brightness\n"+r"[10$^{-21} $erg s$^{-1}$ cm $^{-2}$ arcsec$^{-2}$]",fontsize=fontsize+2)
cax.xaxis.set_ticks_position("top")
# cb.ax.tick_params(direction='out')
# cb.locator_params(nbins=5)

# cb.ax.set_label(r"Surface Brightness [10$^{-21} $erg s$^{-1}$ cm $^{-2}$ arcsec$^{-2}$]")#,y=0.94)
#
ticks= np.linspace(0,Ngas,5)
labels = np.asarray(np.linspace(0,40,5),int)
ax.set_xticks(ticks=ticks)
ax.set_xticklabels(labels)
ax.set_yticks(ticks=ticks)
ax.set_yticklabels(labels)
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.set_xlim(0,Ngas)
ax.set_xlabel("cMpc/h",fontsize=fontsize+1)
ax.set_ylabel("cMpc/h",fontsize=fontsize+1)

plt.tight_layout()
plt.savefig("figures/recomb_only.png",bbox_inches='tight',dpi=300)


# for i in range(1,4): #looping over rows
#     for j in range(1,5): #looping over columns
#         if ((i==1)&(j==1)): #this will be the SB CDF plot
#             pass
#         else:
#             ax[i][j].set_xticks(ticks=ticks)
#             ax[i][j].set_xticklabels(labels)
#             ax[i][j].set_yticks(ticks=ticks)
#             ax[i][j].set_yticklabels(labels)
#             ax[i][j].scatter(x=[Ngas/2],y=[Ngas/2],marker="+",color="red")
#             ax[i][j].tick_params(bottom=True, top=True, left=True, right=True)
#             ax[i][j].set_xlim(0,Ngas)




#
# # slice stack is a file created on expanse: FlexRT/postprocessing/extract_results.py
# Delta_slice, xHII_slice, F_lyC_inc_slice, vIF_slice, F_lya_slice = np.load("data/march/APR_slice_stack.npy")
# xHI_slice  = 1-xHII_slice
# z=5.7
# omega_b = 0.048
# rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
# rhoBaryons_0 = rho_crit_z*omega_b
# Y=0.24
# m_H = 1.672622e-24 #g
#
# n=1
# m=1
# bottom = 0.1; left=0.06
# top=1.-bottom; right = 1.-left
# figheight = 4.5 #inch
# figwidth = figheight*1.5
#
# fig, ax = plt.subplots(nrows=n, ncols=m, figsize=(figwidth, figheight),sharex=True,sharey=True,dpi=800)
# plt.subplots_adjust(bottom=bottom,top=top, left=left, right=right)
#
# ticks = np.linspace(0,Ngas,5)
# labels = np.asarray(np.linspace(0,40,5),int)
# ax.set_xticks(ticks=ticks)
# ax.set_xticklabels(labels)
# ax.set_yticks(ticks=ticks)
# ax.set_yticklabels(labels)
# ax.xaxis.set_ticks_position('both')
# ax.yaxis.set_ticks_position('both')
# ax.xaxis.set_ticks_position('both')
# ax.yaxis.set_ticks_position('both')
# ax.set_xlim(0,Ngas)
# ax.set_ylim(0,Ngas)
# ax.set_ylabel("cMpc/h")
#
# #plotting density field. Note that density is in log-space
# im = ax.imshow(np.log10(Delta_slice),cmap=cm.inferno,vmin=-0.5,vmax=1.5,origin='lower')
# ax_divider = make_axes_locatable(ax)
# cax = ax_divider.append_axes("top", size="7%", pad="2%")
# cb = fig.colorbar(im,cax=cax, orientation="horizontal",pad=-0.02)
# cax.xaxis.set_ticks_position("top")
#
# #neutral region is chosen to be where gas is less than 10% ionized (more than 90% neutral)
# xHII_slice = 1-xHI_slice
# xHII_slice[xHII_slice>0.1]=np.nan #pixels where gas is more than 10% ionized are set to be nan (unshaded)
# ax.imshow(xHII_slice,cmap=cm.binary,vmin=0,vmax=1,origin='lower',alpha=0.5)
# pad_title = 50
# ax.set_title(r"Overdensity, log$_{10}\Delta$", pad=pad_title)
