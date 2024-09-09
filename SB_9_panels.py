import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib.gridspec import GridSpec
from functions import recenter_image

dataCubeDir = "data/"

#path-lists for the SB data files
pathList_beforeRT_collExc = [
    dataCubeDir+"march/SB_a0.5_f7.3.npy",
    dataCubeDir+"march/SB_a1.5_f7.3.npy",
    dataCubeDir+"march/SB_a1.5_f1.5.npy"]

pathlist_beforeRT_Recomb = [
    dataCubeDir+"march/RR_a1.5_f1.5.npy",
    dataCubeDir+"march/RR_a1.5_f1.5.npy",
    dataCubeDir+"march/RR_a1.5_f1.5.npy"]

pathlist_afterRT_both = [
    dataCubeDir+"march/out_cube_Nphot600000_hard_moreRapid_wRecomb.npy",
    dataCubeDir+"march/out_cube_Nphot600000_soft_moreRapid_wRecomb.npy",
    dataCubeDir+"march/out_cube_Nphot600000_soft_lessRapid_wRecomb.npy"]

#combining the three lists: (1) coll. excitations only (2) recombinations only (3) combined maps with Lya RT applied in post-processing
all_paths = [pathList_beforeRT_collExc,pathlist_beforeRT_Recomb,pathlist_afterRT_both]

# setting parameters related to NB816 filter
Ngas = 400
nb_depth = 26 #mpc/h # assuming NB816 filter at z=5.7
sim_depth = 40 #mpc/h
index = int(nb_depth/sim_depth*Ngas+0.5) #used to integrate the SB's over the filter from "a" to "a+index" where "a" is the initial pixel.
sr_to_arcsec2 = 4.25e+10
cell_size = 40/Ngas*1000/0.68/(1+5.7)*3.086e+21 #in proper centimeters

# The datacubes are integrated along axis "orientation" (in the pre-RT data format)
orientation = "1" #note that this is for book-keeping. If you want to integrate over another axis, you must do it manually.

# loading and recentering neutral fraction and density data files
avg_neutral_fraction = np.load("data/march/FEB_avg_neutral_fraction.npy")
avg_Delta = np.load("data/march/FEB_avg_Delta.npy")
xx_recenter, yy_recenter = 8/40,28/40
avg_neutral_fraction = recenter_image(avg_neutral_fraction,xx_recenter, yy_recenter,Ngas)
avg_Delta = recenter_image(avg_Delta,xx_recenter, yy_recenter,Ngas)
# I manually recenter the image by a factor of 8 Mpc rightward and 28 mpc upward (equivalently 12 Mpc downward)
# This translation centers the map on the highly neutral region
# note that the origin is set to the bottom left corner of the maps

#### Make sure everything above this line is correct before running ###
#setting plotting parameters
fontsize=15
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
color_list = ["r","b","g"]

#these lists will eventually contain the 2D SB maps
SB_mocks_collExc = []
SB_mocks_recomb = []
SB_mocks_afterRT_combined = []
flattened_SB_mocks = np.array([])
flattened_afterRT_SB_mocks = np.array([])
im_list = []
vmax = 8
multiplier=1e21 #multiplied to the SB so we get numbers between 1 & 10

alpha_list = ["0.5","1.5","1.5"]
reion_hist_list = [r"FASTER-v$_\mathrm{IF}$",r"FASTER-v$_\mathrm{IF}$",r"SLOWER-v$_\mathrm{IF}$"]
N_cubes_per_row = 3
ncols = 1 + N_cubes_per_row #adding 1 because the first column shows the CDF, denstiy field, and ionization field.

#Making the 9 panel plot
fig = plt.figure(figsize=(13,10))
#colorbars are added as their own plot so 3+1=4 rows and 4+1=5 columns. Gridspec is ideal for this kind of thing
gs = GridSpec(4, 5, width_ratios=[0.05,1,1, 1, 1],height_ratios=[0.1,1,1,1])
ax = [[fig.add_subplot(gs[i, j]) for j in range(5)] for i in range(4)]

#don't show the top row and first column subplots.
ax[0][0].remove()
ax[0][1].remove()
ax[0][2].remove()
ax[0][3].remove()
ax[0][4].remove()
ax[1][0].remove()
ax[2][0].remove()
ax[3][0].remove()

#neutral fraction panel
im31 = ax[3][1].imshow(avg_neutral_fraction,cmap=cm.binary,vmax=0.8,vmin=0,origin='lower')
cb31 = plt.colorbar(im31,
    cax=plt.subplot(gs[3,0]),
    orientation='vertical',
    ticklocation = 'left')
cb31.ax.tick_params(direction='out')
cb31.ax.locator_params(nbins=5)
cb31.set_label(r"$\left< x_\mathrm{HI} \right>$")#,x=0.94)

#average gas density panel
im20 = ax[2][1].imshow(avg_Delta,cmap=cm.binary,vmax=3,vmin=0,origin='lower')
cb20 = plt.colorbar(im20,
    cax=plt.subplot(gs[2,0]),
    orientation='vertical',
    ticklocation = 'left')
cb20.ax.tick_params(direction='out')
cb20.ax.locator_params(nbins=5)
cb20.set_label(r"$\left< \Delta \right>$")#,y=0.94)

ticks= np.linspace(0,Ngas,5)
labels = np.asarray(np.linspace(0,40,5),int)
for i in range(1,4): #looping over rows
    for j in range(1,5): #looping over columns
        if ((i==1)&(j==1)): #this will be the SB CDF plot
            pass
        else:
            ax[i][j].set_xticks(ticks=ticks)
            ax[i][j].set_xticklabels(labels)
            ax[i][j].set_yticks(ticks=ticks)
            ax[i][j].set_yticklabels(labels)
            ax[i][j].scatter(x=[Ngas/2],y=[Ngas/2],marker="+",color="red")
            ax[i][j].tick_params(bottom=True, top=True, left=True, right=True)
            ax[i][j].set_xlim(0,Ngas)

start_idx = 15 #starting at 15. Found to look slightly better than starting at 0
for j in range(N_cubes_per_row): #looping over the 3 SB columns
    SB_collExc_intrinsic = np.load(all_paths[0][j])
    SB_collExc_intrinsic = np.sum(SB_collExc_intrinsic[:,start_idx:start_idx+index,:],axis=1)
    SB_collExc_intrinsic = recenter_image(SB_collExc_intrinsic,xx_recenter, yy_recenter,Ngas)

    SB_Recomb_intrinsic = np.load(all_paths[1][j])/sr_to_arcsec2 #need to convert from 1/sr to 1/arcsec2
    SB_Recomb_intrinsic = np.sum(SB_Recomb_intrinsic[:,start_idx:start_idx+index,:],axis=1)
    SB_Recomb_intrinsic = recenter_image(SB_Recomb_intrinsic,xx_recenter, yy_recenter,Ngas)

    SB_combined_afterRT = np.load(all_paths[2][j])/sr_to_arcsec2  #need to convert from 1/sr to 1/arcsec2
    SB_combined_afterRT = np.sum(SB_combined_afterRT[:,:,start_idx:start_idx+index],axis=2)
    SB_combined_afterRT = recenter_image(SB_combined_afterRT,xx_recenter, yy_recenter,Ngas)

    #textboxed for the three models
    #plus 2 because we are ignoring the first two columns
    ax[1][j+2].text(x=0.02,y=0.91,s=r"$\alpha=${}, {}".format(alpha_list[j],reion_hist_list[j]),
        bbox=dict(facecolor='white', alpha=0.8),transform=ax[1][j+2].transAxes,color=color_list[j],fontsize=fontsize-1)

    im = ax[1][j+2].imshow(SB_collExc_intrinsic*multiplier,origin='lower',vmin=0,vmax=vmax,cmap=cm.inferno)
    ax[2][j+2].imshow((SB_collExc_intrinsic+SB_Recomb_intrinsic)*multiplier,origin='lower', vmin=0,vmax=vmax,cmap=cm.inferno)
    ax[3][j+2].imshow(SB_combined_afterRT*multiplier,origin='lower',vmin=0,vmax=vmax,cmap=cm.inferno)
    im_list.append(im) #only for the 1th row because that's from where we'll be calculating the CDF

    SB_mocks_collExc.append(SB_collExc_intrinsic)
    SB_mocks_afterRT_combined.append(SB_combined_afterRT)

    #total SB for each component before RT (for reference and debugging)
    print("total (coll. exc., recomb.)=({:.1e},{:.1e})".format(np.sum(SB_collExc_intrinsic),np.sum(SB_Recomb_intrinsic)))

#top colorbar for the SB maps
cb_top = plt.colorbar(
    im_list[0],
    cax=plt.subplot(gs[0, 2:5]),
    orientation='horizontal', ticklocation = 'top')
cb_top.ax.tick_params(direction='out')
cb_top.ax.locator_params(nbins=5)
cb_top.set_label(r"Surface Brightness [10$^{-21} $erg s$^{-1}$ cm $^{-2}$ arcsec$^{-2}$]",y=0.94)

#CDF's for the intrinsic SB maps
bins = np.logspace(-22,-20,50)
for j in range(N_cubes_per_row):
    SB_mock = SB_mocks_collExc[j]
    counts,_ = np.histogram(SB_mock,bins=bins,density=True)
    cumsum = np.cumsum(counts)/np.sum(counts)
    vline = np.interp(0.99,cumsum,bins[1:])
    ax[1][1].plot(10**np.log10(bins[1:]),cumsum,color=color_list[j])
    ax[1][1].axvline(10**np.log10(vline),color=color_list[j],label="{:.1e}".format(vline))
    print("CDF(0.99)={:.1e}".format(vline))
ax[1][1].legend(loc="lower left", fontsize=10,title='CDF(0.99)',title_fontsize="small")
ax[1][1].set_xscale('log')
ax[1][1].set_xlim(6*10**-23,1*10**-20)
ax[1][1].set_aspect(2)
ax[1][1].set_xlabel("I-front SB before "+r"Ly$\alpha$ RT",labelpad=10)
ax[1][1].xaxis.set_label_position('top')
ax[1][1].set_ylabel("CDF")

ax[3][1].set_xlabel("cMpc/h")
ax[3][2].set_xlabel("cMpc/h")
ax[3][3].set_xlabel("cMpc/h")
ax[3][4].set_xlabel("cMpc/h")

row_labels = ["I-fronts","IF's + recomb.",r"Ly$\alpha$RT on IF's+recomb."]
for i in range(1,4): #labeling each SB row on the right side of figure
    ax[i][4].yaxis.set_label_position('right')
    ax[i][4].tick_params(labelleft=False, labelright=False)
    ax[i][4].set_ylabel(row_labels[i-1],rotation=-90,labelpad=20)

ax[1][1].xaxis.tick_top()
ax[1][1].tick_params(labelleft=True,labelbottom=False,labeltop=True)
[ax[1][j].tick_params(labelleft=False,labelbottom=False) for j in range(2,5)]
[ax[2][j].tick_params(labelleft=False,labelbottom=False) for j in range(1,5)]
[ax[3][j].tick_params(labelleft=False,labelbottom=True) for j in range(1,5)]

# if dpi=200 you get a weird grid structure that you can see if you look closely
fig.savefig(fname="figures/SB_9_panels_orientation{}.png".format(orientation),dpi=300,bbox_inches='tight')
plt.close()
