#importing libraries
import numpy as np
import matplotlib.pyplot as plt
from functions import make_mock_maps,recenter_image, create_circular_tophat_kernel, convolve_with_circular_tophat
from matplotlib import cm

dataCubeDir = "data/march/"
dataCubePath = dataCubeDir+"out_cube_Nphot600000_hard_moreRapid_wRecomb.npy"

#plotting parameters
fontsize=12
lw=1.5
major_tick_size = 5
minor_tick_size = 3
tick_dir = "inout"
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir

#user-input
Ngas=400 #number of cells in one dimension for the inputted data cube

#three observational scenarios
scope_list=[(8.2,200),(4.5,5000)]
seed = 12346
kernels = [5,10,25]
im_list = []

nrows = len(scope_list)
ncols = len(kernels)

fig,ax=plt.subplots(nrows=nrows,ncols=ncols,figsize=(6,4))
for i in range(nrows): #different kernels (2)
    for j in range(ncols): #number of observations being tested (3)
        D_m,texp_hr = scope_list[i]
        np.random.seed(seed)
        noise_only, signal_only = make_mock_maps(dataCubePath=dataCubePath, D_m=D_m,texp_hr=texp_hr,
            Ngas=Ngas,add_noise=True)


        radius = int(np.sqrt(2*np.log(2))*kernels[j]+0.5) #FWHM/2
        FWHM_arcsec = radius*2*3.63
        ax[-1,j].set_xlabel(r"D$_\mathrm{tophat}=$"+f"{FWHM_arcsec:.1f}"+r"$''$")#,rotation=-90,labelpad=15)

        signal_only,pixels_in_kernel = convolve_with_circular_tophat(signal_only,radius=radius)
        noise_per_pixel = np.sqrt(pixels_in_kernel)*np.std(noise_only)

        #plotting S/N maps
        im =  ax[i,j].imshow(signal_only/noise_per_pixel,origin='lower',cmap=cm.nipy_spectral,vmin=0,vmax=16)
        im_list.append(im)

        #formatting
        ax[i,j].set_xlim(0,Ngas-1)
        ax[i,j].set_ylim(0,Ngas-1)
        ticks = np.linspace(0,Ngas-1,4)
        labels = np.asarray(np.linspace(0,24,4),int)
        ax[i,j].set_xticks(ticks=ticks)
        ax[i,j].set_xticklabels(labels)
        ax[i,j].set_yticks(ticks=ticks)
        ax[i,j].set_yticklabels(labels)
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')
        if j==0:
            ax[i,j].tick_params(labelbottom=False,labelleft=False)
        elif j==1:
            ax[i,j].tick_params(labelbottom=False, labelleft=False)
        elif j==2:
            ax[i,j].set_ylabel("arcmin",rotation=-90,labelpad=15)
            ax[i,j].yaxis.set_label_position('right')
            ax[i,j].tick_params(labelbottom=False, labelleft=False,labelright=True)

for j in range(ncols):
    ax[-1,j].tick_params(labelbottom=True)

ax[0,0].set_ylabel(f"{scope_list[0][0]}m, {scope_list[0][1]}hr",fontsize=fontsize+1)
ax[1,0].set_ylabel(f"{scope_list[1][0]}m, {scope_list[1][1]}hr",fontsize=fontsize+1)

fig.subplots_adjust(top=0.8)  # Make space at the top
cbar_ax = fig.add_axes([0.15, 0.85, 0.7, 0.02])  # [left, bottom, width, height]

# Add colorbar
cb = fig.colorbar(im_list[1], cax=cbar_ax, orientation='horizontal')
cb.ax.tick_params(direction='out')
cb.ax.locator_params(nbins=10)
cb.set_label(r"Signal-to-noise")#,y=0.94)

# Move colorbar label and tick labels to the top
cbar_ax.xaxis.set_ticks_position('top')
cbar_ax.xaxis.set_label_position('top')

fig.savefig(fname="figures/SNR_binned_scopes_circularTopHat.png",dpi=400,bbox_inches='tight')
plt.close()
