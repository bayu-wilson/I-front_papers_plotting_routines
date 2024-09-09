#importing libraries
import numpy as np
import matplotlib.pyplot as plt
from functions import make_mock_maps,recenter_image,make_mock_maps2 #, create_circular_tophat_kernel, convolve_with_circular_tophat
from matplotlib import cm
from scipy.ndimage import gaussian_filter
from matplotlib.patches import FancyArrowPatch, ConnectionPatch

#plotting parameters
fontsize=8
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

#three observational scenarios
scope_list = [(8.2,200),(8.2,200),(4.5,5000)]
seed = 12346 #we will initialize each panel with this same seed
kernels = [5,10,25] #in units of pixels
im_list = []

#user-input
Ngas=400 #number of cells in one dimension for the inputted data cube
nb_depth = 26 #mpc/h
sim_depth = 40 #mpc/h
index = int(nb_depth/sim_depth*Ngas+0.5)
sr_to_arcsec2 = 4.25e+10
cell_size = 40/Ngas*1000/0.68/(1+5.7)*3.086e+21

#path for post-LyaRT SB with more rapid, hard alpha model
dataCubeDir = "data/march/"
dataCubePath = dataCubeDir+"out_cube_Nphot600000_hard_moreRapid_wRecomb.npy"
#recombination fields don't include LyaRT because the LyaRT smoothing scale is much less than the smoothing from the convolutions
recomb_file = dataCubeDir+"RR_a1.5_f1.5.npy"

#loading in SB maps
start_idx = 15 #same as in SB_9_panels.py
xx_recenter, yy_recenter = 8/40,28/40 #same as in SB_9_panels.py
SB_Recomb_intrinsic_3D = np.load(recomb_file)/sr_to_arcsec2 #*cell_size
SB_Recomb_intrinsic = np.sum(SB_Recomb_intrinsic_3D[:,start_idx:start_idx+index,:],axis=1)
SB_Recomb_intrinsic = recenter_image(SB_Recomb_intrinsic,xx_recenter, yy_recenter,Ngas)

ncols = len(kernels)*2 #each kernel will have two columns
nrows = len(scope_list)
fig,ax=plt.subplots(ncols=ncols,nrows=nrows,figsize=(8,4))

#formatting for the plot
for i in range(nrows):
    for j in range(ncols):
        ax[i,j].set_xlim(0,Ngas-1)
        ax[i,j].set_ylim(0,Ngas-1)
        ticks = np.linspace(0,Ngas-1,4)
        labels = np.asarray(np.linspace(0,24,4),int)
        ax[i,j].scatter(x=[Ngas/2],y=[Ngas/2],s=25,marker="+",color="red",alpha=0.5)
        ax[i,j].set_xticks(ticks=ticks)
        ax[i,j].set_xticklabels(labels)
        ax[i,j].set_yticks(ticks=ticks)
        ax[i,j].set_yticklabels(labels)
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')
        ax[i,j].tick_params(labelleft=False,labelbottom=False)

### MOCK MAPS AND NOISE+RECOMBS (second and third row)
for i in range(1,nrows): #scopes
    for j in range(0,ncols,2): #kernels
        D_m,texp_hr = scope_list[i]
        np.random.seed(seed)
        noise_only, signal_only = make_mock_maps(dataCubePath=dataCubePath, D_m=D_m,texp_hr=texp_hr,
            Ngas=Ngas,add_noise=True)
        _,recomb_only=make_mock_maps2(I_ifront=SB_Recomb_intrinsic,D_m=D_m,texp_hr=texp_hr,
            Ngas=Ngas,add_noise=False)
        mock_map = signal_only + noise_only
        mock_map = gaussian_filter(mock_map, sigma=kernels[j//2],mode='wrap')  #j//2 because there are double the amount of columns than kernels
        background = recomb_only + noise_only
        background = gaussian_filter(background, sigma=kernels[j//2],mode='wrap')

        vmin,vmax = np.min([np.min(background),np.min(mock_map)]),np.max([np.max(background),np.max(mock_map)])
        ax[i,j].imshow(mock_map, cmap=cm.binary,origin='lower',vmin=vmin,vmax=vmax) #pairs will have the same colorbar!
        ax[i,j+1].imshow(background, cmap=cm.binary,origin='lower',vmin=vmin,vmax=vmax)
        print(f"(D,t) = ({D_m:.1f},{texp_hr:d}), Kernel: {kernels[j//2]:d}, vmin: {int(vmin):d}, vmax: {int(vmax):d}")

### intrinsic maps (first row)
for j in range(0,ncols,2):
    D_m,texp_hr = scope_list[0]
    np.random.seed(seed)
    noise_only, signal_only = make_mock_maps(dataCubePath=dataCubePath, D_m=D_m,texp_hr=texp_hr,
        Ngas=Ngas,add_noise=True)
    intrinsic_map = gaussian_filter(signal_only, sigma=kernels[j//2],mode='wrap')

    ax[0,j].imshow(intrinsic_map, cmap=cm.binary,origin='lower')  #vmin and vmax don't matter. We just want a to visualize the intrinsic signal
    fig.delaxes(ax[0, j+1])

for j in range(0,ncols,2):
    kernel_arcsec = kernels[j//2]*3.63 #pixel scale in arcsec
    x_start = 400  # x-coordinate for the start of the bracket
    y_start = 440    # y-coordinate for the start of the bracket
    x_end = x_start    # x-coordinate for the end of the bracket
    y_end = y_start + 50      # y-coordinate for the end of the bracket

    # Add a curly bracket using annotations
    ax[0,j].annotate('', xy=(x_start, y_start), xytext=(x_end, y_end),
                arrowprops=dict(arrowstyle='-[, widthB=8.4, lengthB=0.4', color='black'), annotation_clip=False)

    # Add labels for the bracket
    ax[0,j].text((x_start + x_end) / 2, y_end + 0.5,  r"$\sigma=$"+f"{kernel_arcsec:.1f}"+r"$''$", ha='center')

    ax[-1,j].set_xlabel("full mock",fontsize=fontsize)
    ax[-1,j+1].set_xlabel("recomb. + noise",fontsize=fontsize)
    ax[1,j+1].tick_params(labeltop=True)

    ax[0,j].tick_params(labelright=True)
    yticks = ax[0,j].get_yticks()
    yticklabels = ax[0,j].get_yticklabels()
    yticklabels[0] = ''
    ax[0,j].set_yticklabels(yticklabels)

ax[0,0].set_ylabel("intrinsic")
ax[1,0].set_ylabel(f"{scope_list[1][0]}m, {scope_list[1][1]}hr")
ax[2,0].set_ylabel(f"{scope_list[2][0]}m, {scope_list[2][1]}hr")

ax[1,-1].yaxis.set_label_position('right')
ax[1,-1].set_ylabel("arcmin",fontsize=fontsize,rotation=-90,labelpad=10)
ax[1,-1].tick_params(labelright=True)
ax[2,-1].tick_params(labelright=True)

fig.savefig(fname="figures/mock18_binned_scopes.png",dpi=400,bbox_inches='tight')
plt.close()
