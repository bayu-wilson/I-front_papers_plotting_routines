#importing libraries
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.ndimage import gaussian_filter#, median_filter
from functions import make_mock_maps2,recenter_image


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

#columns
recomb_file = "data/march/RR_a1.5_f1.5.npy"
IF_files = ["data/march/out_cube_Nphot600000_soft_lessRapid_wRecomb.npy",
            "data/march/out_cube_Nphot600000_soft_moreRapid_wRecomb.npy"]
model_names   = [r"$\alpha=1.5$, SLOWER-v$_\mathrm{IF}$",
                r"$\alpha=1.5$, FASTER-v$_\mathrm{IF}$"]
apertures = [4.5,4.5]
exposure_times = [5000,5000]
sigma_gauss = [25]

Ngas=400
nb_depth = 26 #mpc/h
sim_depth = 40 #mpc/h
index = int(nb_depth/sim_depth*Ngas+0.5)
index_half = int(index/2)
sr_to_arcsec2 = 4.25e+10
cell_size = 40/Ngas*1000/0.68/(1+5.7)*3.086e+21

SB_Recomb_intrinsic_3D = np.load(recomb_file)/sr_to_arcsec2
SB_Recomb_intrinsic = np.sum(SB_Recomb_intrinsic_3D[:,15:15+index,:],axis=1)
SB_Recomb_intrinsic = recenter_image(SB_Recomb_intrinsic,8/40,28/40,Ngas)

#choosing augmentations of the background to make the wide-field maps. should be same as in panels_large_map.py
A= np.sum(SB_Recomb_intrinsic_3D[:,:index,:],axis=1).T
B= np.sum(SB_Recomb_intrinsic_3D[:index,:,:],axis=0).T
C= np.sum(SB_Recomb_intrinsic_3D[:,:,:index],axis=2)

D= np.sum(SB_Recomb_intrinsic_3D[:index,:,:],axis=0).T
E= np.sum(SB_Recomb_intrinsic_3D[:,:,:index],axis=2).T

F= np.sum(SB_Recomb_intrinsic_3D[:,:,-index:],axis=2)
G= np.sum(SB_Recomb_intrinsic_3D[-index:,:,:],axis=0)
H= np.sum(SB_Recomb_intrinsic_3D[:,-index:,:],axis=1)
H = np.flipud(H)

total_recombs = np.block([
    [A,B,C],
    [D,SB_Recomb_intrinsic,E],
    [F,G,H]])

### mask to only choose I-front region. I-front cut-out
y,x = np.ogrid[0:400,0:400]

x1,y1=(100,0)
x2,y2=(400,180)
m1 = (y2-y1)/(x2-x1)
mask= y>m1*(x-x1) #want above this line

mask *= x>190

x1,y1=(136,244)
x2,y2=(344,180)
m1 = (y2-y1)/(x2-x1)
mask *= y<m1*(x-x1)+y1 #want below this line

radius=45
mask2 = (x-145)**2 + (y-220)**2 <= radius**2
mask3 = ((x-190)**2 + (y-174)**2 <= radius**2)
mask4 = ((x-200)**2 + (y-90)**2 <= radius**2)
mask5 = ((x-330)**2 + (y-160)**2 <= radius**2)
mask = mask|mask2|mask3|mask4|mask5
### end of I-front cut-out

nrows = len(sigma_gauss)
ncols = len(IF_files)+1
fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(7,4.5),squeeze=False)
ax = ax.reshape(nrows, ncols)

# vranges = []
im_list = [ [[] for _ in range(ncols)] for _ in range(nrows)]
seed = 12346
# im_noises = []
for i in range(nrows): #there's only one row!
    np.random.seed(seed)
    N_noise,N_recomb=make_mock_maps2(I_ifront=total_recombs,D_m=4.5,texp_hr=5000,
                        Ngas=Ngas*3,add_noise=True)

    mock_obs = N_noise + N_recomb
    mock_obs = gaussian_filter(mock_obs, sigma=sigma_gauss[i],mode='wrap')

    im = ax[i,0].imshow(mock_obs,origin='lower',cmap=cm.binary)
    im_list[i][0] = im

for i in range(nrows): #onle one row!
    for j in range(1,ncols): # obs setup
        SB_collExc_intrinsic_3D = np.load(IF_files[j-1])/sr_to_arcsec2
        SB_collExc_intrinsic = np.sum(SB_collExc_intrinsic_3D[:,:,15:15+index],axis=2)
        SB_collExc_intrinsic = recenter_image(SB_collExc_intrinsic,8/40,28/40,Ngas)#*multiplier
        total_both = np.copy(total_recombs)
        total_both[Ngas:2*Ngas,Ngas:2*Ngas] = SB_Recomb_intrinsic + SB_collExc_intrinsic*mask
        np.random.seed(seed)
        N_noise,N_ifront=make_mock_maps2(I_ifront=total_both,D_m=apertures[j-1],texp_hr=exposure_times[j-1],
                    Ngas=Ngas*3,add_noise=True)
        mock_obs = N_ifront + N_noise
        mock_obs = gaussian_filter(mock_obs, sigma=sigma_gauss[i],mode='wrap')
        im = ax[i,j].imshow(mock_obs,origin='lower',cmap=cm.binary)#,vmin=vmin,vmax=vmax)
        im_list[i][j] = im

#setting color scale range
for i in range(nrows):
    for j in range(ncols):
        #I got these numbers from the values printed out from mock18_binned_scopes!
        im_list[i][j].set_clim(vmin=44, vmax=11582)

for j in range(1,ncols):
    ax[0,j].set_title(f"{model_names[j-1]}",fontsize=fontsize-3)

#formatting
pix_sz_arcsec = 3.7392
pix_arr = np.linspace(0,Ngas*3-1,Ngas*3)
theta_arcmin = pix_sz_arcsec*Ngas*3 / 60
theta_arr = np.linspace(0,theta_arcmin,Ngas*3)
labels = [ 0,  24, 48, 72] #arcmin labels
ticks = [int(np.interp(i,theta_arr,pix_arr)+0.5) for i in labels]

for i in range(nrows):
    sigma_gauss_arcsec = sigma_gauss[i]*pix_sz_arcsec #pixel scale in arcsec
    text1 = r"$\sigma_g=$"+f"{sigma_gauss_arcsec:.1f}"+r"$''$"
    ax[i,-1].set_ylabel(text1,rotation=-90,labelpad=15,fontsize=fontsize-3)
    ax[i,-1].yaxis.set_label_position('right')

    for j in range(ncols):
        ax[i,j].set_xlim(0,Ngas*3-1)
        ax[i,j].set_ylim(0,Ngas*3-1)

        # ticks = np.linspace(0,Ngas*3-1,4)
        # labels = np.asarray(np.linspace(0,24*3,4),int)
        ax[i,j].set_xticks(ticks=ticks)
        ax[i,j].set_xticklabels(labels)
        ax[i,j].set_yticks(ticks=ticks)
        ax[i,j].set_yticklabels(labels)
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')
        ax[i,j].xaxis.set_ticks_position('both')
        ax[i,j].yaxis.set_ticks_position('both')

ax[0,0].set_ylabel("arcmin",fontsize=fontsize-2)
ax[0,0].tick_params(labelleft=True)
ax[0,0].set_title(f"recomb. + noise",fontsize=fontsize-3)

for i in range(nrows): #kernels
    for j in range(ncols):
        ax[i,j].tick_params(labelbottom=False,labelleft=False)

for i in range(nrows):
    ax[i,0].tick_params(labelleft=True)
for j in range(ncols):
    ax[-1,j].tick_params(labelbottom=True)

fig.savefig(fname="figures/panels_large_map_other_models.png",dpi=400,bbox_inches='tight')
plt.close()
