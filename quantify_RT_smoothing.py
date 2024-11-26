import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from functions import recenter_image
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter


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
ls_list = ["-","dashed"]

path_beforeRT = "data/march/SB_a0.5_f7.3.npy"
path_afterRT = "data/sept2024/out_cube_Nphot600000.npy"

# setting parameters related to NB816 filter
Ngas = 400
nb_depth = 26 #mpc/h # assuming NB816 filter at z=5.7
sim_depth = 40 #mpc/h
index = int(nb_depth/sim_depth*Ngas+0.5) #used to integrate the SB's over the filter from "a" to "a+index" where "a" is the initial pixel.
sr_to_arcsec2 = 4.25e+10
cell_size = 40/Ngas*1000/0.68/(1+5.7)*3.086e+21 #in proper centimeters

xx_recenter, yy_recenter = 8/40,28/40
start_idx = 15
SB_beforeRT = np.load(path_beforeRT)
SB_beforeRT= np.sum(SB_beforeRT[:,start_idx:start_idx+index,:],axis=1)
SB_beforeRT = recenter_image(SB_beforeRT,xx_recenter, yy_recenter,Ngas)
# SB_beforeRT = SB_beforeRT*(SB_beforeRT>2e-21)

SB_afterRT = np.load(path_afterRT)
SB_afterRT = np.sum(SB_afterRT[:,:,start_idx:start_idx+index],axis=2)
SB_afterRT = recenter_image(SB_afterRT,xx_recenter, yy_recenter,Ngas)/sr_to_arcsec2

print(np.sum(SB_beforeRT))
print(np.sum(SB_afterRT))
# print(np.min(SB_afterRT[SB_afterRT>1e-25]))
# SB_afterRT = SB_afterRT*(SB_afterRT>2e-21)

### SMOOTHING TESTS BAYU WILSON OCT. 3, 2024
# sigma = 0
# SB_beforeRT = gaussian_filter(SB_beforeRT, sigma=sigma)
# SB_afterRT = gaussian_filter(SB_afterRT, sigma=sigma)
# size = 5
# SB_beforeRT = median_filter(SB_beforeRT,size=size)
# SB_afterRT = median_filter(SB_afterRT,size=size)


vmax = 8
multiplier=1e21 #multiplied to the SB so we get numbers between 1 & 10
#Making the 9 panel plot
fig = plt.figure(figsize=(8,8))
#colorbars are added as their own plot so 3+1=4 rows and 4+1=5 columns. Gridspec is ideal for this kind of thing
gs = GridSpec(3, 2, height_ratios=[0.05,1,1],width_ratios=[1,1])
ax = [[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(3)]
ax[0][0].remove()
ax[0][1].remove()
# ax[2][0].remove()
# ax[2][1].remove()

# fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(8,4))

im1 = ax[1][0].imshow(SB_beforeRT*multiplier,origin='lower',vmin=0,vmax=vmax,cmap='inferno')
im2 = ax[1][1].imshow(SB_afterRT*multiplier,origin='lower',vmin=0,vmax=vmax,cmap='inferno')

ax[1][0].text(x=0.02,y=0.91,s=r"before Ly$\alpha$ RT",
    bbox=dict(facecolor='white', alpha=0.8),transform=ax[1][0].transAxes,color=color_list[0],fontsize=fontsize-1)
ax[1][1].text(x=0.02,y=0.91,s=r"after Ly$\alpha$ RT",
    bbox=dict(facecolor='white', alpha=0.8),transform=ax[1][1].transAxes,color=color_list[1],fontsize=fontsize-1)
#top colorbar for the SB maps
cb_top = plt.colorbar(im1,
    cax=plt.subplot(gs[0, :]),
    orientation='horizontal', ticklocation = 'top')
cb_top.ax.tick_params(direction='out')
cb_top.ax.locator_params(nbins=5)
cb_top.set_label(r"Surface Brightness [10$^{-21} $erg s$^{-1}$ cm $^{-2}$ arcsec$^{-2}$]",y=0.94)

ticks= np.linspace(0,Ngas,5)
labels = np.asarray(np.linspace(0,40,5),int)
for i in range(2): #looping over columns
    ax[1][i].set_xticks(ticks=ticks)
    ax[1][i].set_xticklabels(labels)
    ax[1][i].set_yticks(ticks=ticks)
    ax[1][i].set_yticklabels(labels)
    ax[1][i].scatter(x=[Ngas/2],y=[Ngas/2],marker="+",color="red")
    ax[1][i].tick_params(bottom=True, top=True, left=True, right=True)
    ax[1][i].set_xlim(0,Ngas)

ax[1][0].set_ylabel("cMpc/h")


#power spectrum
def pk_2d(signal):
    fft2_signal = np.fft.fft2(signal)
    power_spectrum = np.abs(fft2_signal)**2
    power_spectrum /= np.prod(signal.shape)
    power_spectrum_shifted = np.fft.fftshift(power_spectrum)
    return power_spectrum_shifted

center = Ngas // 2
y, x = np.indices((Ngas, Ngas))
k_vals_grid = np.sqrt((x - center)**2 + (y - center)**2).flatten()
pk_beforeRT = pk_2d(SB_beforeRT*multiplier).flatten()
pk_afterRT = pk_2d(SB_afterRT*multiplier).flatten()
kbin_centers = np.arange(1,Ngas//2,1)
kbin_edges = np.arange(0.5,Ngas//2+0.5,1)
kbin_indices = np.digitize(k_vals_grid, kbin_edges)
binned_pk_beforeRT = np.zeros(len(kbin_edges) - 1)
binned_pk_afterRT = np.zeros(len(kbin_edges) - 1)

for i in range(1, len(kbin_edges)):
    in_bin = (kbin_indices == i)  # Boolean array of values in the current bin
    if np.sum(in_bin) > 0:
        binned_pk_beforeRT[i - 1] = np.mean(pk_beforeRT[in_bin])
        binned_pk_afterRT[i - 1] = np.mean(pk_afterRT[in_bin])

# ax3 = fig.add_subplot(gs[2, :])

ax[2][0].loglog(kbin_centers*2*np.pi/40, binned_pk_beforeRT,label=r"before Ly$\alpha$ RT",color=color_list[0],ls=ls_list[0])
ax[2][0].loglog(kbin_centers*2*np.pi/40, binned_pk_afterRT,label=r"after Ly$\alpha$ RT",color=color_list[1],ls=ls_list[1])
# ax[2][0].loglog(40/kbin_centers, binned_pk_beforeRT,label=r"before Ly$\alpha$ RT",color=color_list[0],ls=ls_list[0])
# ax[2][0].loglog(40/kbin_centers, binned_pk_afterRT*1.5,label=r"after Ly$\alpha$ RT",color=color_list[1],ls=ls_list[1])
# ax[2][0].semilogx(40/kbin_centers, binned_pk_beforeRT/binned_pk_afterRT/1.5,label=r"$\frac{\mathrm{before}}{\mathrm{after}}$ RT",color="k")
ax[2][0].set_xlabel("wavenumber [1/(cMpc/h)]")
ax[2][0].set_ylabel("2D power spectrum, P(k)")
ax[2][0].legend(fontsize=fontsize-2)

LyaRT_labels = ["beforeRT","afterRT"]
bins = np.logspace(-22,-20,50)
for i,SB in enumerate([SB_beforeRT,SB_afterRT]):
    counts,_ = np.histogram(SB,bins=bins,density=True)
    cumsum = np.cumsum(counts)/np.sum(counts)
    vline = np.interp(0.99,cumsum,bins[1:])
    vline2 = np.interp(0.90,cumsum,bins[1:])
    print(f"{LyaRT_labels[i]}: 0.99->{vline:.3e}")
    print(f"{LyaRT_labels[i]}: 0.90->{vline2:.3e}")
    print(f"{LyaRT_labels[i]}: avg(SB) = {np.mean(SB):.3e}")
    # print(np.sum(counts))

    ### PDF
    #ax[2][1].semilogx(bins[1:],counts/np.sum(counts),color=color_list[i],ls=ls_list[i],label=LyaRT_labels[i])

    ### CDF
    ax[2][1].plot(10**np.log10(bins[1:]),cumsum,color=color_list[i],ls=ls_list[i],label=LyaRT_labels[i])
    ax[2][1].axvline(10**np.log10(vline),label="{:.3e}".format(vline),color=color_list[i],ls=ls_list[i])
# ax[2][1].set_ylim(0.87,1.0)

"""
#BAYU TEST, 10/3/24
sigma = 0.5  # Standard deviation for the Gaussian kernel
smoothed_data = gaussian_filter(SB_afterRT, sigma=sigma)
counts,_ = np.histogram(smoothed_data,bins=bins,density=True)
ax[2][1].semilogx(bins[1:],counts/np.sum(counts),color='purple',ls='dotted',label=f"sigma={sigma:.1f}")
"""

ax[2][1].legend(loc="upper left", fontsize=10),#title='CDF(0.99)',title_fontsize="small")
ax[2][1].set_xscale('log')
# ax[2][1].set_xlim(6*10**-23,10**-20)
# ax[2][1].set_aspect(2)
ax[2][1].set_xlabel("I-front SB")#,labelpad=10)
# ax[2][1].xaxis.set_label_position('top')
# ax[2][1].set_ylabel("SB PDF")
ax[2][1].set_ylabel("SB CDF")


# ax[2][1].axvline(1.15e-21)

# ax[2][1].loglog(kbin_centers*2*np.pi/40, binned_pk_beforeRT,label="beforeRT")
# ax[2][1].loglog(kbin_centers*2*np.pi/40, binned_pk_afterRT,label="afterRT")
# ax[2][1].set_xlabel("wavenumber [1/(cMpc/h)]")
# ax[2][1].set_ylabel("power spectrum, P(k)")
# ax[2][1].legend()

ax

plt.tight_layout()
plt.savefig("figures/quantify_RT_smoothing_Nov24.png",dpi=350)
plt.close()
# plt.show()
