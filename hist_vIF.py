import numpy as np
import matplotlib.pyplot as plt

fontsize=20
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
# plt.rcParams['text.usetex'] = True

Ngas_FRS = 200 #number of RT cells in Full Reion Sim. (FRS)
Ngas_NIS = 400 #number of RT cells in Neutral Island Sim. (NIS)
vIF_FRS = np.load("data/FullReionSim/vIF_masked_z=5.7.npy") #=

vIF_min = 1e2
vIF_max = 1e5

NIS_files = ["data/march/vIF_a1.5_f1.5.npy","data/march/vIF_a1.5_f7.3.npy"]
color_list_NIS = ["red","orange","green"]
label_list_NIS = [r"SLOWER-v$_\mathrm{IF}$",r"FASTER-v$_\mathrm{IF}$"]

vIF_FRS = vIF_FRS[(vIF_FRS>vIF_min)&(vIF_FRS<vIF_max)]

fig, ax = plt.subplots(figsize=(7,7))
vIF_bin_edges = np.logspace(2,5,50)

### vIF PDF for full reionization simulation
weights_FRS = np.ones_like(vIF_FRS)/float(len(vIF_FRS))
bin_values,_,_ = ax.hist(vIF_FRS,bins=vIF_bin_edges,histtype="step",weights=weights_FRS,color="blue",label="Reference Sim.")
index_max = np.where(bin_values==max(bin_values))[0][0] #finding index where PDF peaks
logvIF_bin_width = np.mean(np.log10(vIF_bin_edges[1:])-np.log10(vIF_bin_edges[:-1]))
max_peak = 10**(np.log10(vIF_bin_edges[:-1][index_max])+logvIF_bin_width/2) #vIF at the center of the vIF bin
print(f"FRS: {max_peak:.3f} km/s") #

### vIF PDF for  neutral island simulations
NIS_norm = [1,1]#[1.05,0.97]
for i,file_i in enumerate(NIS_files):
    vIF_NIS = np.load(file_i).flatten()
    vIF_NIS = vIF_NIS[(vIF_NIS>vIF_min)&(vIF_NIS<vIF_max)]
    weights_NIS = np.ones_like(vIF_NIS)/float(len(vIF_NIS))*NIS_norm[i]#*1.31
    bin_values,_,_ = ax.hist(vIF_NIS,bins=vIF_bin_edges,histtype="step",weights=weights_NIS,
                            color=color_list_NIS[i],label=f"{label_list_NIS[i]}")
    index_max = np.where(bin_values==max(bin_values))[0][0]
    logvIF_bin_width = np.mean(np.log10(vIF_bin_edges[1:])-np.log10(vIF_bin_edges[:-1]))
    max_peak = 10**(np.log10(vIF_bin_edges[:-1][index_max])+logvIF_bin_width/2)
    print(f"The I-front speed mode is at {max_peak:.3f} km/s" )

ax.set_xscale("log")
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True,
                     direction='in')
ax.tick_params(which="minor", labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=True, top=True, left=True, right=True,
                     direction='in')

ax.set_xlim(1e2,1e5)
ax.set_ylabel(r"dP/ln v$_\mathrm{IF}$")
ax.set_xlabel(r"v$_\mathrm{IF}$ [km s$^{-1}$]")
legend = ax.legend(ncols=3,loc='upper center',bbox_to_anchor=(0.5, 1.10),fontsize = 12)
fig.tight_layout()
# plt.show()
plt.savefig("figures/hist_vIF.png",dpi=250)
print("Figure located here: figures/hist_vIF.png")
plt.close()
