## Official plotting routines for my 2024 paper titled "[_Imaging reionization’s last phases with I-front Lyman-α emissions_](https://arxiv.org/abs/2406.14625v1)" 
Written and documented by Bayu Wilson. Most of the corresponding figures in `figures/` have the same name as the plotting routine.


| **Plotting Routine**    | **Description**                                                                                         |
|------------------|---------------------------------------------------------------------------------------------------------|
| `hist_vIF.py`     | Comparing I-front speed PDF's between full reion. and neutral island simulations   |
| `gas_slice_shaded.py`     |  A slice through a gas density field with a shaded neutral region |
| `zoom_slice.py`| A zoom-in on a neutral region of a slice. Other panels show I-front speed, LyC Flux, and Lya emitted from the I-front|
| `SB_9_panels.py`| 9 panels of Lya SB plots and then another column showing the SB CDF, average density, and average neutral fraction|
| `mock18_binned_scopes.py`| Mock narrowband images of neutral islands. Compare fields with neutral islands to fields of only background |
| `SNR_binned_scopes.py`| SNR maps smoothed with circular aperature. The radius is set by the relationship between FWHM and variance|
| `panels_large_map.py`| The effect of a wider (by 3 times) field of view for neutral island visibility using most optimistic model|
| `panels_large_map_other_models.py`|  Similar as previous but using the other two less optimistic models |
| `quantify_RT_smoothing.py`| This routine was added due to the referee report regarding a requrest for a better quantification of the RT smoothing scale. Note that the corresponding figure does NOT appear in the paper as of Nov. 2024.|
| `recomb_only.py`| Surface brightness maps with only recombinations (no I-front emission) |

| **Other**    | **Description**                                                                                         |
|------------------|---------------------------------------------------------------------------------------------------------|
| `functions.py`     | Helper functions for the plotting routines  |

For more info on the version of the FlexRT code that I used to get these simulations, see [this](https://github.com/bayu-wilson/FlexRT_imaging_reion/tree/main) repository.
