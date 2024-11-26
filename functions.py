import numpy as np
import scipy.ndimage as ndimage

def recenter_image(X,xx,yy,Ngas):
    """
    Used to translate the image by xx*Ngas pixels rightward and yy*Ngas pixels upward with periodic boundary conditions.
    Note that the origin in the maps are usually set to the bottom left corner
    """
    X_shift = np.zeros_like(X)
    recenter = [int(yy*Ngas),int(xx*Ngas)] #30/40,7/40
    for ii in range(Ngas): #y
        for jj in range(Ngas): #x
            X_shift[ii,jj] = X[ii-recenter[0],jj-recenter[1]] #[y,x]
    return X_shift

def make_mock_maps(dataCubePath,D_m,texp_hr,Ngas,add_noise):
	"""
    Produces mock maps given observational parameters and SB data cube path.
    dataCubePath (str): path to the data cube.
	D_m (int): aperture diameter in meters
	texp_hr (float): Exposure time in hours
	Ngas (int): Ngas^3 is the total number of cells in the datacube

    Returns:
    N_noise: counts from poisson noise per pixel in shape [Ngas,Ngas]
    N_ifront: counts from the signal per pixel in shape [Ngas,Ngas]
	"""
	# Constants
	c              = 2.998e10
	h              = 6.63e-27
	cm2_per_m2     = 100**2
	Ang_per_micron = 1e4
	sr_to_arcsec2 = 4.25e+10

	# Filter characteristics
	lambda_Ang     = 8160.
	lambda_cm      = lambda_Ang * 1e-8
	dlambda_Ang    = 100   # filter width (Ang)

	# Telescope characteristics
	throughput     = 0.2   # This is a rough guess
	A_cm2          = np.pi * ((D_m/2)**2) * cm2_per_m2

	# Integration time
	texp_sec       = texp_hr * 3600

	# Sky brightness
	# Sky radiance (photon/s/m2/micron/arcsec2)
	# 300 = Cerro Paranal new moon, airmass=1 near 8160 Ang
	if add_noise:
		sky_radiance   = 300#/3
	else:
		sky_radiance = 0
	# Photon energy
	photon_erg     = h * c / lambda_cm
	# Specific sky intensity (erg/s/cm2/Angstrom/arcsec2)
	I_lam_sky      = sky_radiance * photon_erg / cm2_per_m2 / Ang_per_micron
	# Conversion from f_lam to f_nu
	def flam2fnu(flam,wave_Ang):
	    return flam*(1e-8)*(wave_Ang**2)/c
	I_nu_sky       = flam2fnu(I_lam_sky,lambda_Ang)

	HSC_dtheta_arcsec = 0.168     # Solid angle
	HSC_Omega_arcsec2 = HSC_dtheta_arcsec**2

	# Total number of sky photons collected
	HSC_Ndot_sky = throughput * (I_lam_sky/photon_erg) * dlambda_Ang * A_cm2 * HSC_Omega_arcsec2

	# Observed HSC NB0816: Ndot_sky ~ 3 counts/s/pixel = 140 counts/s/arcsec2

	# Read in SB cube (erg/s/cm2/arcsec2)
	SB_cube = np.load(dataCubePath)/sr_to_arcsec2
	sz      = np.shape(SB_cube)
	xsz     = sz[0]
	ysz     = sz[1]

	# Pixels size
	pix_sz_hinvckpc = 100
	h               = 0.68
	z               = 5.7
	pkpc_per_arcsec = 5.87
	pix_sz_ckpc     = pix_sz_hinvckpc / h
	pix_sz_pkpc     = pix_sz_ckpc / (1+z)
	pix_sz_arcsec   = pix_sz_pkpc / pkpc_per_arcsec   
	pix_sz_arcsec2  = pix_sz_arcsec**2

	# Comoving length corresponding to filter
	nb_depth = 26. * (dlambda_Ang/100) #* 40/26

	# Integrate I-front intensity over filter width
	sim_depth = 40 #mpc/h
	index = int(nb_depth/sim_depth*Ngas+0.5)
	I_ifront = np.sum(SB_cube[:,:,15:15+index],axis=2)
	I_ifront = recenter_image(I_ifront,8/40,28/40,Ngas)

	# Total number of I-front photons collected per pixel
	Ndot_ifront = throughput * (I_ifront / photon_erg) * A_cm2 * pix_sz_arcsec2
	N_ifront    = Ndot_ifront * texp_sec

	# Sky photons collected per pixel
	Ndot_sky    = throughput * (I_lam_sky/photon_erg) * dlambda_Ang * A_cm2 * pix_sz_arcsec2
	N_sky       = Ndot_sky * texp_sec

	N_noise = np.random.poisson(lam=N_sky,size=(Ngas,Ngas)) - N_sky #BAYU MARCH 11, 2024
	# observed_image = N_ifront + N_noise
	return N_noise,N_ifront

def make_mock_maps2(I_ifront,D_m,texp_hr,Ngas,add_noise):
	"""
    This version is same as make_mock_maps except that it does not load in the 3D datafile. It takes 2D ``I_ifront'' directly.
    I_ifront (float): signal intensity. dims are [Ngas,Ngas]
	D_m (int): aperture diameter in meters
	texp_hr (float): Exposure time in hours
	Ngas (int): Ngas^3 is the total number of cells in the datacube

    Returns:
    N_noise: counts from poisson noise per pixel in shape [Ngas,Ngas]
    N_ifront: counts from the signal per pixel in shape [Ngas,Ngas]
	"""
	# Constants
	c              = 2.998e10
	h              = 6.63e-27
	cm2_per_m2     = 100**2
	Ang_per_micron = 1e4
	sr_to_arcsec2 = 4.25e+10

	# Filter characteristics
	lambda_Ang     = 8160.
	lambda_cm      = lambda_Ang * 1e-8
	dlambda_Ang    = 100   # filter width (Ang)

	# Telescope characteristics
	throughput     = 0.2   # This is a rough guess
	A_cm2          = np.pi * ((D_m/2)**2) * cm2_per_m2

	# Integration time
	#texp_hr        = texp_hr #100 hours
	texp_sec       = texp_hr * 3600

	# Sky brightness
	# Sky radiance (photon/s/m2/micron/arcsec2)
	# 300 = Cerro Paranal new moon, airmass=1 near 8160 Ang
	if add_noise:
		sky_radiance   = 300#/3
	else:
		sky_radiance = 0
	# Photon energy
	photon_erg     = h * c / lambda_cm
	# Specific sky intensity (erg/s/cm2/Angstrom/arcsec2)
	I_lam_sky      = sky_radiance * photon_erg / cm2_per_m2 / Ang_per_micron
	# Conversion from f_lam to f_nu
	def flam2fnu(flam,wave_Ang):
	    return flam*(1e-8)*(wave_Ang**2)/c
	I_nu_sky       = flam2fnu(I_lam_sky,lambda_Ang)

	HSC_dtheta_arcsec = 0.168     # Solid angle
	HSC_Omega_arcsec2 = HSC_dtheta_arcsec**2

	# Total number of sky photons collected
	HSC_Ndot_sky = throughput * (I_lam_sky/photon_erg) * dlambda_Ang * A_cm2 * HSC_Omega_arcsec2
	# Observed HSC NB0816: Ndot_sky ~ 3 counts/s/pixel = 140 counts/s/arcsec2

	# Pixels size
	pix_sz_hinvckpc = 100 #40 / 400 *1000 = 100
	h               = 0.68
	z               = 5.7
	pkpc_per_arcsec = 5.87
	pix_sz_ckpc     = pix_sz_hinvckpc / h
	pix_sz_pkpc     = pix_sz_ckpc / (1+z)
	pix_sz_arcsec   = pix_sz_pkpc / pkpc_per_arcsec   # 3.73
	pix_sz_arcsec2  = pix_sz_arcsec**2

	# Total number of I-front photons collected per pixel
	Ndot_ifront = throughput * (I_ifront / photon_erg) * A_cm2 * pix_sz_arcsec2
	N_ifront    = Ndot_ifront * texp_sec

	# Sky photons collected per pixel
	Ndot_sky    = throughput * (I_lam_sky/photon_erg) * dlambda_Ang * A_cm2 * pix_sz_arcsec2
	N_sky       = Ndot_sky * texp_sec

    # calculating noise
	N_noise = np.random.poisson(lam=N_sky,size=(Ngas,Ngas)) - N_sky #BAYU MARCH 11, 2024
	# observed_image = N_ifront + N_noise
	return N_noise, N_ifront

def create_circular_tophat_kernel(radius):
    """Create a circular tophat kernel with a given radius."""
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    kernel =  np.asarray(mask,int)
    return kernel

def convolve_with_circular_tophat(image, radius):
    """Convolve a 2D image with a circular tophat kernel of given radius."""
    kernel = create_circular_tophat_kernel(radius)
    pixels_in_kernel = np.sum(kernel)
    convolved_image = ndimage.convolve(image, kernel, mode='wrap')
    return convolved_image,pixels_in_kernel

# no longer being used but could be useful function.
# def rebin_matrix(original_matrix, M):
#     N = original_matrix.shape[0] #number of cells in one direction
#     SF = N // M #how many times does M go into N. SF stands for scaling factor
#
#     rebinned_matrix = np.zeros((M, M))
#
#     for i in range(M):
#         for j in range(M):
#             chunk_sum = np.sum(original_matrix[i*SF:(i+1)*SF, j*SF:(j+1)*SF])
#             rebinned_matrix[i, j] = chunk_sum
#     return rebinned_matrix
