import numpy as np
from astroquery.sdss import SDSS

'''
# these are utilities that were originally written for the purpose of finding
# strange variability, but can likely be used for other things
'''

'''
# Finds the minimum contiguous sum of an array
# arr: The input array
# n  : the size of the contiguous sum (how many elements of arr in the sum)
'''
def getMinPartition(arr, n):
	mat = np.tile(arr, (n, 1))
	
	for index,row in enumerate(mat):
		mat[index]         = np.roll(mat[index], index)
		mat[index][:index] = np.nan
#	print(mat)
	return np.nanargmin(mat.sum(axis=0)) - n + 1

'''
# Downloads a spectrum file from the SDSS database using Astroquery
# plate: the spectrum's plate
# mjd  : the spectrum's modified julian date, the day it was taken
# fiber: the spectrum's fiber
'''
def getSpectrum(plate, mjd, fiber):
	if int(plate) == 10658:
		return None, None, None

	print(str(plate) + "-"+str(mjd)+"-"+str(fiber))
	spectrumFits = SDSS.get_spectra(plate=int(plate), mjd=int(mjd), fiberID=int(fiber))
	# astroquery is kinda dumb, so it sends us a "list" of one element (hence
	# the [0]. The data we care about is in [1] of the hdu, .data specifies we're
	# not looking at the header. Flux is kind of self-expanatory, but loglam is
	# the wavelength in angstroms log 10 (don't ask why, I don't get it)
	# oh btw 'error' is measured in variance. SDSS documentation gives a formula
	# for converting that to a typical error measurement, TODO: attach link to
	# that
	flux = spectrumFits[0][1].data["flux"]
	wavelength = np.power(10, spectrumFits[0][1].data["loglam"])
	ivar = spectrumFits[0][1].data["ivar"]

	# ivar values of 0 are invalid
	validPixels = np.where(ivar>0)

	error = np.reciprocal(np.sqrt(ivar[validPixels]))
	return flux[validPixels], wavelength[validPixels], error

'''
# Finds the intersections in wavelength between two observations
# often, repeated observations of an object cover different wavelengths. This
# function matches them together
'''
def intersectSpectra(spec1flux, spec1wl, spec1error, spec2flux, spec2wl, spec2error):

	# TODO: Maybe round the wavelengths to make this more flexible

	# basically, spec1 and spec2 might have different wavelengths read. So we
	# have to do this black magic to make sure the spectra actually 'line up'
	spec1indices = np.in1d(spec1wl, spec2wl) #apparently in1d is not
	spec2indices = np.in1d(spec2wl, spec1wl) #commutative
	spec1flux = spec1flux[spec1indices]
	spec1wl   = spec1wl  [spec1indices]
	spec1error = spec1error[spec1indices]
	spec2flux = spec2flux[spec2indices]
	spec2wl   = spec2wl  [spec2indices]
	spec2error = spec2error[spec2indices]
	return spec1flux, spec1wl, spec1error, spec2flux, spec2wl, spec2error

