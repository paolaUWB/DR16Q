import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import scipy.constants as sc
from scipy.optimize import curve_fit

MIN_TROUGH_CUTOFF = 0.0

# pulled from http://astronomy.nmsu.edu/drewski/tableofemissionlines.html
# all in angstroms
CIIRest    = 1335.120
CIVRest    = 1549.0524
SiIVRest   = 1396.747
OIRest     = 1302.168

'''
# draws a rectangular patch on a graph. May be deprecated
'''
def displayPatch(figureAxis, x, y, start, length):
	wlSlice = x[start:start+length]
	fluxSlice= y[start:start+length]
	figureAxis.add_patch(Rectangle((wlSlice[0],min(fluxSlice)),wlSlice[-1]-wlSlice[0],max(fluxSlice)-min(fluxSlice),linewidth=1,edgecolor='r',facecolor='gray', alpha=0.5))

'''
# Assuming that this absorption line is SiIV, this function shows where the 
# CIV, CII, and OI lines would be using axvspan. With visual inspection this can be used to
# determine if an absorption line is SiIV or CIV

# figureAxis: The axis that is being modified
# SiIVMin: The beginning of the absorption line, in angstroms
# SiIVMax: The end of the absorption line, in angstroms

# NOTE: SiIVMin < SiIVMax is assumed. undefined behavior otherwise
'''
def plotPotentialLines(figureAxis, SiIVMin, SiIVMax):
	SiIVMinZ = (SiIVMin/SiIVRest) - 1.
	SiIVMaxZ = (SiIVMax/SiIVRest) - 1.

	CIVMin = (SiIVMinZ + 1.) * CIVRest 
	CIVMax = (SiIVMaxZ + 1.) * CIVRest

	CIIMin = (SiIVMinZ + 1) * CIIRest
	CIIMax = (SiIVMaxZ + 1) * CIIRest

	OIMin = (SiIVMinZ + 1) * OIRest
	OIMax = (SiIVMaxZ + 1) * OIRest

	figureAxis.axvspan( CIVMin, CIVMax, color='grey', alpha=0.2)
	figureAxis.axvspan( CIIMin, CIIMax, color='blue', alpha=0.2)
	figureAxis.axvspan( OIMin, OIMax, color='yellow', alpha=0.2)

def powerlaw(wavelength, b, c):
	""" Calculates the power law. 

	Parameters
	----------
	wavelength: array
	Comes from RangesData().    
	b: int
	Initial parameter of powerlaw. 
	c: float
	Initial parameter of powerlaw.

	Returns
	-------
	array
		Power law value in the form of an array.
	"""
	return b * (np.power(wavelength, c))



def findParams(flux, wl, z):
	RLF1 = (1250, 1350)
	RLF2 = (1700, 1800)
	RLF3 = (1950, 2200)

	indices1 = np.logical_and(wl > RLF1[0]*(z+1), wl < RLF1[1]*(z+1))
	indices2 = np.logical_and(wl > RLF2[0]*(z+1), wl < RLF2[1]*(z+1))
	indices3 = np.logical_and(wl > RLF3[0]*(z+1), wl < RLF3[1]*(z+1))
	
	finalIndices = np.logical_or(np.logical_or(indices1, indices2), indices3)
	pars, covar = curve_fit(powerlaw, wl[finalIndices], flux[finalIndices], p0=[1250, -0.5], maxfev=10000)
	return pars



# Generates a graph from the provided information
# f1, wl1, e1: The flux, wavelength, and error values from the first observation
# f2, wl2, e2: The flux, wavelength, and error values from the second observation
# troughStart: The beginning of the trough, measured in angstroms
# bottomWavelength: Where the trough bottom is
# z: The redshift of the quasar itself. This primarily comes from the expansion
#    of the universe
def graphMain(f1, wl1, e1, f2, wl2, e2, troughStart, troughEnd, bottomWavelength, z, spec1, spec2):
	outputFigure = plt.figure(1)
	outputAxis   = outputFigure.add_subplot()
	outputAxis.grid(visible=True)

	# calculate the wavelength range
	startVelocity = 0 / sc.speed_of_light * (10**3)
	endVelocity   = 70000 / sc.speed_of_light* (10**3)
	endWL = np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + z, 2) * (1 - startVelocity) / (1 + startVelocity)) )
	startWL   = np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + z, 2) * (1 - endVelocity) / (1 + endVelocity)) )

	# It's possible that the two observations covered different wavelengths and
	# have different indices, so we use 2 ranges
	goodRange1 = np.where((wl1 > startWL) & (wl1 < endWL))
	goodRange2 = np.where((wl2 > startWL) & (wl2 < endWL))


	plotPotentialLines(outputAxis, troughStart, troughEnd)
	outputAxis.set_xlim(left = (wl1[goodRange1])[0], right = (wl1[goodRange1])[-1] )

	# add red line for arbitrary cutoff
	outputAxis.plot(wl1[goodRange1], np.ones(len(wl1[goodRange1])) * MIN_TROUGH_CUTOFF, color='red', linestyle='dashed')



	outputAxis.plot(wl1[goodRange1], f1[goodRange1], color='red', label=str(spec1.plate)+"-"+str(spec1.mjd)+"-"+str(spec1.fiber))

	outputAxis.plot(wl2[goodRange2], f2[goodRange2], color='navy', label=str(spec2.plate)+"-"+str(spec2.mjd)+"-"+str(spec2.fiber))

	p1 = findParams(f1, wl1, z)
	outputAxis.plot(wl1[goodRange1], powerlaw(wl1[goodRange1], *p1), color='xkcd:light orange', label='curve fit of ' + str(spec1.plate)+"-"+str(spec1.mjd)+"-"+str(spec1.fiber))

	p2 = findParams(f2, wl2, z)
	outputAxis.plot(wl2[goodRange1], powerlaw(wl2[goodRange1], *p2), color='xkcd:azure', label='curve fit of ' + str(spec2.plate)+"-"+str(spec2.mjd)+"-"+str(spec2.fiber))



	plt.axvline(x=troughStart)
	plt.axvline(x=troughEnd)

	worstError = np.max(np.row_stack((e1, e2)), axis=0)
	outputAxis.plot(wl1, worstError, color='black')

	return outputFigure
