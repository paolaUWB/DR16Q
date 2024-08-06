import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from matplotlib.backends.backend_pdf import PdfPages
import sys
import re
import os

from matplotlib.patches import Rectangle
from astroquery.sdss import SDSS
from utility_functions  import read_spectra
from collections import defaultdict, namedtuple
from ASTROsoft import DER_SNR
import time
import scipy.constants as sc


MIN_TROUGH_CUTOFF = 0.0

# pulled from http://astronomy.nmsu.edu/drewski/tableofemissionlines.html
# all in angstroms
CIIRest    = 1335.120
CIVRest    = 1549.0524
SiIVRest   = 1396.747
OIRest     = 1302.168


start_time = time.time()

Spectra = namedtuple("Spectra", ["plate", "mjd", "fiber"])

pp = PdfPages('test.pdf')

def toDered(inputSpec):
	return f"spec-{inputSpec.plate}-{inputSpec.mjd}-{inputSpec.fiber}-dered.dr16"

def extractIntegerArrayFromCSV(inputString):
	tmp = inputString.strip('[]').split()
	return list(map(int, tmp))

def extractSpectraTupleFromString(inputString):
	nums = re.findall('\d+', inputString)
	return Spectra(plate=nums[0], mjd=nums[1], fiber=nums[2])

def getMinPartition(arr, n):
	mat = np.tile(arr, (n, 1))
	
	for index,row in enumerate(mat):
		mat[index]         = np.roll(mat[index], index)
		mat[index][:index] = np.nan
#	print(mat)
	return np.nanargmin(mat.sum(axis=0)) - n + 1
	#return np.argmin(arr)


def displayPatch(figureAxis, x, y, start, length):
	wlSlice = x[start:start+length]
	fluxSlice= y[start:start+length]
	figureAxis.add_patch(Rectangle((wlSlice[0],min(fluxSlice)),wlSlice[-1]-wlSlice[0],max(fluxSlice)-min(fluxSlice),linewidth=1,edgecolor='r',facecolor='gray', alpha=0.5))

def getSpectrum(plate, mjd, fiber):
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
	ivar = spectrumFits[0][1].data["flux"]

	# ivar values of 0 are invalid
	validPixels = np.where(ivar>0)

	error = np.reciprocal(np.sqrt(ivar[validPixels]))
	return flux[validPixels], wavelength[validPixels], error

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

def plotPotentialLines(figureAxis, SiIVMin, SiIVMax, hubble_redshift):
	SiIVMinZ = np.power( (1 + hubble_redshift)/(SiIVMin/SiIVRest), 2)
	SiIVMinZ = (SiIVMin/SiIVRest) - 1.#(SiIVMinZ - 1)/(SiIVMinZ + 1)

	SiIVMaxZ = np.power( (1 + hubble_redshift)/(SiIVMax/SiIVRest), 2)
	SiIVMaxZ = (SiIVMax/SiIVRest) - 1.# (SiIVMaxZ - 1)/(SiIVMaxZ + 1)

	CIVMin = (SiIVMinZ + 1.) * CIVRest #np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMinZ) / (1 + SiIVMinZ)) ) 
	CIVMax = (SiIVMaxZ + 1.) * CIVRest #np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMaxZ) / (1 + SiIVMaxZ)) )

	CIIMin = (SiIVMinZ + 1) * CIIRest#np.sqrt( ( np.power(CIIRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMinZ) / (1 + SiIVMinZ)) ) 
	CIIMax = (SiIVMaxZ + 1) * CIIRest#np.sqrt( ( np.power(CIIRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMaxZ) / (1 + SiIVMaxZ)) )

	OIMin = (SiIVMinZ + 1) * OIRest#np.sqrt( ( np.power(OIRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMinZ) / (1 + SiIVMinZ)) ) 
	OIMax = (SiIVMaxZ + 1) * OIRest#np.sqrt( ( np.power(OIRest, 2) * np.power(1 + hubble_redshift, 2) * (1 - SiIVMaxZ) / (1 + SiIVMaxZ)) )

	print("SiIVmin: " + str(SiIVMin))
	print("SiIVmax: " + str(SiIVMax))


	print("civmin: " + str(CIVMin))
	print("civmax: " + str(CIVMax))

	figureAxis.axvspan( CIVMin, CIVMax, color='grey', alpha=0.2)
	figureAxis.axvspan( CIIMin, CIIMax, color='blue', alpha=0.2)
	figureAxis.axvspan( OIMin, OIMax, color='yellow', alpha=0.2)

# assume indices correspond to spec1
# returns:
# (isWeird, (I_01, I_02, I_1, I_2), figure)
#  boolean  float float float float object
def findStats(spec1, spec2, minWLindex, maxWLindex, objectName, indexNum, z):
#	plt.clf()
	print("oomf")
	spec1flux, spec1wl, spec1error = getSpectrum(spec1.plate, spec1.mjd, spec1.fiber)
	spec2flux, spec2wl, spec2error = getSpectrum(spec2.plate, spec2.mjd, spec2.fiber)

	#velocities = wavelength_to_velocity(z, spec1wl)
	startVelocity = 0 / sc.speed_of_light * (10**3)
	endVelocity   = 70000 / sc.speed_of_light* (10**3)

	endWL = np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + z, 2) * (1 - startVelocity) / (1 + startVelocity)) )
	startWL   = np.sqrt( ( np.power(CIVRest, 2) * np.power(1 + z, 2) * (1 - endVelocity) / (1 + endVelocity)) )

	print("start: " + str(startWL))
	print("end: " + str(endWL))


	goodRange1 = np.where((spec1wl > startWL) & (spec1wl < endWL))
	goodRange2 = np.where((spec2wl > startWL) & (spec2wl < endWL))

	# we have a pretty good shot of finding strange variability so we build 
	# graph. To reduce code coupling, we make a figure object here which the
	# caller can use
	outputFigure = plt.figure(1)
	outputAxis   = outputFigure.add_subplot()
	outputAxis.grid(visible=True)

	outputAxis.set_xlim(left = (spec1wl[goodRange1])[0], right = (spec1wl[goodRange1])[-1] )

	# add red line for arbitrary cutoff
	outputAxis.plot(spec1wl[goodRange1], np.ones(len(spec1wl[goodRange1])) * MIN_TROUGH_CUTOFF, color='red', linestyle='dashed')

	outputAxis.plot(spec1wl[goodRange1], spec1flux[goodRange1], color='lightcoral', label=str(spec1.plate)+"-"+str(spec1.mjd)+"-"+str(spec1.fiber))
	outputAxis.plot(spec2wl[goodRange2], spec2flux[goodRange2], color='navy', label=str(spec2.plate)+"-"+str(spec2.mjd)+"-"+str(spec2.fiber))

	maxWL = spec1wl[maxWLindex]
	minWL = spec1wl[minWLindex]

	plt.axvline(x=maxWL)
	plt.axvline(x=minWL)

	spec1flux, spec1wl, spec1error, spec2flux, spec2wl, spec2error = intersectSpectra(spec1flux, spec1wl, spec1error, spec2flux, spec2wl, spec2error)


	worstError = np.max(np.row_stack((spec1error, spec2error)), axis=0)
	plt.plot(spec1wl, worstError, color='black', label='error')

	spec1SNR = np.mean(spec1flux/spec1error)
	spec2SNR = np.mean(spec2flux/spec2error)

	if min(spec1SNR, spec2SNR) < 20:
		outputFigure.clear()
		print("snr too low")
		return (False, (None, None, None, None))

	outputAxis.set_title(f"{indexNum}: {objectName}, z={z}, SNR={min(spec1SNR, spec2SNR)}", loc="left")

	

	troughStartIdx = np.argmin(np.abs(minWL - spec1wl))
	troughEndIdx = np.argmin(np.abs(maxWL - spec1wl))

	totalEvaluated = 0
	totalDiff      = 0
	MIN_PERCENTAGE = 0.9
	# make sure it's variable enough
	for i in range(1, 30):
		lowerIdx = troughStartIdx - i
		upperIdx = troughStartIdx + i

		if np.abs(spec1flux[lowerIdx] - spec2flux[lowerIdx]) > max(spec1error[lowerIdx], spec2error[lowerIdx]):
			totalDiff += 1
		if np.abs(spec1flux[upperIdx] - spec2flux[upperIdx]) > 2 * max(spec1error[upperIdx], spec2error[upperIdx]):
			totalDiff += 1
		totalEvaluated += 2

	if totalDiff < MIN_PERCENTAGE * totalEvaluated:
		outputFigure.clear()
		print("spectra not variable enough. rejected")
		return (False, (None, None, None, None))


	plotPotentialLines(outputAxis, minWL, maxWL, z)

	# first we need to get indices within a certain range
	# according to SDSS documentation, a variance of 0 is invalid. 1/0 in numpy
	# gives infinity, so anything with infinity error is invalid
	spec1indices = np.logical_and(spec1wl<maxWL, spec1wl>minWL, spec1error!=float('inf'))
	spec2indices = np.logical_and(spec2wl<maxWL, spec2wl>minWL, spec2error!=float('inf'))

	# TODO: make this line less shit
	spec1flux, spec1wl, spec1error, spec2flux, spec2wl, spec2error = intersectSpectra(spec1flux[spec1indices], spec1wl[spec1indices], spec1error[spec1indices], spec2flux, spec2wl, spec2error)
	
	if len(spec1flux) == 0:
		outputFigure.clear()
		print("ERROR: EMPTY FLUX ARRAY. TODO: TROUBLESHOOT")
		return (False, (None, None, None, None))

	# fluxes before the trough. TODO: Find a more sophisitcated way to calculate
	# this
	I_01 = np.mean([spec1flux[0], spec1flux[-1]])
	I_02 = np.mean([spec2flux[0], spec2flux[-1]])

	# we now make sure that the continuum is variable enough for this to be considered
	# strange variability
	

	'''
	sigma_01 = np.sqrt(np.power(spec1error[0], 2) + np.power(spec1error[-1], 2))/2
	sigma_02 = np.sqrt(np.power(spec2error[0], 2) + np.power(spec2error[-1], 2))/2
	print(sigma_01)
	print(sigma_02)
	if np.abs(I_01-I_02) < max(sigma_01, sigma_02):
		outputFigure.clear()
		return (False, (None, None, None, None))
	'''


		# this is really arbitrary, but I feel like we should look at the lowest ~40%
	# of points.
	width = np.ceil(spec1flux.size * 0.40).astype('int')
	spec1lowest = getMinPartition(spec1flux, width)
	spec2lowest = getMinPartition(spec2flux, width)
#	print("spec1lowest: " + str(spec1lowest+width))
#	print("spec2lowest: " + str(spec2lowest))
	if spec1lowest + width <= spec2lowest or spec2lowest + width <= spec1lowest:
		outputFigure.clear()
		print("no intersection of minima")
		return (False, (None, None, None, None))
	
	possibleMatchBegin = np.max((spec1lowest, spec2lowest))
	possibleMatchEnd = np.min((spec1lowest, spec2lowest)) + width

	spec1flux = spec1flux   [possibleMatchBegin:possibleMatchEnd]
	spec1wl   = spec1wl     [possibleMatchBegin:possibleMatchEnd]
	spec1error = spec1error [possibleMatchBegin:possibleMatchEnd]
	spec2flux = spec2flux   [possibleMatchBegin:possibleMatchEnd]
	spec2wl   = spec2wl     [possibleMatchBegin:possibleMatchEnd]
	spec2error = spec2error [possibleMatchBegin:possibleMatchEnd]

	# this line is kind of a doozy, but we basically look at both spectra's error
	# #'s, and we pick the smaller error for each index. It might make more sense to
	# pick the biggest S/N ratio, and multiply that by the smallest value to get a
	# more conservative estimate
	possibleMatchError = np.min(np.row_stack((spec1error, spec2error)), axis=0)
	possibleMatchDiff  = np.abs(spec1flux - spec2flux)
	
#	print(possibleMatchDiff)
#	print(possibleMatchError)

	# TODO: figure out if it should be < possibleMatchError or if it should be
	#       possibleMatchError/2
	closeEnough = possibleMatchDiff < possibleMatchError
	highEnough  = np.logical_and(spec1flux > spec1error + MIN_TROUGH_CUTOFF, spec2flux > spec2error + MIN_TROUGH_CUTOFF)
	# at least some of the points must be close enough, and all must be high enough
	if np.all(closeEnough == False) or np.any(highEnough == False):
		outputFigure.clear()
		return (False, (None, None, None, None))

	valid = np.logical_and(closeEnough, highEnough)
	print("valid: " + str(valid))

	# low error good
	# low difference good
	# so we want to minimize diff*error. There should probably be a more soph-
	# isticated way to do this...
#	print("valid: " + str(valid))
	closestMatchIdx    = np.argmin(possibleMatchDiff[valid])#*possibleMatchError[valid])
	if min((spec1flux[valid])[closestMatchIdx], (spec2flux[valid])[closestMatchIdx]) < MIN_TROUGH_CUTOFF:
		outputFigure.clear()
		return (False, (None, None, None, None))

	displayPatch(outputAxis, spec1wl, spec1flux, 0, len(spec1wl))
	#if int(spec1.mjd)==56574 and int(spec1.plate)==6598 and int(spec1.fiber)==271:
	#plt.show()
	#else:
	outputFigure.legend()
	pp.savefig()
	plt.clf()
	#outputFigure.show()
	return (True, (I_01, I_02, (spec1flux[valid])[closestMatchIdx], (spec2flux[valid])[closestMatchIdx]))



if len(sys.argv) != 4:
	print('Incorrect # of arguments:\n \npython ProgramName absorptionFile.csv duplicates.csv SDSSfolder')



absorptionCsvFilename = sys.argv[1]
duplicateCsvFilename = sys.argv[2]
SDSSfolder = sys.argv[3]
SDSSfolder2= sys.argv[4]
df = pd.read_csv(absorptionCsvFilename)
duplicatesDF = pd.read_csv(duplicateCsvFilename)

# kind of janky but we filter out everything that doesn't have a defined vmin
# .str.len() > 2 basically converts into string and checks length. If it's an
# empty array that's [] (2 characters) so by enasuring there's > 2, we know
# it's not empty
df = df[df['VMINS'].str.len()>2]

spectraWithAbsorption = [] # all spectra with significant absorption
spectraWithVMIN = defaultdict(list) # vmin(s) for each absorbing spectra
spectraWithVMAX = defaultdict(list) # vmax(s) for each absorbing spectra
spectraRedshifts= defaultdict(list) # z-values
# we basically just extract each spectra's MJD, fiber, and plate. Then, we get
# the absorption code's calculated start/end BAL region for this
for index, row in df.iterrows():
	filename = row['NORM SPECTRA FILE NAME']
	nums = re.findall('\d+', filename)
	key = Spectra(plate=nums[0], mjd=nums[1], fiber=nums[2])
	vmins = row['VMINS']
	vmaxs = row['VMAXS']
	vmins = extractIntegerArrayFromCSV(vmins)
	vmaxs = extractIntegerArrayFromCSV(vmaxs)

	spectraWithAbsorption.append(key)
	spectraWithVMIN[key] = vmins
	spectraWithVMAX[key]  = vmaxs
	spectraRedshifts[key] = row['z']

duplicateSpectra = []
objectNames =      []
specNum     = 0
# we then look for repeated observations of BALs
for index, duplicateRow in duplicatesDF.iterrows():
	duplicates = duplicateRow['duplicate_spectra'].split(",")
	print(duplicates)
	duplicates = list(map(extractSpectraTupleFromString, duplicates))
	print(duplicates)

	tmp = [] # if this set of duplicates has at least 2 spectra with absorption
	         # we save it to duplicateSpectra
	for spectra in duplicates:
		#if os.access(os.path.join(SDSSfolder, toDered(spectra)), os.R_OK) or os.access(os.path.join(SDSSfolder2, toDered(spectra)), os.R_OK):
		tmp.append(spectra)
	if len(tmp)>= 2:
		specNum = specNum + 1
		duplicateSpectra.append(tmp)
		objectNames.append(duplicateRow['sdssName'])
		print(str(duplicateRow['sdssName']) + ": " + str(tmp))

print(sum(len(l) for l in duplicateSpectra))
outputDF = pd.DataFrame(columns=['objectName', 'spec1', 'spec2', 'I01', 'I02', 'I1', 'I2'])
idx=0
foundStrange = False

print("specNum: " + str(specNum))
# now we look for strange variability in all the duplicate spectra
for specListIndex, spectraList in enumerate(duplicateSpectra):
	print ("new set")
	# daily reminder that it goes vmax -> vmin due to how blueshifting works
	maxWLIndices = spectraWithVMIN[spectraList[0]]
	minWLIndices = spectraWithVMAX[spectraList[0]]
	redshift     = spectraRedshifts[spectraList[0]]

	generatedFigures = [] # temporary storage before saving plots to the pdf

	# we assume that spectraList[0] has the BAL. We test that spectra against
	# every duplicate observation spectraList[1:]
	for i in range(1, len(spectraList)):
		print("looking for strange variability...")
		for minWLIndex, maxWLIndex in zip(minWLIndices, maxWLIndices):
			tmp = findStats(spectraList[0], spectraList[i], minWLIndex, maxWLIndex, objectNames[specListIndex], idx, redshift)
			# i have no fucking clue why concat is the "standard" now, but this
			# is what the pandas dev team wants apparently
			if tmp[0] == True:
				foundStrange = True
				idx=idx+1
				outputDF = pd.concat([outputDF, pd.DataFrame([{'objectName': objectNames[specListIndex],'spec1': spectraList[0], 'spec2': spectraList[1], 'I01': tmp[1][0], 'I02': tmp[1][1], 'I1': tmp[1][2], 'I2': tmp[1][3]}])])
	

outputDF = outputDF.reset_index()
outputDF.to_csv('out.csv')
print(str(outputDF))
pp.close()
print("--- %s seconds ---" % (time.time() - start_time))
