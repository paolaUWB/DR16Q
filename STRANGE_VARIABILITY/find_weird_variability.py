import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import re
import os

from utility_functions  import read_spectra
from collections import defaultdict, namedtuple
from strange_utilities import getMinPartition, getSpectrum, intersectSpectra
from graph_utilities import graphMain
import time
import scipy.constants as sc

MIN_TROUGH_CUTOFF = 0.0
SNR_CUTOFF = 6
GRAPH_RESULTS = True

# pulled from http://astronomy.nmsu.edu/drewski/tableofemissionlines.html
# all in angstroms
CIIRest    = 1335.120
CIVRest    = 1549.0524
SiIVRest   = 1396.747
OIRest     = 1302.168

start_time = time.time()

Spectra = namedtuple("Spectra", ["plate", "mjd", "fiber"])

pp = PdfPages('test.pdf')

'''
# Takes a string like 6742-124359-0943 and converts it to a string tuple
'''
def extractSpectraTupleFromString(inputString):
	nums = re.findall('\d+', inputString)
	return Spectra(plate=nums[0], mjd=nums[1], fiber=nums[2])

'''
# Looks to the left and to the right of lineStartIdx and lineEndIdx respectively
# and checks to make sure that flux1 and flux2's difference is greater than those
# pixel's error values. Points is the # of points checked, percentage is the min
# percentage that has to be different
'''
def isSpectrumVariable(flux1: np.ndarray, 
                       flux2: np.ndarray, 
					   error: np.ndarray, 
					   lineStartIdx: float, 
					   lineEndIdx: float, 
					   points: int, 
					   percentage: float) -> bool:
	totalEvaluated = 0
	totalDiff      = 0
	MIN_PERCENTAGE = 0.6
	# make sure it's variable enough
	for i in range(1, int(np.ceil(points/2))):
		lowerIdx = lineStartIdx - i
		upperIdx = lineEndIdx + i

		if np.abs(flux1[lowerIdx] - flux2[lowerIdx]) >= error[lowerIdx]:
			totalDiff += 1
		if np.abs(flux1[upperIdx] - flux2[upperIdx]) >= error[lowerIdx]:
			totalDiff += 1
		totalEvaluated += 2
	return totalDiff > (totalEvaluated * percentage)

# assume indices correspond to spec1
# returns:
#       (isWeird, (I_01, I_01, I_1, I_2), bottomWavelength)
#types:  boolean  float float float float float
# if isWeird is False, all other values default to 0.0
def findStats(spec1flux: np.ndarray, 
              spec1wl: np.ndarray, 
			  spec1error: np.ndarray, 
			  spec2flux: np.ndarray, 
			  spec2wl: np.ndarray, 
			  spec2error: np.ndarray, 
			  minWL: float, 
			  maxWL: float) -> tuple[bool, tuple[float, float, float, float], float]:
	troughStartIdx = np.argmin(np.abs(minWL - spec1wl))
	troughEndIdx = np.argmin(np.abs(maxWL - spec1wl))

	#	worstCaseError = np.max(np.row_stack((spec1error, spec2error)), axis=0)
	#	isVariable = isSpectrumVariable(spec1flux, spec2flux, worstCaseError, troughStartIdx, troughEndIdx, 30, 0.6) 

	# first we need to get indices within a certain range
	# according to SDSS documentation, a variance of 0 is invalid. 1/0 in numpy
	# gives infinity, so anything with infinity error is invalid
	indices = np.logical_and(spec1wl<maxWL, spec1wl>minWL)

	spec1flux, spec1wl, spec1error = spec1flux[indices], spec1wl[indices], spec1error[indices]
	spec2flux, spec2wl, spec2error = spec2flux[indices], spec2wl[indices], spec2error[indices]

	if len(spec1flux) == 0:
		print("ERROR: EMPTY FLUX ARRAY. TODO: TROUBLESHOOT")
		return (False, (0, 0, 0, 0), 0)

	# fluxes before the trough. TODO: Find a more sophisitcated way to calculate
	# this
	I_01 = np.mean([spec1flux[0], spec1flux[-1]])
	I_02 = np.mean([spec2flux[0], spec2flux[-1]])

	# this is really arbitrary, but I feel like we should look at the lowest ~40%
	# of points.
	width = np.ceil(spec1flux.size * 0.40).astype('int')
	spec1lowest = getMinPartition(spec1flux, width)
	spec2lowest = getMinPartition(spec2flux, width)
#	print("spec1lowest: " + str(spec1lowest+width))
#	print("spec2lowest: " + str(spec2lowest))
	if spec1lowest + width <= spec2lowest or spec2lowest + width <= spec1lowest:
		print("no intersection of minima")
		return (False, (0, 0, 0, 0), 0)
	
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
	
	# TODO: figure out if it should be < possibleMatchError or if it should be
	#       possibleMatchError/2
	closeEnough = possibleMatchDiff < possibleMatchError
	highEnough  = np.logical_and(spec1flux > spec1error + MIN_TROUGH_CUTOFF, spec2flux > spec2error + MIN_TROUGH_CUTOFF)
	# at least some of the points must be close enough, and all must be high enough
	if np.all(closeEnough == False) or np.any(highEnough == False):
		return (False, (0, 0, 0, 0), 0)

	valid = np.logical_and(closeEnough, highEnough)

	# low error good
	# low difference good
	# so we want to minimize diff*error. There should probably be a more soph-
	# isticated way to do this...
	closestMatchIdx    = np.argmin(possibleMatchDiff[valid])
	if min((spec1flux[valid])[closestMatchIdx], (spec2flux[valid])[closestMatchIdx]) < MIN_TROUGH_CUTOFF:
		return (False, (0, 0, 0, 0), 0)


	return (True, (I_01, I_02, (spec1flux[valid])[closestMatchIdx], (spec2flux[valid])[closestMatchIdx]), (spec1wl[valid])[closestMatchIdx])



if len(sys.argv) != 3:
	print('Incorrect # of arguments:\n \npython ProgramName absorptionFile.csv duplicates.csv')



absorptionCsvFilename = sys.argv[1]
duplicateCsvFilename = sys.argv[2]
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
	print(str(row))
	filename = row['NORM SPECTRA FILE NAME']
	nums = re.findall('\d+', filename)
	key = Spectra(plate=nums[0], mjd=nums[1], fiber=nums[2])
	vmins = row['VMINS']
	vmaxs = row['VMAXS']
	vmins = np.fromstring(vmins.strip("["), dtype=float, sep=', ')
	vmaxs = np.fromstring(vmaxs.strip("["), dtype=float, sep=', ')

	spectraWithAbsorption.append(key)
	spectraWithVMIN[key] = vmins
	spectraWithVMAX[key]  = vmaxs
	spectraRedshifts[key] = row['z']

duplicateSpectra = []
objectNames =      []
specNum     = 0
# we then look for repeated observations of BALs
for index, duplicateRow in duplicatesDF.iterrows():
	print("row")
	duplicates = duplicateRow['duplicate_spectra'].split(",")
	print(duplicates)
	duplicates = list(map(extractSpectraTupleFromString, duplicates))
	print(duplicates)

	tmp = [] # if this set of duplicates has at least 2 spectra with absorption
	         # we save it to duplicateSpectra
	for spectra in duplicates:
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
		for minWL, maxWL in zip(minWLIndices, maxWLIndices):
			f1, wl1, e1 = getSpectrum(spectraList[0].plate, spectraList[0].mjd, spectraList[0].fiber)
			f2, wl2, e2 = getSpectrum(spectraList[i].plate, spectraList[i].mjd, spectraList[i].fiber)
			f1, wl1, e1, f2, wl2, e2 =  intersectSpectra(f1, wl1, e1, f2, wl2, e2)

			SNR1 = round(np.mean(f1/e1), 2)
			SNR2 = round(np.mean(f2/e2), 2)
			if min(SNR1, SNR2) < SNR_CUTOFF:
				print("snr too low")
				continue

			minError = np.min(np.row_stack((e1, e2)), axis=0)
			lineStartIdx = int(np.argmin(np.abs(minWL - wl1)))
			lineEndIdx   = int(np.argmin(np.abs(maxWL - wl1)))
			if not isSpectrumVariable(f1, f2, minError, lineStartIdx, lineEndIdx, 30, 0.6):
				print("spectrum not variable enough")
				continue
			

			tmp = findStats(f1, wl1, e1, f2, wl2, e2, minWL, maxWL)

			# i have no fucking clue why concat is the "standard" now, but this
			# is what the pandas dev team wants apparently
			if tmp[0] == True:
				if GRAPH_RESULTS:
					plt.clf()
					generatedFig = graphMain(f1, wl1, e1, f2, wl2, e2, minWL, maxWL, tmp[2], redshift)
					generatedFig.gca().set_title(f"{idx}: {objectNames[specListIndex]} z={redshift} snr={str(round(min(SNR1, SNR2), 2))}")
					pp.savefig()
				foundStrange = True
				idx=idx+1
				outputDF = pd.concat([outputDF, pd.DataFrame([{'objectName': objectNames[specListIndex],'spec1': spectraList[0], 'spec2': spectraList[1], 'I01': tmp[1][0], 'I02': tmp[1][1], 'I1': tmp[1][2], 'I2': tmp[1][3]}])])
	

outputDF = outputDF.reset_index()
outputDF.to_csv('out.csv')
print(str(outputDF))
pp.close()
print("--- %s seconds ---" % (time.time() - start_time))
