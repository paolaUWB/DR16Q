# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:20:48 2024

@author: annam
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import stats
import pandas as pd
from sample_bisecting_function import bisector
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
################## Changeable Variables ######################################################################################################

# want to plot bisectors of EHVO sample
plot_bisectors = True
#plot_bisectors = False

hist_parent = True
#hist_parent = 'Rankineall'
#hist_parent = 'parent bal sep'

colormap = 'CIV Blueshift'
colormap = 'CIV Distance'
colormap = 'CIV EW'
############################################################################################################################################################################################################


#XQR-30 Properties File
xqrdata = pd.read_csv('https://raw.githubusercontent.com/paolaUWB/DR16Q/Anna/CROSS_CORRELATION/XQR30_Property_Data.csv')

#Rankines info file PRH selected the good cases
infoRankineparent = 'https://raw.githubusercontent.com/paolaUWB/DR16Q/master/CROSS_CORRELATION/DR16parent_DR14RankineInfo.csv'
infoRankineEHVO = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv"
infoRankineparent_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9parent_DR14RankineInfo.csv"
infoRankineEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv"
RH_parent_16= 'https://raw.githubusercontent.com/paolaUWB/DR16Q/master/CROSS_CORRELATION/DR16_parent_sample.csv'

#Extracting values from Rankines info file for EHVO Parent
dfRPA = pd.read_csv(infoRankineparent, header=None)
dfRHV = pd.read_csv(infoRankineEHVO, header=0)
dfRPA_9 = pd.read_csv(infoRankineparent_DR9, header=None)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=0)

#Extracting Values from XQR-30 Properties Files
xqr_loglbol = np.array(xqrdata["Log_Lbol"])
xqr_logmbh = np.array(xqrdata["Log_bhm_civ"])
xqr_edd = np.array(xqrdata["edd_civ"])
xqr_logedd = np.log10(xqr_edd)
xqr_EHVOstatus = np.array(xqrdata["status"])
xqr_bs = np.array(xqrdata["CIV_Emission_Blueshift"])

#XQR30: Separating Quasars by type using the midZstatus column of the csv
xqr_EHVO = (xqr_EHVOstatus == "EHVO")
xqr_nonBAL = (xqr_EHVOstatus == "nBAL")
xqr_BAL = (xqr_EHVOstatus == 'BAL')
xqr_unknown = (xqr_EHVOstatus == "N/A")

def safe_log10(arr):
  arr = np.array(arr)
  arr[arr <= 0] = np.nan
  return np.log10(arr)


#EHVO Physical Properties-----------------------------------------------------------
xqr_lbol_EHVO = xqr_loglbol[xqr_EHVO]
xqr_mbh_EHVO = xqr_logmbh[xqr_EHVO]
xqr_edd_EHVO = xqr_logedd[xqr_EHVO]
xqr_bs_EHVO = xqr_bs[xqr_EHVO]



#Parent Physical Properties-----------------------------------------------------------

xqr_nonBALlbol = xqr_loglbol[xqr_nonBAL]
xqr_nonBALmbh = xqr_logmbh[xqr_nonBAL]
xqr_nonBALedd = xqr_logedd[xqr_nonBAL]
xqr_bs_nonBAL = xqr_bs[xqr_nonBAL]
xqr_BALlbol = xqr_loglbol[xqr_BAL]
xqr_BALmbh = xqr_logmbh[xqr_BAL]
xqr_BALedd = xqr_logedd[xqr_BAL]
xqr_bs_BAL = xqr_bs[xqr_BAL]


xqr_lbol_Parent = np.concatenate([xqr_nonBALlbol, xqr_BALlbol])
xqr_mbh_Parent = np.concatenate([xqr_nonBALmbh, xqr_BALmbh])
xqr_edd_Parent = np.concatenate([xqr_nonBALedd, xqr_BALedd])
xqr_bs_Parent = np.concatenate([xqr_bs_nonBAL, xqr_bs_BAL])

xqr_bs_both = np.concatenate([xqr_bs_Parent, xqr_bs_EHVO])

###DR16###
#for parent sample
Parentin16Rank_mbh=dfRPA[dfRPA.columns[16]].to_numpy()
Parentin16Rank_lbol=dfRPA[dfRPA.columns[17]].to_numpy()
Parentin16Rank_redd=dfRPA[dfRPA.columns[18]].to_numpy()
Parentin16Rank_CivBlue=dfRPA[dfRPA.columns[5]].to_numpy()
Parentin16Rank_CivEW=dfRPA[dfRPA.columns[6]].to_numpy()
Parentin16Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()
Parentin16Rank_Civdist=dfRPA[dfRPA.columns[14]].to_numpy()



  #for EHVOs
EHVOin16Rank_mbh=dfRHV[dfRHV.columns[16]].to_numpy()
EHVOin16Rank_lbol=dfRHV[dfRHV.columns[17]].to_numpy()
EHVOin16Rank_redd=dfRHV[dfRHV.columns[18]].to_numpy()
EHVOin16Rank_CivBlue=dfRHV[dfRHV.columns[5]].to_numpy()
EHVOin16Rank_CivEW=dfRHV[dfRHV.columns[6]].to_numpy()
EHVOin16Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy() #LEF: Why is this dfRPA
EHVOin16Rank_Civdist=dfRHV[dfRHV.columns[14]].to_numpy()


###DR9###

#for parent sample
Parentin9Rank_mbh=dfRPA_9[dfRPA_9.columns[16]].to_numpy()
Parentin9Rank_lbol=dfRPA_9[dfRPA_9.columns[17]].to_numpy()
Parentin9Rank_redd=dfRPA_9[dfRPA_9.columns[18]].to_numpy()
Parentin9Rank_CivBlue=dfRPA_9[dfRPA_9.columns[5]].to_numpy()
Parentin9Rank_CivEW=dfRPA_9[dfRPA_9.columns[6]].to_numpy()
Parentin9Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()
Parentin9Rank_Civdist=dfRPA[dfRPA_9.columns[14]].to_numpy()

#for EHVOs
EHVOin9Rank_mbh=dfRHV_9[dfRHV_9.columns[16]].to_numpy()
EHVOin9Rank_lbol=dfRHV_9[dfRHV_9.columns[17]].to_numpy()
EHVOin9Rank_redd=dfRHV_9[dfRHV_9.columns[18]].to_numpy()
EHVOin9Rank_CivBlue=dfRHV_9[dfRHV_9.columns[5]].to_numpy()
EHVOin9Rank_CivEW=dfRHV_9[dfRHV_9.columns[6]].to_numpy()
EHVOin9Rank_HeIIEW=dfRPA[dfRPA.columns[7]].to_numpy()
EHVOin9Rank_Civdist=dfRHV_9[dfRHV_9.columns[14]].to_numpy()


###Rankine all info for non-BAL and BAL sample

rankineAll = os.getcwd() + "/../Rankine20_CIV_BAL.fits"

#Reading RankineAll fit table for Emission Parameters
rankine_all = fits.open(rankineAll) # reading in Rankine20_CIV_BAL.fits data
data_rankAll = rankine_all[1].data
CivBlue_rankAll = data_rankAll['CIV_blue']
CivEW_rankAll = data_rankAll['CIV_EW']
HeIIBlue_rankAll = data_rankAll['HeII_blue']
HeIIEW_rankAll = data_rankAll['HeII_EW']
CivDistance_rankAll = data_rankAll['CIV_distance']
SNR_rankAll = data_rankAll['S/N']
bi_bi = data_rankAll['BI_BI']
vmin_rankAll = data_rankAll['BI_VMIN']
vmax_rankAll = data_rankAll['BI_VMAX']
good = data_rankAll['good']
LogLbol_rankAll = data_rankAll['LOG_LBOL']
sdss_name_rankAll = data_rankAll['SDSS_NAME']
z_rankAll = data_rankAll['z']
edd_rankAll = data_rankAll['LOG_REDD']
MBH_rankAll = data_rankAll['LOG_BHM_CIV']


good_rankAll=np.where(good==1)
snr_10 = np.where((SNR_rankAll > 10)&(good==1))
bal_rankAll = np.where((bi_bi > 0)&(SNR_rankAll > 10)&(good == 1))
nonbal_rankAll = np.where((bi_bi==0)&(SNR_rankAll > 10)&(good==1))

CivBlue_rankAll_SNR = CivBlue_rankAll[snr_10]
CivBlue_rankAll_good = CivBlue_rankAll[good_rankAll]
CivEW_rankAll_good = CivEW_rankAll[good_rankAll]
CivEW_rankAll_SNR = CivEW_rankAll[snr_10]
CivBlue_rankAll_bal = CivBlue_rankAll[bal_rankAll]
CivBlue_rankAll_nonbal = CivBlue_rankAll[nonbal_rankAll]
CivEW_rankAll_bal = CivEW_rankAll[bal_rankAll]
CivEW_rankAll_nonbal = CivEW_rankAll[nonbal_rankAll]
HeIIBlue_rankAll_bal = HeIIBlue_rankAll[bal_rankAll]
HeIIBlue_rankAll_nonbal = HeIIBlue_rankAll[nonbal_rankAll]
HeIIEW_rankAll_bal = HeIIEW_rankAll[bal_rankAll]
HeIIEW_rankAll_nonbal = HeIIEW_rankAll[nonbal_rankAll]
CivDistance_rankAll_bal = CivDistance_rankAll[bal_rankAll]
CivDistance_rankAll_nonbal = CivDistance_rankAll[nonbal_rankAll]
LogLbol_rankAll_bal = LogLbol_rankAll[bal_rankAll]
LogLbol_rankAll_nonbal = LogLbol_rankAll[nonbal_rankAll]
vmin_rankAll_bal = vmin_rankAll[bal_rankAll]
sdss_rankAll_bal = sdss_name_rankAll[bal_rankAll]
z_rankAll_bal = z_rankAll[bal_rankAll]
vmax_rankAll_bal = vmax_rankAll[bal_rankAll]
edd_rankAll_bal = edd_rankAll[bal_rankAll]
edd_rankAll_nonbal = edd_rankAll[nonbal_rankAll]
MBH_rankAll_bal = MBH_rankAll[bal_rankAll]
MBH_rankAll_nonbal = MBH_rankAll[nonbal_rankAll]



MBH_rankAll_nonbal = MBH_rankAll_nonbal[np.isfinite(MBH_rankAll_nonbal) & (CivBlue_rankAll_nonbal > 500)]
MBH_rankAll_bal = MBH_rankAll_bal[np.isfinite(MBH_rankAll_bal) & (CivBlue_rankAll_bal > 500)]

edd_rankAll_nonbal = edd_rankAll_nonbal[np.isfinite(edd_rankAll_nonbal) & (CivBlue_rankAll_nonbal > 500)]
edd_rankAll_bal = edd_rankAll_bal[np.isfinite(edd_rankAll_bal) & (CivBlue_rankAll_bal > 500)]


CivBlueCut_nonbal = np.where(CivBlue_rankAll_nonbal>500)
CivBlueCut_bal = np.where(CivBlue_rankAll_bal>500)


###Combined DR9 and DR16 Arrays###

bi_bi0 = dfRPA[dfRPA.columns[10]].to_numpy()
SNR = dfRPA[dfRPA.columns[15]].to_numpy()
bi_bi0_9 = dfRPA_9[dfRPA_9.columns[10]].to_numpy()
SNR_9 = dfRPA_9[dfRPA_9.columns[15]].to_numpy()
#non_bal = np.where((bi_bi0==0)&(SNR > 10))
#pos_bal = np.where((bi_bi0>0) & (SNR > 10))
#non_bal_9 = np.where((bi_bi0_9==0)&(SNR_9 > 10))
#pos_bal_9 = np.where((bi_bi0_9>0) & (SNR_9 > 10))

non_bal = np.where((bi_bi0==0))
pos_bal = np.where((bi_bi0>0))
non_bal_9 = np.where((bi_bi0_9==0))
pos_bal_9 = np.where((bi_bi0_9>0))

#Parent Sample
parentboth_mbh = np.concatenate([Parentin9Rank_mbh, Parentin16Rank_mbh])
parentboth_lbol = np.concatenate([Parentin9Rank_lbol, Parentin16Rank_lbol])
parentboth_redd = np.concatenate([Parentin9Rank_redd, Parentin16Rank_redd])
parentboth_CivBlue = np.concatenate([Parentin9Rank_CivBlue, Parentin16Rank_CivBlue])
parentboth_CivEW = np.concatenate([Parentin9Rank_CivEW, Parentin16Rank_CivEW])
parentboth_HeIIEW = np.concatenate([Parentin9Rank_HeIIEW, Parentin16Rank_HeIIEW])
parentboth_logHeIIEW = np.log10(parentboth_HeIIEW)
parentboth_CivDist = np.concatenate([Parentin9Rank_Civdist, Parentin16Rank_Civdist])

parentboth_nonbal_mbh = np.concatenate([Parentin9Rank_mbh[non_bal_9], Parentin16Rank_mbh[non_bal]])
parentboth_nonbal_lbol = np.concatenate([Parentin9Rank_lbol[non_bal_9], Parentin16Rank_lbol[non_bal]])
parentboth_nonbal_redd = np.concatenate([Parentin9Rank_redd[non_bal_9], Parentin16Rank_redd[non_bal]])

parentboth_bal_mbh = np.concatenate([Parentin9Rank_mbh[pos_bal_9], Parentin16Rank_mbh[pos_bal]])
parentboth_bal_lbol = np.concatenate([Parentin9Rank_lbol[pos_bal_9], Parentin16Rank_lbol[pos_bal]])
parentboth_bal_redd = np.concatenate([Parentin9Rank_redd[pos_bal_9], Parentin16Rank_redd[pos_bal]])

#print('parent non BAL: ', len(parentboth_bal_lbol))
#print('parent BAL: ', len(parentboth_nonbal_lbol))
#print('parent total: ', len(parentboth_lbol))


#EHVOs
EHVOboth_mbh = np.concatenate([EHVOin9Rank_mbh, EHVOin16Rank_mbh])
EHVOboth_lbol = np.concatenate([EHVOin9Rank_lbol, EHVOin16Rank_lbol])
EHVOboth_redd = np.concatenate([EHVOin9Rank_redd, EHVOin16Rank_redd])
EHVOboth_CivBlue = np.concatenate([EHVOin9Rank_CivBlue, EHVOin16Rank_CivBlue])
EHVOboth_CivEW = np.concatenate([EHVOin9Rank_CivEW, EHVOin16Rank_CivEW])
EHVOboth_HeIIEW = np.concatenate([EHVOin9Rank_HeIIEW, EHVOin16Rank_HeIIEW])
EHVOboth_logHeIIEW = np.log10(EHVOboth_HeIIEW)
EHVOboth_CivDist = np.concatenate([EHVOin9Rank_Civdist, EHVOin16Rank_Civdist])

#'Good' Values, Coatman et al. 2017 cutoff for Black hole mass (and thus eddington ratio) plots
bb9 = np.where(Parentin9Rank_CivBlue > 500)
bb9_EHVO = np.where(EHVOin9Rank_CivBlue > 500)
bb16 = np.where(Parentin16Rank_CivBlue > 500)
bb16_EHVO = np.where(EHVOin16Rank_CivBlue > 500)
bb_both = np.where(parentboth_CivBlue > 500)
bbboth_EHVO = np.where(EHVOboth_CivBlue > 500) # using this also removes the outlier!!

bbXQR_nonBAL = np.where(xqr_bs_nonBAL > 500)
bbXQR_BAL = np.where(xqr_bs_BAL > 500)
bbXQR_EHVO = np.where(xqr_bs_EHVO > 500)

highbsEHVO = np.where(EHVOboth_CivBlue > 2500) #not used


# LEF: ARE THESE SECTIONS ACTUALLY USED AND NEEDED?
#Combined EHVO and Parent Blueshift Arrays for Color Map
c_DR16_CivBlue = np.concatenate([Parentin16Rank_CivBlue[bb16], EHVOin16Rank_CivBlue[bb16_EHVO]])
c_DR9_CivBlue = np.concatenate([Parentin9Rank_CivBlue[bb9], EHVOin9Rank_CivBlue])
c_both_CivBlue = np.concatenate([parentboth_CivBlue[bb_both], EHVOboth_CivBlue[bbboth_EHVO]])

#COMBINED EHVO and Parent CIV Equivalent Width Arrays for Color Map
c_DR16_CivEW = np.concatenate([Parentin16Rank_CivEW[bb16], EHVOin16Rank_CivEW[bb16_EHVO]])
c_DR9_CivEW = np.concatenate([Parentin9Rank_CivEW[bb9], EHVOin9Rank_CivEW])
c_both_CivEW = np.concatenate([parentboth_CivEW[bb_both], EHVOboth_CivEW[bbboth_EHVO]])

#COMBINED EHVO and Parent He II Equivalent Width Arrays for Color Map
c_DR16_HeIIEW = np.concatenate([Parentin16Rank_HeIIEW[bb16], EHVOin16Rank_HeIIEW[bb16_EHVO]])
c_DR9_HeIIEW = np.concatenate([Parentin9Rank_HeIIEW[bb9], EHVOin9Rank_HeIIEW])
c_both_HeIIEW = np.concatenate([parentboth_HeIIEW[bb_both], EHVOboth_HeIIEW[bbboth_EHVO]])
c_both_logHeIIEW = safe_log10(c_both_HeIIEW)

#COMBINED EHVO and Parent CIV Distance Arrays for Color Map
c_DR16_CivDist = np.concatenate([Parentin16Rank_Civdist[bb16], EHVOin16Rank_Civdist[bb16_EHVO]])
c_DR9_CivDist = np.concatenate([Parentin9Rank_Civdist[bb9], EHVOin9Rank_Civdist[bb9_EHVO]])
c_both_CivDist = np.concatenate([parentboth_CivDist[bb_both], EHVOboth_CivDist[bbboth_EHVO]])


#Color Map
cm = plt.get_cmap('viridis')
cmap = plt.get_cmap('viridis')


if colormap == 'CIV Blueshift':
    cbar_label = r'$\mathrm{C\,IV}$ blueshift [km s$^{-1}$]'
    c_values = EHVOboth_CivBlue[bbboth_EHVO]
    
    C_values = parentboth_CivBlue[bb_both]
    vmax=max(np.concatenate([EHVOboth_CivBlue[bbboth_EHVO], parentboth_CivBlue[bb_both]]))
    vmin=500 
    
elif colormap == 'CIV Distance':
    cbar_label = r'$\mathrm{C\,IV}$ Distance'
    c_values = EHVOboth_CivDist[bbboth_EHVO]
    
    C_values = parentboth_CivDist[bb_both]
    vmax=max(np.concatenate([EHVOboth_CivDist[bbboth_EHVO], parentboth_CivDist[bb_both]]))
    vmin=min(np.concatenate([EHVOboth_CivDist[bbboth_EHVO], parentboth_CivDist[bb_both]]))
    
elif colormap == 'CIV EW':
    cbar_label = r'$\log_{10}(\mathrm{C\,IV}\ \mathrm{EW\ [\AA]})$'
    c_values = np.log10(EHVOboth_CivEW[bbboth_EHVO])
    
    C_values = np.log10(parentboth_CivEW[bb_both])
    vmax=max(np.concatenate([np.log10(EHVOboth_CivEW[bbboth_EHVO]), np.log10(parentboth_CivEW[bb_both])]))
    vmin=min(np.concatenate([np.log10(EHVOboth_CivEW[bbboth_EHVO]), np.log10(parentboth_CivEW[bb_both])]))
    
else:
    print('**Select colormap option at the top settings!**')
    cbar_label = '???'

from scipy.stats import mannwhitneyu

'''
def compare_ehvo_hexbins(cm, x_ehvo, y_ehvo, verbose=True):
    """
    Compare hexbin C IV blueshift values in regions with and without EHVOs.

    Parameters
    ----------
    cm : matplotlib.collections.PolyCollection
        The hexbin object returned by ax.hexbin().
    x_ehvo, y_ehvo : array-like
        Coordinates of EHVO objects in the same coordinate system as the hexbin.

    Returns
    -------
    results : dict
        Statistics comparing hexagons containing EHVOs vs. hexagons without EHVOs.
    """

    hexpositions = cm.get_offsets()
    hexvalues = cm.get_array()

    # Remove hexagons whose value is NaN
    valid = np.isfinite(hexvalues)

    hexpositions = hexpositions[valid]
    hexvalues = np.asarray(hexvalues)[valid]

    # ---------------------------------------------------------
    # Find the closest hexbin to every EHVO point
    # ---------------------------------------------------------

    ehvo_hex_indices = []

    for x, y in zip(x_ehvo, y_ehvo):

        distances = (
            (hexpositions[:, 0] - x)**2 +
            (hexpositions[:, 1] - y)**2
        )

        idx = np.argmin(distances)

        ehvo_hex_indices.append(idx)

    # Only keep unique hexagons
    ehvo_hex_indices = np.unique(ehvo_hex_indices)

    # ---------------------------------------------------------
    # Split hexagons into EHVO and non-EHVO regions
    # ---------------------------------------------------------

    ehvo_hex_mask = np.zeros(len(hexvalues), dtype=bool)
    ehvo_hex_mask[ehvo_hex_indices] = True

    values_ehvo = hexvalues[ehvo_hex_mask]
    values_no_ehvo = hexvalues[~ehvo_hex_mask]

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------

    median_ehvo = np.median(values_ehvo)
    median_no_ehvo = np.median(values_no_ehvo)

    mean_ehvo = np.mean(values_ehvo)
    mean_no_ehvo = np.mean(values_no_ehvo)

    difference = median_ehvo - median_no_ehvo

    percent_difference = (
        difference / median_no_ehvo * 100
    )

    

    results = {
        'n_ehvo_hex': len(values_ehvo),
        'n_no_ehvo_hex': len(values_no_ehvo),

        'median_ehvo': median_ehvo,
        'median_no_ehvo': median_no_ehvo,

        'mean_ehvo': mean_ehvo,
        'mean_no_ehvo': mean_no_ehvo,

        'median_difference': difference,
        'percent_difference': percent_difference,

    }

    if verbose:
        print("\n------------------------------")
        print("EHVO HEXBIN COMPARISON")
        print("------------------------------")
        print(f"EHVO hexagons:     {len(values_ehvo)}")
        print(f"Non-EHVO hexagons: {len(values_no_ehvo)}")
        print()
        print(f"Median EHVO:       {median_ehvo:.1f}")
        print(f"Median non-EHVO:   {median_no_ehvo:.1f}")
        print(f"Difference:        {difference:.1f}")
        print(f"Percent difference: {percent_difference:.1f}%")
        print("------------------------------")
        print()

    return results
'''

def compare_ehvo_hexbins(hb, x_ehvo, y_ehvo, verbose=True):
    """
    Compare hexbin values in regions with and without EHVOs.

    Parameters
    ----------
    hb : matplotlib.collections.PolyCollection
        The collection returned by ax.hexbin().
    x_ehvo, y_ehvo : array-like
        Coordinates of EHVO objects in the same coordinate system
        as the hexbin.

    Returns
    -------
    results : dict
        Statistics comparing EHVO and non-EHVO hexagons.
    """
    import numpy as np

    # Extract hexagon centers and associated values
    hexpositions = hb.get_offsets()
    hexvalues = np.asarray(hb.get_array())

    # Remove invalid hexagon values and positions
    valid = (
        np.isfinite(hexvalues)
        & np.isfinite(hexpositions).all(axis=1)
    )
    hexpositions = hexpositions[valid]
    hexvalues = hexvalues[valid]

    x_ehvo = np.asarray(x_ehvo)
    y_ehvo = np.asarray(y_ehvo)

    # Remove invalid EHVO coordinates
    valid_ehvo = np.isfinite(x_ehvo) & np.isfinite(y_ehvo)
    x_ehvo = x_ehvo[valid_ehvo]
    y_ehvo = y_ehvo[valid_ehvo]

    if len(hexvalues) == 0:
        raise ValueError("No valid hexbin values available.")

    if len(x_ehvo) == 0:
        raise ValueError("No valid EHVO coordinates provided.")

    # Find the nearest hexagon center for every EHVO
    ehvo_hex_indices = []

    for x, y in zip(x_ehvo, y_ehvo):
        distances = (
            (hexpositions[:, 0] - x)**2
            + (hexpositions[:, 1] - y)**2
        )
        ehvo_hex_indices.append(np.argmin(distances))

    # Keep unique EHVO-associated hexagons
    ehvo_hex_indices = np.unique(ehvo_hex_indices)

    # Split hexagons into EHVO and non-EHVO groups
    ehvo_hex_mask = np.zeros(len(hexvalues), dtype=bool)
    ehvo_hex_mask[ehvo_hex_indices] = True

    values_ehvo = hexvalues[ehvo_hex_mask]
    values_no_ehvo = hexvalues[~ehvo_hex_mask]

    if len(values_ehvo) == 0 or len(values_no_ehvo) == 0:
        raise ValueError(
            "Both EHVO and non-EHVO hexagons are required "
            "for comparison."
        )

    # Summary statistics
    median_ehvo = np.median(values_ehvo)
    median_no_ehvo = np.median(values_no_ehvo)

    mean_ehvo = np.mean(values_ehvo)
    mean_no_ehvo = np.mean(values_no_ehvo)

    difference = median_ehvo - median_no_ehvo

    percent_difference = (
        difference / median_no_ehvo * 100
        if median_no_ehvo != 0
        else np.nan
    )

    results = {
        "n_ehvo_hex": len(values_ehvo),
        "n_no_ehvo_hex": len(values_no_ehvo),
        "median_ehvo": median_ehvo,
        "median_no_ehvo": median_no_ehvo,
        "mean_ehvo": mean_ehvo,
        "mean_no_ehvo": mean_no_ehvo,
        "median_difference": difference,
        "percent_difference": percent_difference,
    }

    if verbose:
        print("\n------------------------------")
        print("EHVO HEXBIN COMPARISON")
        print("------------------------------")
        print(f"EHVO hexagons:      {len(values_ehvo)}")
        print(f"Non-EHVO hexagons:  {len(values_no_ehvo)}")
        print()
        print(f"Median EHVO:        {median_ehvo:.3f}")
        print(f"Median non-EHVO:    {median_no_ehvo:.3f}")
        print(f"Difference:         {difference:.3f}")
        print(f"Percent difference: {percent_difference:.3f}%")
        print("------------------------------\n")

    return results



if hist_parent == True:
    
    #Histogram
    #MBH EHVO
    c_both_mbh = np.concatenate([parentboth_mbh[bb_both], EHVOboth_mbh[bbboth_EHVO]])
    knuthbinMBH = (stats.knuth_bin_width(c_both_mbh))
    bin_MBH = np.arange(min(c_both_mbh), max(c_both_mbh), knuthbinMBH)
    
    #Lbol EHVO
    c_both_lbol = np.concatenate([parentboth_lbol[bb_both], EHVOboth_lbol[bbboth_EHVO]])
    knuthbinLbol = (stats.knuth_bin_width(c_both_lbol))
    bin_Lbol = np.arange(min(c_both_lbol), max(c_both_lbol), knuthbinLbol)
    
    #Edd EHVO
    c_both_edd = np.concatenate([parentboth_redd[bb_both], EHVOboth_redd[bbboth_EHVO]])
    knuthbinEdd = (stats.knuth_bin_width(c_both_edd))
    bin_Edd = np.arange(min(c_both_edd), max(c_both_edd), knuthbinEdd)
    
    #Histogram Weights
    
    weightsEHVOmbh = [np.ones_like(EHVOboth_mbh[bbboth_EHVO])/float(len(EHVOboth_mbh[bbboth_EHVO]))]
    weightsEHVOlbol = [np.ones_like(EHVOboth_lbol[bbboth_EHVO])/float(len(EHVOboth_lbol[bbboth_EHVO]))]
    weightsEHVOedd = [np.ones_like(EHVOboth_redd[bbboth_EHVO])/float(len(EHVOboth_redd[bbboth_EHVO]))]
    weightsParentmbh = [np.ones_like(parentboth_mbh[bb_both])/float(len(parentboth_mbh[bb_both]))]
    weightsParentlbol = [np.ones_like(parentboth_lbol[bb_both])/float(len(parentboth_lbol[bb_both]))]
    weightsParentedd = [np.ones_like(parentboth_redd[bb_both])/float(len(parentboth_redd[bb_both]))]


elif hist_parent == 'Rankineall':
    #Histogram
    #MBH EHVO
    c_all_mbh = np.concatenate([MBH_rankAll_nonbal,MBH_rankAll_bal, EHVOboth_mbh[bbboth_EHVO]])
    '''
    print("NaN:", np.sum(np.isnan(c_all_mbh)))
    print("Inf:", np.sum(np.isinf(c_all_mbh)))
    print("N:", len(c_all_mbh))
    '''
    #c_all_mbh = c_all_mbh[np.isfinite(c_all_mbh)]
    print("NaN:", np.sum(np.isnan(c_all_mbh)))
    print("Inf:", np.sum(np.isinf(c_all_mbh)))
    print("N:", len(c_all_mbh))
    knuthbinMBH = (stats.knuth_bin_width(c_all_mbh))
    bin_MBH = np.arange(min(c_all_mbh), max(c_all_mbh), knuthbinMBH)
    


    #Lbol EHVO
    c_all_lbol = np.concatenate([LogLbol_rankAll_nonbal,LogLbol_rankAll_bal, EHVOboth_lbol])
    #c_all_lbol = c_all_lbol[np.isfinite(c_all_lbol)]
    print("NaN:", np.sum(np.isnan(c_all_lbol)))
    print("Inf:", np.sum(np.isinf(c_all_lbol)))
    print("N:", len(c_all_lbol))
    knuthbinLbol = (stats.knuth_bin_width(c_all_lbol))
    bin_Lbol = np.arange(min(c_all_lbol), max(c_all_lbol), knuthbinLbol)

    #Edd EHVO
    c_all_edd = np.concatenate([edd_rankAll_nonbal, edd_rankAll_bal, EHVOboth_redd[bbboth_EHVO]])
    #c_all_edd = c_all_edd[np.isfinite(c_all_edd)]
    print("NaN:", np.sum(np.isnan(c_all_edd)))
    print("Inf:", np.sum(np.isinf(c_all_edd)))
    print("N:", len(c_all_edd))
    knuthbinEdd = (stats.knuth_bin_width(c_all_edd))
    bin_Edd = np.arange(min(c_all_edd), max(c_all_edd), knuthbinEdd)

    #Histogram Weights

    weightsEHVOmbh = [np.ones_like(EHVOboth_mbh[bbboth_EHVO])/float(len(EHVOboth_mbh[bbboth_EHVO]))]
    weightsEHVOlbol = [np.ones_like(EHVOboth_lbol[bbboth_EHVO])/float(len(EHVOboth_lbol[bbboth_EHVO]))]
    weightsEHVOedd = [np.ones_like(EHVOboth_redd[bbboth_EHVO])/float(len(EHVOboth_redd[bbboth_EHVO]))]
    
    weights_nonbal_mbh = [np.ones_like(MBH_rankAll_nonbal)/float(len(MBH_rankAll_nonbal))]
    weights_nonbal_lbol = [np.ones_like(LogLbol_rankAll_nonbal[CivBlueCut_nonbal])/float(len(LogLbol_rankAll_nonbal[CivBlueCut_nonbal]))]
    weights_nonbal_edd = [np.ones_like(edd_rankAll_nonbal)/float(len(edd_rankAll_nonbal))]

    weights_bal_mbh = [np.ones_like(MBH_rankAll_bal)/float(len(MBH_rankAll_bal))]
    weights_bal_lbol = [np.ones_like(LogLbol_rankAll_bal[CivBlueCut_bal])/float(len(LogLbol_rankAll_bal[CivBlueCut_bal]))]
    weights_bal_edd = [np.ones_like(edd_rankAll_bal)/float(len(edd_rankAll_bal))]

elif hist_parent == 'parent bal sep':
    #Histogram
    '''
    #MBH EHVO
    c_all_mbh = np.concatenate([parentboth_nonbal_mbh, parentboth_bal_mbh, EHVOboth_mbh[bbboth_EHVO]])
    knuthbinMBH = (stats.knuth_bin_width(c_all_mbh))
    bin_MBH = np.arange(min(c_all_mbh), max(c_all_mbh), knuthbinMBH)
    
    #Lbol EHVO
    c_all_lbol = np.concatenate([parentboth_nonbal_lbol, parentboth_bal_lbol, EHVOboth_lbol[bbboth_EHVO]])
    knuthbinLbol = (stats.knuth_bin_width(c_all_lbol))
    bin_Lbol = np.arange(min(c_all_lbol), max(c_all_lbol), knuthbinLbol)

    #Edd EHVO
    c_all_edd = np.concatenate([parentboth_nonbal_redd, parentboth_bal_redd, EHVOboth_redd[bbboth_EHVO]])
    knuthbinEdd = (stats.knuth_bin_width(c_all_edd))
    bin_Edd = np.arange(min(c_all_edd), max(c_all_edd), knuthbinEdd)
    '''
    
    #Histogram
    #MBH EHVO
    c_all_mbh = np.concatenate([parentboth_mbh[bb_both], EHVOboth_mbh[bbboth_EHVO]])
    knuthbinMBH = (stats.knuth_bin_width(c_all_mbh))
    bin_MBH = np.arange(min(c_all_mbh), max(c_all_mbh), knuthbinMBH)
    
    #Lbol EHVO
    c_all_lbol = np.concatenate([parentboth_lbol[bb_both], EHVOboth_lbol[bbboth_EHVO]])
    knuthbinLbol = (stats.knuth_bin_width(c_all_lbol))
    bin_Lbol = np.arange(min(c_all_lbol), max(c_all_lbol), knuthbinLbol)
    
    #Edd EHVO
    c_all_edd = np.concatenate([parentboth_redd[bb_both], EHVOboth_redd[bbboth_EHVO]])
    knuthbinEdd = (stats.knuth_bin_width(c_all_edd))
    bin_Edd = np.arange(min(c_all_edd), max(c_all_edd), knuthbinEdd)

    #Histogram Weights
    weightsEHVOmbh = [np.ones_like(EHVOboth_mbh[bbboth_EHVO])/float(len(EHVOboth_mbh[bbboth_EHVO]))]
    weightsEHVOlbol = [np.ones_like(EHVOboth_lbol[bbboth_EHVO])/float(len(EHVOboth_lbol[bbboth_EHVO]))]
    weightsEHVOedd = [np.ones_like(EHVOboth_redd[bbboth_EHVO])/float(len(EHVOboth_redd[bbboth_EHVO]))]
    
    weights_nonbal_mbh = [np.ones_like(parentboth_nonbal_mbh)/float(len(parentboth_nonbal_mbh))]
    weights_nonbal_lbol = [np.ones_like(parentboth_nonbal_lbol)/float(len(parentboth_nonbal_lbol))]
    weights_nonbal_edd = [np.ones_like(parentboth_nonbal_redd)/float(len(parentboth_nonbal_redd))]

    weights_bal_mbh = [np.ones_like(parentboth_bal_mbh)/float(len(parentboth_bal_mbh))]
    weights_bal_lbol = [np.ones_like(parentboth_bal_lbol)/float(len(parentboth_bal_lbol))]
    weights_bal_edd = [np.ones_like(parentboth_bal_redd)/float(len(parentboth_bal_redd))]



fig = plt.figure(figsize=(16.5, 10.5), dpi=500)

gs = plt.GridSpec(3, 4)

#LBOL vs MBH
ax1 = fig.add_subplot(gs[1, 0])

#EDD vs MBH
ax2 = fig.add_subplot(gs[2, 0])

#EDD vs LBOL
ax3 = fig.add_subplot(gs[2, 1])

#MBH Histogram
ax4 = fig.add_subplot(gs[0, 0])

#LBOL Histogram
ax5 = fig.add_subplot(gs[1, 1])

#EDD Histogram
ax6 = fig.add_subplot(gs[2, 2])

#Legend for Histograms
ax7 = fig.add_subplot(gs[0,1])

#--------------LBOL vs MBH--------------------
print('----------------LBOL vs MBH--------------------')
cm = ax1.hexbin(parentboth_mbh[bb_both], parentboth_lbol[bb_both],  C = C_values, gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=vmin, vmax=vmax, edgecolor = 'None', zorder= 99)
# Adding stat analysis of hexbins where EHVOs overlap compared to not:
stats_Lbol_MBH = compare_ehvo_hexbins(
    cm,
    EHVOboth_mbh[bbboth_EHVO],
    EHVOboth_lbol[bbboth_EHVO]
)

cbarax = fig.add_axes([-0.03, 0.055, 0.02, 0.93])
cbar = plt.colorbar(cm, cax=cbarax)
cbarax.yaxis.set_ticks_position('left')
cbarax.yaxis.set_label_position('left')

sc = ax1.scatter(EHVOboth_mbh[bbboth_EHVO], EHVOboth_lbol[bbboth_EHVO], c = c_values, zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+2020 and CS+in prep)')
cbarax.set_ylabel(cbar_label, rotation = 90, fontname='serif', fontsize=16)

if plot_bisectors == True:
    above_Lbol_MBH, below_Lbol_MBH = bisector(EHVOboth_mbh[bbboth_EHVO], EHVOboth_lbol[bbboth_EHVO], ax = ax1)
elif plot_bisectors == False:
    above_Lbol_MBH, below_Lbol_MBH = bisector(EHVOboth_mbh[bbboth_EHVO], EHVOboth_lbol[bbboth_EHVO], plot=False)


ax1.set_xlim(8.4,10.5)
ax1.set_ylim(46.1,48.1)
#ax1.legend(loc='upper right')
ax1.set_xticklabels([])
ax1.set_ylabel("log$_{10}$($L_{bol}$/erg s$^{-1}$)", fontname='serif', fontsize=14)
ax1.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#-------------------EDD vs MBH------------------
print('-------------------EDD vs MBH------------------')
cm = ax2.hexbin(parentboth_mbh[bb_both], parentboth_redd[bb_both],  C = C_values, gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=vmin, vmax=vmax, edgecolor = 'None', zorder= 99)

# Adding stat analysis of hexbins where EHVOs overlap compared to not:
stats_REDD_MBH = compare_ehvo_hexbins(
    cm,
    EHVOboth_mbh[bbboth_EHVO],
    EHVOboth_redd[bbboth_EHVO]
)


sc = ax2.scatter(EHVOboth_mbh[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], c = c_values, zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+2020 and CS+in prep)')


if plot_bisectors == True:
    above_REDD_MBH, below_REDD_MBH = bisector(EHVOboth_mbh[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], ax = ax2)
elif plot_bisectors == False:
    above_REDD_MBH, below_REDD_MBH = bisector(EHVOboth_mbh[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], plot=False)

    
ax2.set_xlabel("log$_{10}$($M_{BH}/M_\odot$)", fontname='serif', fontsize=14)
ax2.set_xlim(8.4,10.5)
ax2.set_ylim(-1.5,0.5)
ax2.legend(loc='upper right',prop={'size': 11})
ax2.set_ylabel(r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',fontname='serif', fontsize=14)
ax2.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------EDD vs LBOL---------------------
print('---------------EDD vs LBOL---------------------')
cm = ax3.hexbin(parentboth_lbol[bb_both], parentboth_redd[bb_both],  C = C_values, gridsize = 35, cmap = 'viridis', reduce_C_function = np.median, mincnt = 4, vmin=vmin, vmax=vmax, edgecolor = 'None', zorder= 99)

# Adding stat analysis of hexbins where EHVOs overlap compared to not:
stats_REDD_Lbol = compare_ehvo_hexbins(
    cm,
    EHVOboth_lbol[bbboth_EHVO],
    EHVOboth_redd[bbboth_EHVO]
)

sc = ax3.scatter(EHVOboth_lbol[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], c = c_values, zorder = 100, cmap = 'viridis', edgecolor = 'k', s = 30, vmin=cbar.vmin, vmax=cbar.vmax, label = 'EHVO (RH+2020 and CS+in prep)')

if plot_bisectors == True:
    above_REDD_Lbol, below_REDD_Lbol = bisector(EHVOboth_lbol[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], ax = ax3)
elif plot_bisectors == False:
    above_REDD_Lbol, below_REDD_Lbol = bisector(EHVOboth_lbol[bbboth_EHVO], EHVOboth_redd[bbboth_EHVO], plot=False)

ax3.set_xlabel("log$_{10}$($L_{bol}$/erg s$^{-1}$)", fontname='serif', fontsize=14)
#Removing yticks for second plot
ax3.set_yticklabels([])
ax3.set_xlim(46,48.5)
ax3.set_ylim(-1.5,0.5)
#ax3.legend(loc='upper right')
ax3.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------MBH Histogram---------------------
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.set_xlim(8.4,10.5)

if hist_parent == True:
    ax4.hist(parentboth_mbh[bb_both], bins = bin_MBH, weights= weightsParentmbh, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')
    ax4.hist(EHVOboth_mbh[bbboth_EHVO], bins = bin_MBH, weights = weightsEHVOmbh, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'Rankineall':
    ax4.hist(MBH_rankAll_nonbal, bins = bin_MBH, weights= weights_nonbal_mbh, histtype = 'step', color = 'blue', label = 'Rankine+2020 non-BALs')
    ax4.hist(MBH_rankAll_bal, bins = bin_MBH, weights= weights_bal_mbh, histtype = 'step', color = 'green', label = 'Rankine+2020 BALS')

    ax4.hist(EHVOboth_mbh[bbboth_EHVO], bins = bin_MBH, weights = weightsEHVOmbh, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'parent bal sep':
    ax4.hist(parentboth_nonbal_mbh, bins = bin_MBH, weights= weights_nonbal_mbh, histtype = 'step', color = 'blue')
    ax4.hist(parentboth_bal_mbh, bins = bin_MBH, weights= weights_bal_mbh, histtype = 'step', color = 'green')

    ax4.hist(EHVOboth_mbh[bbboth_EHVO], bins = bin_MBH, weights = weightsEHVOmbh, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')

    
ax4.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)


#---------------LBOL Histogram--------------------
ax5.set_xticklabels([])
ax5.set_yticklabels([])

if hist_parent == True:
    ax5.hist(parentboth_lbol[bb_both], bins = bin_Lbol, weights = weightsParentlbol, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')
    ax5.hist(EHVOboth_lbol[bbboth_EHVO], bins = bin_Lbol, weights = weightsEHVOlbol, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'Rankineall':
    ax5.hist(LogLbol_rankAll_nonbal[CivBlueCut_nonbal], bins = bin_Lbol, weights= weights_nonbal_lbol, histtype = 'step', color = 'blue', label = 'Rankine+2020 non-BALs')
    ax5.hist(LogLbol_rankAll_bal[CivBlueCut_bal], bins = bin_Lbol, weights= weights_bal_lbol, histtype = 'step', color = 'green', label = 'Rankine+2020 BALS')

    ax5.hist(EHVOboth_lbol[bbboth_EHVO], bins = bin_Lbol, weights = weightsEHVOlbol, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'parent bal sep':
    ax5.hist(parentboth_nonbal_lbol, bins = bin_Lbol, weights= weights_nonbal_lbol, histtype = 'step', color = 'blue')
    ax5.hist(parentboth_bal_lbol, bins = bin_Lbol, weights= weights_bal_lbol, histtype = 'step', color = 'green')

    ax5.hist(EHVOboth_lbol[bbboth_EHVO], bins = bin_Lbol, weights = weightsEHVOlbol, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')


ax5.set_xlim(46,48.5)
#ax5.legend(loc='upper right',prop={'size': 6})
# bbox_to_anchor=(0.87, 1.3)
ax5.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#---------------EDD Histogram----------------------
ax6.set_yticklabels([])
ax6.set_xlabel(r'Eddington Ratio $\log_{10}(L_{\mathrm{bol}}/L_{\mathrm{edd}})$',fontname='serif', fontsize=14)

if hist_parent == True:
    ax6.hist(parentboth_redd[bb_both], bins = bin_Edd, weights = weightsParentedd, histtype = 'step', color = 'blue', label = 'Parent RH+ 2020 and RH+ in prep)')
    ax6.hist(EHVOboth_redd[bbboth_EHVO], bins = bin_Edd, weights = weightsEHVOedd, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'Rankineall':
    ax6.hist(edd_rankAll_nonbal, bins = bin_Edd, weights= weights_nonbal_edd, histtype = 'step', color = 'blue', label = 'Rankine+2020 non-BALs')
    ax6.hist(edd_rankAll_bal, bins = bin_Edd, weights= weights_bal_edd, histtype = 'step', color = 'green', label = 'Rankine+2020 BALS')

    ax6.hist(EHVOboth_redd[bbboth_EHVO], bins = bin_Edd, weights = weightsEHVOedd, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')
elif hist_parent == 'parent bal sep':
    ax6.hist(parentboth_nonbal_redd, bins = bin_Edd, weights= weights_nonbal_edd, histtype = 'step', color = 'blue')
    ax6.hist(parentboth_bal_redd, bins = bin_Edd, weights= weights_bal_edd, histtype = 'step', color = 'green')

    ax6.hist(EHVOboth_redd[bbboth_EHVO], bins = bin_Edd, weights = weightsEHVOedd, histtype = 'step', color = 'red', label='EHVO (RH+ 2020 and RH+ in prep)')


ax6.set_xlim(-1.5,0.5)
ax6.tick_params(top=True, labeltop=False, bottom=True, right = True, labelright = False)

#Legend

if hist_parent == True:
    parent_line = plt.Line2D([0], [0], color='blue', label='Parent (RH+2020 and CS+in prep)')
elif hist_parent == 'Rankineall':
    nonbal_line = plt.Line2D([0], [0], color='blue', label='Rankine+2020 non-BALs')
    bal_line = plt.Line2D([0], [0], color='green', label='Rankine+2020 BALs')
elif hist_parent == 'parent bal sep':
    nonbal_line = plt.Line2D([0], [0], color='blue', label='Parent (RH+2020 and CS+in prep) BALs')
    bal_line = plt.Line2D([0], [0], color='green', label='Parent (RH+2020 and CS+in prep) non-BALs')
EHVO_line = plt.Line2D([0], [0], color='red', label='EHVO (RH+2020 and CS+in prep)')
fitline_proxy = plt.Line2D([0], [0], color='k', linestyle='-', label='EHVO Fit Line')
sepline_proxy = plt.Line2D([0], [0], color='red', linestyle='--', label='EHVO Bisector')
ax7.set_xticks([])
ax7.set_yticks([])
if plot_bisectors == True:
    if hist_parent == True:
        ax7.legend(handles=[parent_line, EHVO_line, fitline_proxy, sepline_proxy], loc='lower left', prop={'size': 12})
    elif hist_parent == 'Rankineall' or hist_parent == 'parent bal sep':
        ax7.legend(handles=[nonbal_line, bal_line, EHVO_line, fitline_proxy, sepline_proxy], loc='lower left', prop={'size': 12})
    
elif plot_bisectors == False:
    if hist_parent == True:
        ax7.legend(handles=[parent_line, EHVO_line], loc='lower left', prop={'size': 12})
    elif hist_parent == False:
        ax7.legend(handles=[nonbal_line, bal_line, EHVO_line], loc='lower left', prop={'size': 12})

    
    ax7.spines['top'].set_visible(False)
ax7.spines['right'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
ax7.spines['left'].set_visible(False)
#



#plt.subplots_adjust(wspace=0.16, hspace=0.16)
plt.tight_layout()
plt.show()
#filepath = "C:/Users/annam/Desktop/Research/Emission_PhysicalProperties/Plots/CORNERPLOTSV3.pdf"
#plt.savefig(filepath, dpi=400, bbox_inches='tight')


# 2D K-S test ---------------------------------------------------------------------
print()
print()

from scipy.stats import ks_2samp

CIVBlue_EHVOs = EHVOboth_CivBlue[bbboth_EHVO]

# -- Lbol vs MBH [Top left] ----------------------------------------------
statistic, p_value = ks_2samp(CIVBlue_EHVOs[above_Lbol_MBH], CIVBlue_EHVOs[below_Lbol_MBH])
print(f'Sample Size Above = {len(CIVBlue_EHVOs[above_Lbol_MBH])}, Below = {len(CIVBlue_EHVOs[below_Lbol_MBH])}')
print(f'Lbol vs MBH [Top left] p-value = {p_value} & stat = {statistic}')

# -- EDD vs MBH [Top left] ----------------------------------------------
statistic, p_value = ks_2samp(CIVBlue_EHVOs[above_REDD_MBH], CIVBlue_EHVOs[below_REDD_MBH])
print(f'Sample Size Above = {len(CIVBlue_EHVOs[above_REDD_MBH])}, Below = {len(CIVBlue_EHVOs[below_REDD_MBH])}')
print(f'EDD vs MBH [Lower left] p-value = {p_value} & stat = {statistic}')

# -- EDD vs Lbol [Top left] ----------------------------------------------
statistic, p_value = ks_2samp(CIVBlue_EHVOs[above_REDD_Lbol], CIVBlue_EHVOs[below_REDD_Lbol])
print(f'Sample Size Above = {len(CIVBlue_EHVOs[above_REDD_Lbol])}, Below = {len(CIVBlue_EHVOs[below_REDD_Lbol])}')
print(f'EDD vs Lbol [Lower right] p-value = {p_value} & stat = {statistic}')
print('p-values indicate no statistically significant difference in CIV Blueshift values between samples split by bisector.')


