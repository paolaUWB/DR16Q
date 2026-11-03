#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 14:11:15 2026

@author: lilianaflores
"""

import os
import sys
import pandas as pd


#EHVO info files w/o BI_EHVO
infoRankineEHVO = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR16EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv"
infoRankineEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9EHVO_DR14RankineInfo_wEHVOspeed_vmax_reordered.csv"

infoBIEHVO = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/EHVO_DR16_list_wBIspeed.csv"
infoBIEHVO_DR9 = "https://raw.githubusercontent.com/paolaUWB/DR16Q/Liliana/CROSS_CORRELATION/DR9Q_EHVO_list_amy_wBIspeed.csv"


#Extracting values 
dfRHV = pd.read_csv(infoRankineEHVO, header=0)
dfRHV_9 = pd.read_csv(infoRankineEHVO_DR9, header=0)

dfBI = pd.read_csv(infoBIEHVO, header=0)
dfBI_9 = pd.read_csv(infoBIEHVO_DR9, header=0)

###DR16###
EHVO_SDSSName=dfRHV[dfRHV.columns[0]].to_numpy()

EHVOBI_SDSSName=dfBI[dfBI.columns[1]].to_numpy()

###DR9###
EHVO_SDSSName9=dfRHV_9[dfRHV_9.columns[0]].to_numpy()

EHVOBI_SDSSName=dfBI_9[dfBI_9.columns[2]].to_numpy()


### DR16 ###

# Remove J from BI SDSSName for matching
EHVOBI_SDSSName = dfBI[dfBI.columns[1]].str.replace("J", "", regex=False).to_numpy()

# Reorder BI rows to match RankineEHVO order
dfBI_reordered = (
    dfBI.assign(_match_name=dfBI[dfBI.columns[1]].str.replace("J", "", regex=False))
        .set_index("_match_name")
        .loc[EHVO_SDSSName]
        .reset_index(drop=True)
)


### DR9 ###

# Remove J from BI SDSSName for matching
EHVOBI_SDSSName9 = dfBI_9[dfBI_9.columns[2]].str.replace("J", "", regex=False).to_numpy()

# Reorder BI rows to match RankineEHVO order
dfBI_9_reordered = (
    dfBI_9.assign(_match_name=dfBI_9[dfBI_9.columns[2]].str.replace("J", "", regex=False))
          .set_index("_match_name")
          .loc[EHVO_SDSSName9]
          .reset_index(drop=True)
)


# Save reordered files
dfBI_reordered.to_csv("EHVO_DR16_list_wBIspeed_reordered.csv", index=False)
dfBI_9_reordered.to_csv("DR9Q_EHVO_list_amy_wBIspeed_reordered.csv", index=False)



### DR16 CHECK ###

print("DR16:")
print("RankineEHVO rows:", len(dfRHV))
print("BI rows before:", len(dfBI))
print("BI rows after:", len(dfBI_reordered))

# Compare names after removing J from the BI names
names_16 = dfBI_reordered[dfBI_reordered.columns[1]].str.replace(
    "J", "", regex=False
).to_numpy()

print("All names match in order:", (names_16 == EHVO_SDSSName).all())
print("Number of matching rows:", (names_16 == EHVO_SDSSName).sum())


### DR9 CHECK ###

print("\nDR9:")
print("RankineEHVO rows:", len(dfRHV_9))
print("BI rows before:", len(dfBI_9))
print("BI rows after:", len(dfBI_9_reordered))

names_9 = dfBI_9_reordered[dfBI_9_reordered.columns[2]].str.replace(
    "J", "", regex=False
).to_numpy()

print("All names match in order:", (names_9 == EHVO_SDSSName9).all())
print("Number of matching rows:", (names_9 == EHVO_SDSSName9).sum())
